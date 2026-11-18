#!/usr/bin/env nextflow
/*
 * nf-viral-consensus
 * Reference-based consensus genomes, variant calls, coverage QC, and a
 * phylogeny from paired-end Illumina sequencing of viral samples.
 * Author: Patricia K. Kamanda
 */

nextflow.enable.dsl = 2

include { FASTP                          } from './modules/fastp'
include { BOWTIE2_BUILD; BOWTIE2_ALIGN   } from './modules/bowtie2'
include { IVAR_TRIM                      } from './modules/ivar_trim'
include { SAMTOOLS_SORT                  } from './modules/samtools_sort'
include { SAMTOOLS_SORT as SORT_TRIMMED } from './modules/samtools_sort'
include { SAMTOOLS_STATS                 } from './modules/samtools_stats'
include { BCFTOOLS_CALL                  } from './modules/bcftools'
include { BCFTOOLS_CONSENSUS             } from './modules/bcftools'
include { ANNOTATE_VARIANTS              } from './modules/annotate'
include { MAFFT; IQTREE                  } from './modules/phylogeny'
include { SUMMARY; MULTIQC               } from './modules/summary'

def helpMessage() {
    log.info """
    nf-viral-consensus

    Usage:
      nextflow run main.nf --input samplesheet.csv -profile docker

    Required:
      --input               CSV with columns sample,fastq_1,fastq_2 (paired-end Illumina reads)

    Reference (defaults to SARS-CoV-2 MN908947.3, bundled in assets/reference):
      --fasta               Reference genome FASTA
      --gff                 Gene annotation (GFF3) used to label variants with their gene
      --primer_bed          Primer BED for amplicon data; primers are soft-clipped with iVar trim

    Thresholds:
      --min_depth           Minimum depth to call a base and a variant (default: ${params.min_depth})
      --min_variant_qual    Minimum variant QUAL (default: ${params.min_variant_qual})
      --min_genome_fraction Fraction of the genome at min_depth for a sample to PASS (default: ${params.min_genome_fraction})
      --consensus_min_af    Minimum variant allele fraction applied to the consensus (default: ${params.consensus_min_af})
      --skip_phylogeny      Skip the alignment and tree (default: ${params.skip_phylogeny})
    """.stripIndent()
}

def locateFastq(String path, sheet_dir) {
    def f = (path.startsWith('/') || path.contains('://')) ? file(path) : sheet_dir.resolve(path)
    if (!path.contains('://') && !f.exists()) {
        error "FASTQ file not found: ${f}"
    }
    return f
}

workflow {

    if (params.help) {
        helpMessage()
        return
    }
    if (!params.input) {
        error "Please provide a samplesheet with --input"
    }

    def sheet_dir = file(params.input).parent
    ch_reads = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2) {
                error "Every samplesheet row needs sample, fastq_1, and fastq_2: ${row}"
            }
            tuple([ id: row.sample ], [ locateFastq(row.fastq_1, sheet_dir), locateFastq(row.fastq_2, sheet_dir) ])
        }

    fasta = file(params.fasta, checkIfExists: true)
    gff   = file(params.gff, checkIfExists: true)

    // 1. Read trimming and quality control
    FASTP(ch_reads)

    // 2. Alignment to the reference
    BOWTIE2_BUILD(fasta)
    BOWTIE2_ALIGN(FASTP.out.reads, BOWTIE2_BUILD.out.index)
    SAMTOOLS_SORT(BOWTIE2_ALIGN.out.sam, '.sorted')

    // 3. Primer soft-clipping for amplicon data
    if (params.primer_bed) {
        IVAR_TRIM(SAMTOOLS_SORT.out.bam, file(params.primer_bed, checkIfExists: true))
        SORT_TRIMMED(IVAR_TRIM.out.bam, '.primertrim.sorted')
        ch_bam = SORT_TRIMMED.out.bam
    } else {
        ch_bam = SAMTOOLS_SORT.out.bam
    }

    // 4. Coverage and alignment statistics
    SAMTOOLS_STATS(ch_bam)

    // 5. Variant calling, filtering, and normalization
    BCFTOOLS_CALL(ch_bam, fasta)

    // 6. Consensus genome with low-coverage positions masked as N
    primers = params.primer_bed ? file(params.primer_bed, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")
    BCFTOOLS_CONSENSUS(BCFTOOLS_CALL.out.vcf.join(SAMTOOLS_STATS.out.depth), fasta, primers)

    // 7. Variant table with gene labels
    ANNOTATE_VARIANTS(BCFTOOLS_CALL.out.vcf, gff, primers)

    // 8. Multiple alignment and maximum-likelihood tree
    if (!params.skip_phylogeny) {
        MAFFT(BCFTOOLS_CONSENSUS.out.fasta.map { meta, f -> f }.collect(), fasta)
        IQTREE(MAFFT.out.alignment)
    }

    // 9. Run-level summary tables
    SUMMARY(
        SAMTOOLS_STATS.out.coverage.map { meta, f -> f }.collect(),
        SAMTOOLS_STATS.out.depth.map { meta, f -> f }.collect(),
        BCFTOOLS_CONSENSUS.out.fasta.map { meta, f -> f }.collect(),
        ANNOTATE_VARIANTS.out.table.map { meta, f -> f }.collect()
    )

    // 10. MultiQC report
    ch_multiqc = FASTP.out.json.map { meta, f -> f }
        .mix(BOWTIE2_ALIGN.out.log.map { meta, f -> f })
        .mix(SAMTOOLS_STATS.out.flagstat.map { meta, f -> f })
        .mix(BCFTOOLS_CALL.out.stats.map { meta, f -> f })
        .collect()
    MULTIQC(ch_multiqc, file("${projectDir}/assets/multiqc_config.yml"))
}
