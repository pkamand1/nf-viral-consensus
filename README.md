# nf-viral-consensus

A Nextflow (DSL2) pipeline that turns paired-end Illumina sequencing of viral samples into the outputs a genomic surveillance program needs: trimmed and quality-checked reads, alignments to a reference genome, filtered and normalized variant calls labelled with the gene in which they fall, consensus genomes in which poorly covered positions are masked rather than guessed, a per-sample pass or fail decision on genome completeness, and a maximum-likelihood phylogeny of the samples that pass. It handles both metagenomic (SISPA or shotgun) and tiled-amplicon data; for amplicon data, primers are soft-clipped with iVar and calls near primer sites are flagged for review. The SARS-CoV-2 reference (MN908947.3) is bundled as the default, and any other viral reference can be supplied with `--fasta` and `--gff`.

## Workflow

| Step | Tool | Output |
| --- | --- | --- |
| 1. Read trimming and QC | fastp | trimmed reads, fastp reports |
| 2. Alignment | Bowtie2 (`--very-sensitive`), samtools sort and index | sorted BAM |
| 3. Primer trimming (amplicon data only) | iVar trim | primer-clipped BAM |
| 4. Coverage | samtools flagstat, coverage, depth | per-position depth, coverage summary |
| 5. Variant calling | bcftools mpileup, call (haploid), norm, fill-tags, filter | normalized, filtered VCF with allele fractions |
| 6. Consensus | bcftools consensus | consensus FASTA, low-coverage mask |
| 7. Variant table | `bin/annotate_variants.py` | variants with depth, allele fraction, gene, and review flags |
| 8. Phylogeny | MAFFT (reference-anchored), IQ-TREE | alignment and tree of samples passing QC |
| 9. Summary and report | `bin/summarize_samples.py`, MultiQC | sample summary table, combined variant table, MultiQC report |

## Consensus rules

A position is masked as `N` when its depth falls below `--min_depth` (default 10). A variant is applied to the consensus only if it passes the call filters (QUAL of at least 20 and depth of at least `--min_depth`) and is carried by at least `--consensus_min_af` of the reads (default 0.75); deletions are written as `-` so that every consensus keeps reference coordinates. For amplicon data, insertions and deletions at or within 10 bp of a primer site are reported in the variant table but not applied, because amplicon junctions are a well-known source of artifactual indels; single-nucleotide variants at these sites are applied, since primer sequence has already been clipped from the reads. A sample passes QC when at least `--min_genome_fraction` of the genome (default 0.90) is covered at `--min_depth`, and only passing samples enter the phylogeny.

## Quick start

```
git clone https://github.com/pkamand1/nf-viral-consensus.git
cd nf-viral-consensus
nextflow run main.nf -profile test,docker            # synthetic truth set, bundled
nextflow run main.nf -profile test_sispa,docker      # real SARS-CoV-2 SISPA reads (downloads about 10 MB)
nextflow run main.nf -profile test_amplicon,docker   # real SARS-CoV-2 ARTIC V3 amplicon reads
```

On your own data:

```
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results
nextflow run main.nf -profile docker --input samplesheet.csv --primer_bed scheme.primer.bed
nextflow run main.nf -profile apptainer,slurm --input samplesheet.csv --fasta ref.fasta --gff ref.gff
```

The samplesheet has the columns `sample,fastq_1,fastq_2`. Paths may be absolute, may be URLs, or may be relative to the folder that contains the samplesheet.

## Validation

### Synthetic truth set (`-profile test`)

`assets/test/make_test_data.py` simulates three samples from the SARS-CoV-2 reference with known variants, at about 60x coverage with a 0.2% read error rate: sample_A carries 12 SNVs and a 6 bp deletion in the spike gene, sample_B carries the same variants plus 4 more (a descendant lineage), and sample_C carries an unrelated set of 10 SNVs and has no reads over positions 10,001 to 10,500. The truth is recorded in `assets/test/truth_variants.tsv`.

| Check | Result |
| --- | --- |
| Truth variants recovered | 40 of 40 |
| False-positive calls | 0 |
| Truth variants correctly present in the consensus genomes | 40 of 40 |
| sample_C dropout region | masked (742 N in the consensus, against about 80 at the genome ends in the other samples) |
| Tree | sample_A and sample_B grouped together with bootstrap support of 100 |

### Real SARS-CoV-2 data from nf-core/test-datasets

The `test_sispa` and `test_amplicon` profiles read small subsets of public SARS-CoV-2 sequencing runs directly from the [nf-core/test-datasets](https://github.com/nf-core/test-datasets/tree/viralrecon) repository (MIT License).

| Sample | Data type | Mean depth | Genome at 10x or more | Consensus N | SNVs | Indels | QC |
| --- | --- | --- | --- | --- | --- | --- | --- |
| SRR11140744 | SISPA | 147.1 | 97.5% | 746 | 1 | 1 | PASS |
| SRR11140746 | SISPA | 104.3 | 97.8% | 663 | 1 | 1 | PASS |
| SRR11140748 | SISPA | 124.6 | 99.3% | 216 | 1 | 1 | PASS |
| SRR11140750 | SISPA | 4.6 | 13.0% | 26,018 | 0 | 0 | FAIL |
| sample1 | ARTIC V3 amplicon | 474.4 | 97.8% | 649 | 6 | 1 | PASS |
| sample2 | ARTIC V3 amplicon | 369.5 | 90.5% | 2,850 | 6 | 2 | PASS |

The three passing SISPA samples share C17373T and a 3 bp deletion at 20,297, and SRR11140750, with less than 5x mean depth, is correctly failed and kept out of the tree. Amplicon sample1 carries the C241T, C3037T, C14408T, and A23403G (spike D614G) substitutions characteristic of early B.1 lineages. Both amplicon samples show insertions at positions 23,812 to 23,813, at the junction between ARTIC amplicons 78 and 79 and inside the primer site of 79_LEFT; these are flagged in the variant table and kept out of the consensus under the rule described above, which is exactly the situation that rule is designed for.

## Output

```
results/
├── qc/fastp/          fastp reports
├── alignment/         sorted (and primer-trimmed) BAM files, Bowtie2 logs, flagstat, coverage, depth
├── variants/          filtered VCFs with allele fractions, bcftools stats, annotated variant tables
├── consensus/         consensus FASTA and low-coverage mask (BED) per sample
├── phylogeny/         reference-anchored alignment, included samples, IQ-TREE tree and report
├── summary/
│   ├── sample_summary.tsv   depth, completeness, consensus N count, variant counts, PASS or FAIL
│   └── all_variants.tsv     every variant with depth, allele fraction, gene, and review flag
├── multiqc/           MultiQC report
└── pipeline_info/     Nextflow execution report, timeline, and trace
```

The `review_flag` column marks calls that deserve a second look before they are reported: `primer_site` or `near_primer` (amplicon data), `indel`, `frameshift_length` for indels whose length is not a multiple of three, and `below_consensus_af` for minority variants that were not applied to the consensus.

## Software environments

| Tool | Version |
| --- | --- |
| fastp | 1.4.0 |
| Bowtie2 | 2.5.5 |
| samtools | 1.24 |
| bcftools | 1.25 |
| iVar | 1.4.4 |
| MAFFT | 7.525 |
| IQ-TREE | 3.1.4 |
| MultiQC | 1.35 |

The `docker`, `singularity`, and `apptainer` profiles use per-tool BioContainers images pinned to these versions; the `conda` profile builds `envs/environment.yml`; and `containers/Dockerfile` and `containers/apptainer.def` build a single image holding every tool, selected with `--single_container`. The `slurm` and `condor` profiles submit tasks to those schedulers. The validation runs above were made with the `conda` profile.

## Reference files

`assets/reference/` holds the SARS-CoV-2 Wuhan-Hu-1 reference genome (GenBank MN908947.3), its gene annotation, and the ARTIC Network V3 primer scheme, as distributed in nf-core/test-datasets.

## Limitations

The variant table labels each variant with its gene but does not translate codons; amino acid consequences require a tool such as SnpEff or Nextclade with a matching database. Lineage assignment (Pangolin or Nextclade) is not included. The pipeline is designed for paired-end Illumina data and a single, unsegmented reference; segmented genomes should be supplied as a multi-sequence FASTA and checked carefully.

## Author

Patricia K. Kamanda. Released under the MIT License.
