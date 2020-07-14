import os
import pandas as pd
import snakemake

# Author: Logan Blaine (Pellman Lab)
# Credits: Cheng-Zhong Zhang, Etai Jacob

### Global constants and functions

# TODO these could also be speficied in a "default" config file
DEFAULT_THREADS = int(config["DEFAULT_THREADS"])
MAX_THREADS = int(config["MAX_THREADS"])
JAVA_OPTIONS=f'-Xmx12G -Xms12G -XX:ActiveProcessorCount=4 -Djava.io.tmpdir={config["tmp_dir"]}'
GATK = config["gatk_cmd"] #f'gatk --java-options "{JAVA_OPTIONS}"'
PICARD_OPTIONS = f'--TMP_DIR {config["tmp_dir"]} --MAX_RECORDS_IN_RAM {config["max_records"]} --USE_JDK_DEFLATER --USE_JDK_INFLATER --COMPRESSION_LEVEL 5'
GATK_FILTERS = ("-RF MappingQualityReadFilter --minimum-mapping-quality 30 "
                "-RF OverclippedReadFilter --filter-too-short 50")

# NOTE: requires 'samples' to be defined in the config specified at runtime
sample_sheet = pd.read_table(config['samples'], sep=None, engine='python')
#snakemake.utils.validate(samples, "samples.schema.yaml")
all_samples = set(sample_sheet['sample'])
all_bams = [f'processed_bams/{sample}.bam' for sample in all_samples]

# TODO Remove this if making 'group' an optional field
all_groups = list(set(sample_sheet['group']))

# NOTE SLURM needs log dir to exist or jobs will silently fail
os.makedirs("logs/cluster", exist_ok=True)

###  Functions for finding job input to avoid mixing rules and code
def get_sample_from_prefix(wildcards):
    entries = sample_sheet.query(f'prefix=="{wildcards.sample}"')
    return(entries['sample'].iloc[0])


def get_merged_bams_for_sample(wildcards):
    entries = sample_sheet.query(f'sample=="{wildcards.sample}"')
    return([f'merged_bams/{p}.bam' for p in entries['prefix']])


def get_samples_for_group(wildcards):
    names = sample_sheet.query(f'group=="{wildcards.group}"')
    return([f'processed_bams/{sample}.bam' for sample in set(names['sample'])])


# NOTE Put quick tasks here to avoid wasting cluster queue priority
localrules: filter_structural_variants, refilter_vcf, refilter_vcf_final

# Input-only rules to specify targets for common workflows
rule align:
    input:
        expand("processed_bams/{sample}.bam", sample=all_samples)


rule counts:
    input:
        expand("read_depth/{sample}.counts.tsv", sample=all_samples),
        expand("allelic_depth/{sample}.AD.tsv", sample=all_samples)


rule svs:
    input:
        expand("svaba/{sample}.somatic.joint.sv.vcf", sample=all_samples)
        # expand("svaba/{group}.somatic.sv.counts.csv", group=all_groups)


rule metrics:
    input:
        expand("metrics/{sample}.alignment_summary_metrics",
               sample=all_samples),


rule discordant:
    input:
        expand("discordant_reads/{sample}.discordant.bam", sample=all_samples)


# TODO refactor subworkflows into separate files eventually

### PREPROCESSING STEPS ###

rule bwa_map:
    input:
        config['reference'],
        fq1 = "data/{sample}.unmapped.1.fastq.gz",
        fq2 = "data/{sample}.unmapped.2.fastq.gz"
    output:
        temp("merged_bams/{sample}.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    threads: MAX_THREADS
    params:
        sm = get_sample_from_prefix,
        bwa_threads = lambda wildcards, threads: max(1, threads - 1)
    shell:
        "bwa mem -M -Y -t {params.bwa_threads} {input} 2> {log} "
        " | samtools addreplacerg - "
        " -r ID:{wildcards.sample} -r PU:{wildcards.sample} "
        " -r SM:{params.sm} -r LB:{params.sm} -r PL:illumina -o {output}"

rule mark_duplicates:
    input:
        bam = get_merged_bams_for_sample
    output:
        bam = temp("deduped_bams/{sample}.bam"),
        txt = "metrics/{sample}.dup_metrics.txt"
    params:
        so = "queryname",
        px_dist = 2500,
        bams = lambda wildcards, input: ' '.join([f"-I {b}" for b in input.bam])
    log:
        "logs/gatk/MarkDuplicates/{sample}.log"
    threads: DEFAULT_THREADS
#    group: 'preprocessing'
    shell:
        "{GATK} MarkDuplicates {PICARD_OPTIONS} "
        "--OPTICAL_DUPLICATE_PIXEL_DISTANCE {params.px_dist} "
        "{params.bams} -O {output.bam} "
        "-M {output.txt} -ASO {params.so} 2>{log}"

rule sort_bam:
    input:
        "deduped_bams/{sample}.bam"
    output:
        bam = "processed_bams/{sample}.bam"
    params:
        so = "coordinate"
    log:
        "logs/gatk/SortSam/{sample}.log"
    threads: DEFAULT_THREADS
#    group: 'preprocessing'
    shell:
        "{GATK} SortSam {PICARD_OPTIONS} "
        " -I {input} -O {output.bam} -SO {params.so} "
        "--CREATE_INDEX 2>{log}"


### QUALITY CONTROL
# TODO add a QC check step after metrics to exclude bad samples?

rule collect_metrics:
    input:
        ref = config['reference'],
        bam = "processed_bams/{sample}.bam"
    output:
        "metrics/{sample}.alignment_summary_metrics"
    params:
        "--PROGRAM null"
        "--PROGRAM CollectAlignmentSummaryMetrics",
        "--PROGRAM CollectInsertSizeMetrics",
        "--PROGRAM CollectSequencingArtifactMetrics",
        "--PROGRAM CollectGcBiasMetrics"
    log:
        "logs/gatk/CollectMultipleMetrics/{sample}.log"
    threads: DEFAULT_THREADS
    shell:
        "{GATK} CollectMultipleMetrics {PICARD_OPTIONS} {params} "
        "-I {input.bam} -O metrics/{wildcards.sample} -R {input.ref} 2>{log}"

### COPY NUMBER PIPELINE

rule collect_read_counts:
    input:
        intervals = config['intervals'],
        bam = "processed_bams/{sample}.bam"
    output:
        "read_depth/{sample}.counts.tsv"
    params:
        "--interval-merging-rule OVERLAPPING_ONLY",
        "--format TSV"
    log:
        "logs/gatk/CollectReadCounts/{sample}.log"
    threads: DEFAULT_THREADS
    shell:
        "{GATK} CollectReadCounts -I {input.bam} -L {input.intervals} "
        "{params} {GATK_FILTERS} -O {output} 2>{log}"


rule count_reads_allelic:
    input:
        ref = config['reference'],
        intervals = config['snp_sites'],
        bam = "processed_bams/{sample}.bam"
    output:
        "allelic_depth/{sample}.AD.tsv"
    log:
        "logs/gatk/ASEReadCounter/{sample}.log"
    threads: DEFAULT_THREADS
    shell:
        "{GATK} ASEReadCounter -I {input.bam} -V {input.intervals} "
        "-R {input.ref} -O {output} {GATK_FILTERS} 2>{log}"

#### SV PIPELINE ####

rule call_structural_variants:
    input:
        ref = config['reference'],
        bam = "processed_bams/{sample}.bam"
        # simple = config['simple_repeats'],
        # germline = config['germline_svs']
    threads: DEFAULT_THREADS
    params:
        # bams = lambda wildcards, input: ' '.join([f"-t {b}" for b in input.bam]),
        normal = "-n " + config['normal'],
        flags = "--min-overlap 20 --num-assembly-rounds 1"
    log:
        "svaba/{sample}.log"
    output:
        "svaba/{sample}.svaba.unfiltered.somatic.sv.vcf"
    shell:
        "svaba run -a svaba/{wildcards.sample} -p {threads} "
        "-G {input.ref} -t {input.bam} {params} "
        # "-V {input.germline} -R {input.simple}"

rule filter_structural_variants:
    input:
        "svaba/{sample}.svaba.unfiltered.somatic.sv.vcf"
    threads: 1
    log:
        "logs/bcftools/{sample}.log"
    output:
        "svaba/{sample}.svaba.prefiltered.somatic.sv.vcf"
    params:
        "-i '(SPAN>150000 | SPAN==-1) & AD[0]==0 & SECONDARY==0 & MAPQ>30'"
    shell:
        "bcftools view {input} {params} -o {output} 2>{log}"

## TODO eliminate redundancy in filtering steps
rule recount_somatic_svs:
    output:
        "svaba/{sample}.somatic.sv.counts.csv"
    input:
        "svaba/{sample}.svaba.prefiltered.somatic.sv.vcf",
        "processed_bams/{sample}.bam"
    threads: 1
    script:
        "scripts/recount_svs.py"

rule refilter_vcf:
    input:
        "svaba/{sample}.svaba.prefiltered.somatic.sv.vcf",
        "svaba/{sample}.somatic.sv.counts.csv"
    output:
        "svaba/{sample}.somatic.filtered.sv.vcf"
    threads: 1
    script:
        "scripts/filter_vcf_by_tbl.py"

rule joint_call_somatic_svs:
    output:
        "svaba/{sample}.somatic.joint.sv.counts.csv"
    input:
        "svaba/{sample}.somatic.filtered.sv.vcf",
        all_bams
    threads: DEFAULT_THREADS
    script:
        "scripts/recount_svs.py"

rule refilter_vcf_final:
    input:
        "svaba/{sample}.somatic.filtered.sv.vcf",
        "svaba/{sample}.somatic.joint.sv.counts.csv"
    output:
        "svaba/{sample}.somatic.joint.sv.vcf"
    threads: 1
    script:
        "scripts/filter_vcf_by_tbl.py"

# NOTE Here's a draft of a variant calling rule using HaplotypeCaller
# rule call_short_variants:
#     input:
#         ref = config['reference'],
#         bam = get_samples_for_group
#         # simple = config['simple_repeats'],
#         # germline = config['germline_svs']
#     threads: 8
#     params:
#         bams = lambda wildcards, input: ' '.join([f"-I {b}" for b in input.bam])
#         # normal = "-n " + config['normal'],
#         # flags = "--min-overlap 25"
#     log:
#         "logs/gatk/HaplotypeCaller/{group}.log"
#     output:
#         "short_variants/{group}.vcf.gz"
#     shell:
#         "{GATK} HaplotypeCaller {GATK_FILTERS} "
#         "-R {input.ref} {params} -O {output}"
#         # "-V {input.germline} -R {input.simple}"

# NOTE Legacy preprocessing rules from the GATK best practices workflow
# TODO Add back to pipeline for final production version

# rule fastq_to_ubam:
#     input:
#         fq1 = "data/{sample}.unmapped.1.fastq.gz",
#         fq2 = "data/{sample}.unmapped.2.fastq.gz"
#     output:
#         temp("ubams/{sample}.bam")
#     params:
#         sm = get_sample_from_prefix,
#         platform = "illumina"
#     log:
#         "logs/gatk/FastqToSam/{sample}.log"
#     threads: MAX_THREADS
#     shell:
#         "{GATK} FastqToSam "
#         "-F1 {input.fq1} -F2 {input.fq2} -O {output} "
#         "-SM {params.sm} -LB {params.sm} "
#         "-RG {wildcards.sample} -PU {wildcards.sample} "
#         "-PL {params.platform} 2>{log}"

#
# rule merge_ubam:
#     input:
#         ref = config['reference'],
#         ref_dict = config['reference'].rsplit(".", 1)[0] + ".dict",
#         ubam = "ubams/{sample}.bam",
#         bam = "mapped_reads/{sample}.bam"
#     output:
#         temp("merged_bams/{sample}.bam")
#     log:
#         "logs/gatk/MergeBamAlignment/{sample}.log"
#     params:
#         "-SO unsorted",
#         "-MAX_GAPS -1"
#     threads: MAX_THREADS
#     shell:
#         "{GATK} MergeBamAlignment {PICARD_OPTIONS_RECORDS} {PICARD_OPTIONS} "
#         "-R {input.ref} -O {output} {params} "
#         "-UNMAPPED {input.ubam} -ALIGNED {input.bam} 2>{log}"
