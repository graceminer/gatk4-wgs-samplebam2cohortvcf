configfile: "config/config.yaml"
import glob
import os

# import gvcf
BATCH_SIZE=50
INTERVALS,= glob_wildcards("output/variant_calling_scatter_interval/{interval}.interval_list")

# genotype_gvcfs
CALLSET_NAME="aaa_cohort_2022"

# hard_filter_make_sites_only_vcf
# ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
    # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69  
EXCESS_HET_THRESHOLD= 54.69

# indel realignment
INDEL_MAX_GAUSSIANS=4

INDEL_RECALIBRATION_ANNOTATION_VALUES = ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"]
INDEL_RECALIBRATION_TRANCHE_VALUES= ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"]
AXIOMPOLY_RESOURCE_VCF= "/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"

# snp realignment
SNP_MAX_GAUSSIANS=6

SNP_RECALIBRATION_ANNOTATION_VALUES=["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"]
SNP_RECALIBRATION_TRANCHE_VALUES=["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0"]
HAPMAP_RESOURCE_VCF="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-hapmap_3.3.hg38.vcf.gz"
ONE_THOUSAND_GENOMES_RESOURCE_VCF="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-1000G_phase1.snps.high_confidence.hg38.vcf.gz"
OMNI_RESOURCE_VCF="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-1000G_omni2.5.hg38.vcf.gz"

# apply recalibration vcf
INDEL_FILTER_LEVEL=99.0
SNP_FILTER_LEVEL=99.7

REF_FASTA = config["REF_FASTA"]
HAPLOTYPE_DATABASE_FILE="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.haplotype_database.txt"


# base_recalibration_report
REF_DICT=config["REF_DICT"]
DBSNP_VCF=config["DBSNP_VCF"]

#known_indels_sites_vcfs
KNOWN_SITES_MILLS1000G="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

EVALUATION_INTERVAL_LIST="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/intervals_wgs_evaluation_regions.hg38.interval_list"


    
localrules: main
#default combinatorial function is product, which can be replaced with other functions in this case "zip()" to create a tuple so that expand doesnt create all combos of sample and rg which ends up swapping sample dirs and rg and producing all possible cobinations as opposed to sampleA/sampleA_rg_bam

rule main6:
    input:
        #expand("output/variant_calling_scatter_interval/{intervals}.interval_list", intervals=INTERVALS),
        #expand("output/{callset_name}_{intervals}.txt", callset_name=CALLSET_NAME, intervals=INTERVALS),
        #expand("output/genomicsdb/{interval}", interval=INTERVALS),
        #expand("output/genotype_gvcfs/{interval}.cohort.vcf.gz", interval=INTERVALS),
        #expand("output/interval_vcfs/{interval}.cohort.variant_filtered.vcf.gz", interval=INTERVALS),
        #expand("output/interval_vcfs/{interval}.cohort.variant_filtered_sites_only.vcf.gz", interval=INTERVALS),
        #expand("output/gather_vcfs/cohort.variant_filtered_sites_only.vcf.gz", interval=INTERVALS),
        #"output/gather_vcfs/cohort.indels.recal","output/gather_vcfs/cohort.indels.tranches",
        #"output/gather_vcfs/cohort.snps.recal","output/gather_vcfs/cohort.snps.tranches",
        #expand("output/recalibrated_vcf_interval/cohort.filtered.{interval}.vcf.gz", interval=INTERVALS),
        "output/final_vcf_metrics/cohort.variant_calling_detail_metrics", "output/final_vcf_metrics/cohort.variant_calling_summary_metrics","output/final_vcf_metrics/cohort.fingerprintcheck"
        
rule import_gvcfs:
    input: "output/variant_calling_scatter_interval/{interval}.interval_list"
    #input: lambda wildcards: ['output/variant_calling_scatter_interval/{0}.interval_list'.format(interval) \
#        for interval in {wildcards.interval}]
    output: directory("output/genomicsdb/{interval}")
    params: sample_name_map="output/sample_name_map", interval="output/variant_calling_scatter_interval/{intervals}.interval_list"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        
        set -euo pipefail

        rm -rf {output}
        #mkdir -p {output}
        # We've seen some GenomicsDB performance regressions related to intervals, so we're going to pretend we only have a single interval
        # using the --merge-input-intervals arg
        # There's no data in between since we didn't run HaplotypeCaller over those loci so we're not wasting any compute

        # The memory setting here is very important and must be several GiB lower
        # than the total memory allocated to the VM because this tool uses
        # a significant amount of non-heap memory for native libraries.
        # Also, testing has shown that the multithreaded reader initialization
        # does not scale well beyond 5 threads, so don't increase beyond that.

        gatk --java-options "-Xms8000m -Xmx25000m" \
          GenomicsDBImport \
          --genomicsdb-workspace-path {output} \
          --batch-size {BATCH_SIZE} \
          -L {params.interval} \
          --sample-name-map {params.sample_name_map} \
          --reader-threads 5 \
          --merge-input-intervals \
          --consolidate

        tar -cf {output}.tar {output}
        
        """
rule genotype_gvcfs:
    input: db=rules.import_gvcfs.output, 
    #output: "output/genotype_gvcfs/"+CALLSET_NAME+"_{intervals}.vcf.gz"
    output: "output/genotype_gvcfs/{interval}.cohort.vcf.gz"
    params:interval="output/variant_calling_scatter_interval/{interval}.interval_list"
    benchmark: "benchmarks/genotype_gvcfs/genotype_gvcfs_cohort_{interval}.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        set -euo pipefail

        #tar -xf output/genomicsdb/{wildcards.interval}
        #WORKSPACE=$(basename output/genomicsdb/{wildcards.interval} .tar)
        #tar -xf {input.db}
        #WORKSPACE=$({input.db})

        gatk --java-options "-Xms8000m -Xmx25000m" \
          GenotypeGVCFs \
          -R {REF_FASTA} \
          -O {output} \
          -D {DBSNP_VCF} \
          -G StandardAnnotation -G AS_StandardAnnotation \
          --only-output-calls-starting-in-intervals \
          -V gendb://{input.db} \
          -L {params.interval} \
          --allow-old-rms-mapping-quality-annotation-data \
          --merge-input-intervals
        """
rule hard_filter_make_sites_only_vcf:
    input: rules.genotype_gvcfs.output
    output: filtered="output/interval_vcfs/{interval}.cohort.variant_filtered.vcf.gz", sites_only="output/interval_vcfs/{interval}.cohort.variant_filtered_sites_only.vcf.gz"
    benchmark: "benchmarks/interval_vcfs/{interval}.cohort.hard_filter.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        set -euo pipefail

        gatk --java-options "-Xms3000m -Xmx3250m" \
          VariantFiltration \
          --filter-expression "ExcessHet > {EXCESS_HET_THRESHOLD}" \
          --filter-name ExcessHet \
          -O {output.filtered} \
          -V {input}

        gatk --java-options "-Xms3000m -Xmx3250m" \
          MakeSitesOnlyVcf \
          -I {output.filtered} \
          -O {output.sites_only}

        """
rule gather_vcfs:
    input: expand(rules.hard_filter_make_sites_only_vcf.output.sites_only, interval=INTERVALS)
    output: "output/gather_vcfs/cohort.variant_filtered_sites_only.vcf.gz"
    benchmark:"benchmarks/gather_vcfs/cohort.gather_vcf.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3 htslib/1.9

        set -euo pipefail

        input=$(ls output/interval_vcfs/*.cohort.variant_filtered_sites_only.vcf.gz | sed 's/output/--input output/g')

# --ignore-safety-checks makes a big performance difference so we include it in our invocation.
        # This argument disables expensive checks that the file headers contain the same set of
        # genotyped samples and that files are in order by position of first record.
        gatk --java-options "-Xms6000m -Xmx6500m" \
          GatherVcfsCloud \
          --ignore-safety-checks \
          --gather-type BLOCK \
          $input \
          --output {output}

        tabix {output}
        """
rule indel_variant_recal:
    input: rules.gather_vcfs.output
    output: recal="output/gather_vcfs/cohort.indels.recal",tranches="output/gather_vcfs/cohort.indels.tranches"
    params: tranche=' -tranche '.join(INDEL_RECALIBRATION_TRANCHE_VALUES), an=' -an '.join(INDEL_RECALIBRATION_ANNOTATION_VALUES)
    benchmark: "benchmarks/gather_vcfs/cohort.indel_recal.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        
        set -euo pipefail

        gatk --java-options "-Xms24000m -Xmx25000m" \
          VariantRecalibrator \
          -V {input} \
          -O {output.recal} \
          --tranches-file {output.tranches} \
          --trust-all-polymorphic \
          -tranche {params.tranche} \
          -an {params.an} \
          -mode INDEL \
          --max-gaussians {INDEL_MAX_GAUSSIANS} \
          -resource:mills,known=false,training=true,truth=true,prior=12 {KNOWN_SITES_MILLS1000G} \
          -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {AXIOMPOLY_RESOURCE_VCF} \
          -resource:dbsnp,known=true,training=false,truth=false,prior=2 {DBSNP_VCF}
        """
        
#rule snp_variant_recal_model_if_numgvcf_greaterthan_snp_recal_thresh20000:
  #model_report_filename = callset_name + ".snps.model.report",
  #model_report = "~{model_report_filename}"
  #model_report_arg = if defined(model_report) then "--input-model $MODEL_REPORT --output-tranches-for-scatter" else ""
  #~{model_report_arg} #added to variant recal step after -mode SNP
rule snp_variant_recal_classic_if_numgvcf_lessthan_equalto_snp_recal_thresh20000:
    input: rules.gather_vcfs.output
    output: recal="output/gather_vcfs/cohort.snps.recal",tranches="output/gather_vcfs/cohort.snps.tranches"
    params: tranche=' -tranche '.join(SNP_RECALIBRATION_TRANCHE_VALUES), an=' -an '.join(SNP_RECALIBRATION_ANNOTATION_VALUES)
    benchmark: "benchmarks/gather_vcfs/cohort.snps_recal.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        set -euo pipefail

        #MODEL_REPORT=$(model_report)

        gatk --java-options "-Xms1000m -Xmx25000m" \
          VariantRecalibrator \
          -V {input} \
          -O {output.recal} \
          --tranches-file {output.tranches} \
          --trust-all-polymorphic \
          -tranche {params.tranche} \
          -an {params.an} \
          -mode SNP \
          --max-gaussians {SNP_MAX_GAUSSIANS} \
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {HAPMAP_RESOURCE_VCF} \
          -resource:omni,known=false,training=true,truth=true,prior=12 {OMNI_RESOURCE_VCF} \
          -resource:1000G,known=false,training=true,truth=false,prior=10 {ONE_THOUSAND_GENOMES_RESOURCE_VCF} \
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {DBSNP_VCF}
        """
# check if larger or small callset, for small callsets, gather vcf shards then collect metrics
rule apply_recal:
    input:
        vcf=rules.hard_filter_make_sites_only_vcf.output.filtered,
        indels_tranches=rules.indel_variant_recal.output.tranches,
        indels_recal=rules.indel_variant_recal.output.recal,
        snps_tranches=rules.snp_variant_recal_classic_if_numgvcf_lessthan_equalto_snp_recal_thresh20000.output.tranches,
        snps_recal=rules.snp_variant_recal_classic_if_numgvcf_lessthan_equalto_snp_recal_thresh20000.output.recal
    output: "output/recalibrated_vcf_interval/cohort.filtered.{interval}.vcf.gz"
    benchmark: "benchmarks/recalibrated_vcf_interval/cohort.filtered.{interval}.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        
        set -euo pipefail

        gatk --java-options "-Xms5000m -Xmx6500m" \
          ApplyVQSR \
          -O output/recalibrated_vcf_interval/tmp.indel.recalibrated.vcf \
          -V {input.vcf} \
          --recal-file {input.indels_recal} \
          --tranches-file {input.indels_tranches} \
          --truth-sensitivity-filter-level {INDEL_FILTER_LEVEL} \
          --create-output-variant-index true \
          -mode INDEL

        gatk --java-options "-Xms5000m -Xmx6500m" \
          ApplyVQSR \
          -O {output} \
          -V output/recalibrated_vcf_interval/tmp.indel.recalibrated.vcf \
          --recal-file {input.snps_recal} \
          --tranches-file {input.snps_tranches} \
          --truth-sensitivity-filter-level {SNP_FILTER_LEVEL} \
          --create-output-variant-index true \
          -mode SNP

        rm output/recalibrated_vcf_interval/tmp.indel.recalibrated.vcf
        rm output/recalibrated_vcf_interval/tmp.indel.recalibrated.vcf.idx
        """
# if is small callset ie <= 1000 gvcfs or num_gvcfs <= 1000 then can gather VCF intervals and then collect metrics
# if callset is > 1000 then collect metrics on interval gvcf and gather them later
rule final_gather_vcf:
    input: expand(rules.apply_recal.output, interval=INTERVALS)
    output: "output/recal_final_vcf/cohort.vcf.gz"
    benchmark:"output/recal_final_vcf/benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3 htslib/1.9

        input=$(ls output/recalibrated_vcf_interval/cohort.filtered.*.vcf.gz | sed 's/output/--input output/g')
        
        set -euo pipefail

        # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
        # This argument disables expensive checks that the file headers contain the same set of
        # genotyped samples and that files are in order by position of first record.
        gatk --java-options "-Xms6000m -Xmx6500m" \
          GatherVcfsCloud \
          --ignore-safety-checks \
          --gather-type BLOCK \
          $input \
          --output {output}

        tabix {output}

        """

rule collect_final_variant_calling_metrics:
    input: rules.final_gather_vcf.output
    output: "output/final_vcf_metrics/cohort.variant_calling_detail_metrics","output/final_vcf_metrics/cohort.variant_calling_summary_metrics"
    params: prefix="output/final_vcf_metrics/cohort",
    benchmark:"output/final_vcf_metrics/final_variant_calling_metrics.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3 htslib/1.9

        set -euo pipefail

        gatk --java-options "-Xms6000m -Xmx7000m" \
        CollectVariantCallingMetrics \
        --INPUT {input} \
        --DBSNP {DBSNP_VCF} \
        --SEQUENCE_DICTIONARY {REF_DICT} \
        --OUTPUT {params.prefix} \
        --THREAD_COUNT 8 \
        --TARGET_INTERVALS {EVALUATION_INTERVAL_LIST}
        """

rule crosscheck_fingerprints_solo:
    input: expand(rules.apply_recal.output, interval=INTERVALS)
    output: "output/final_vcf_metrics/cohort.fingerprintcheck"
    params: sample_name_map="output/sample_name_map"
    benchmark:"benchmarks/final_vcf_metrics/cohort.fingerprintcheck.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3 htslib/1.9
        set -eu

        
        num_gvcfs=$(wc -l {params.sample_name_map} | awk '{{print $1}}')
        cpu=$(( $num_gvcfs < 32 ? $num_gvcfs :32))
        memMb=$((1024*$cpu*375/100))
        java=$(($memMb-512))
        java_mem=$(printf "\\"-Xms%sm -Xmx%sm\\"" $java $java)
        echo $java_mem

        
        cat {params.sample_name_map} | awk '{{print $2}}' > output/final_vcf_metrics/gvcf_inputs.list
        ls output/recalibrated_vcf_interval/*.filtered.*.vcf.gz > output/final_vcf_metrics/vcf_inputs.list

        gatk --java-options "-Xms8000m -Xmx8000m"\
          CrosscheckFingerprints \
          --INPUT output/final_vcf_metrics/gvcf_inputs.list \
          --SECOND_INPUT output/final_vcf_metrics/vcf_inputs.list \
          --HAPLOTYPE_MAP {HAPLOTYPE_DATABASE_FILE} \
          --INPUT_SAMPLE_FILE_MAP {params.sample_name_map} \
          --CROSSCHECK_BY SAMPLE \
          --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
          --NUM_THREADS $cpu \
          --OUTPUT {output}
         """

    
    
#######################################################################################################                                                                                                                            


                                                                                                                                                                                                             

