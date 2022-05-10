configfile: "config/config.yaml"

import glob
import os

rule check_unique:
    input: expand(rules.merge_interval_gvcfs.output, sample=SAMPLES)
    output: name_map="output/intervals/sample_name_map",stdout="output/merge_gvcf_metrics/" + CALLSET_NAME + ".check_unique.gvcf.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        
        ls -l output/merge_gvcf/*.g.vcf.gz | awk '{{print $9, $9}}' | awk '{{gsub("/","\t",$1);sub(/\./," ")}}$1'| awk '{{print $3"\t"$5}}' > {output.name_map}
        
        set -euo pipefail
        if [[ $(cut -f 1 {output.name_map} | wc -l) -ne $(cut -f 1 {output.name_map} | sort | uniq | wc -l) ]]
        then
          echo "Samples in the sample_name_map are not unique" 1>&2{output.stdout}
          #exit 1 #snakemake uses bash strict mode, so these exit codes will cause errors
        elif [[ $(cut -f 1 {output.name_map} | wc -l) -lt {SAMPLE_NUM_THRESHOLD} ]]
        then
          printf "%s\n" "There are fewer than {SAMPLE_NUM_THRESHOLD} samples in the sample_name_map" \
          "Having fewer than {SAMPLE_NUM_THRESHOLD} samples means there likely isn't enough data to complete joint calling" 1>&2>{output.stdout}
          #exit 1 #snakemake uses bash strict mode, so these exit codes will cause errors
        else
          echo "All good to go, samples in the sample_name_map are unique, and you have plenty of samples to complete joint calliing" 1>&2>{output.stdout}
        fi
        """

checkpoint split_interval_list:
    input: rules.check_unique.output.name_map
    output: directory("output/intervals/genotyping_interval")
    #params:
    benchmark: "benchmarks/split_interval_list/var_calling_interval_list.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        scatter_count=$(wc -l {input} | awk '{{printf("%.0f", $1*{UNBOUNDED_SCATTER_COUNT_SCALE_FACTOR});}}')
        mkdir -p {output}

        gatk --java-options "-Xms3000m -Xmx3250m" \
          SplitIntervals \
          -L {UNPADDED_INTERVAL_FILE} \
          -O  {output} \
          -scatter $scatter_count \
          -R {REF_FASTA} \
          -mode {SCATTER_MODE} \
          --interval-merging-rule OVERLAPPING_ONLY
        """

rule import_gvcfs:
    input: sample_name_map=rules.check_unique.output.name_map ,intervals="output/intervals/genotyping_interval/{geno_interval}.interval_list"
    #input: lambda wildcards: ['output/variant_calling_scatter_interval/{0}.interval_list'.format(geno_interval) \
#        for interval in {wildcards.geno_interval}]
    output: directory("output/genomicsdb/{geno_interval}")
    #params: sample_name_map=rules.check_unique., interval="output/variant_calling_scatter_interval/{geno_interval}.interval_list"
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
          -L {input.intervals} \
          --sample-name-map {input.sample_name_map} \
          --reader-threads 5 \
          --merge-input-intervals \
          --consolidate

        tar -cf {output}.tar {output}
        
        """
rule genotype_gvcfs:
    input: db=rules.import_gvcfs.output, 
    #output: "output/genotype_gvcfs/"+CALLSET_NAME+"_{geno_interval}.vcf.gz"
    output: "output/genotype_gvcfs/" + CALLSET_NAME + ".{geno_interval}.cohort.vcf.gz"
    params:interval="output/intervals/genotyping_interval/{geno_interval}.interval_list"
    benchmark: "benchmarks/genotype_gvcfs/genotype_gvcfs_cohort_{geno_interval}.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        set -euo pipefail

        #tar -xf output/genomicsdb/{wildcards.geno_interval}
        #WORKSPACE=$(basename output/genomicsdb/{wildcards.geno_interval} .tar)
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
    output:
        filtered="output/interval_vcfs/" + CALLSET_NAME + ".{geno_interval}.cohort.variant_filtered.vcf.gz",
        sites_only="output/interval_vcfs/" + CALLSET_NAME + ".{geno_interval}.cohort.variant_filtered_sites_only.vcf.gz"
    benchmark: "benchmarks/interval_vcfs/{geno_interval}.cohort.hard_filter.benchmark.txt"
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
#   x=expand("output/interval_vcfs/" + CALLSET_NAME + ".{geno_interval}.cohort.variant_filtered_sites_only.vcf.gz",
def aggregate_genotyping_intervals(wildcards):
   checkpoint_output = checkpoints.split_interval_list.get(**wildcards).output[0]
   x=expand(rules.hard_filter_make_sites_only_vcf.output.sites_only,
       geno_interval=glob_wildcards(os.path.join(checkpoint_output,"{geno_interval}.interval_list")).geno_interval)
   print(x)
   return x

rule genotyping_intervals_aggregate:
    input: aggregate_genotyping_intervals
    output: "output/gather_vcfs/" + CALLSET_NAME + ".interval_filtered_sites_only_vcf_list.txt"
    shell: "ls -v {input} > {output}"
        
rule gather_vcfs:
    input: rules.genotyping_intervals_aggregate.output
    output: "output/gather_vcfs/" + CALLSET_NAME + ".cohort.variant_filtered_sites_only.vcf.gz"
    benchmark:"benchmarks/gather_vcfs/cohort.gather_vcf.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3 htslib/1.9

        set -euo pipefail

        #input=$(ls output/interval_vcfs/*.cohort.variant_filtered_sites_only.vcf.gz | sed 's/output/--input output/g')
        input=$(cat {input} | sed 's/output/--input output/g')

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
    output:
        recal="output/gather_vcfs/" + CALLSET_NAME + ".cohort.indels.recal",
        tranches="output/gather_vcfs/" + CALLSET_NAME + ".cohort.indels.tranches"
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
    output:
        recal="output/gather_vcfs/" + CALLSET_NAME + ".cohort.snps.recal",
        tranches="output/gather_vcfs/" + CALLSET_NAME + ".cohort.snps.tranches"
    params:
        tranche=' -tranche '.join(SNP_RECALIBRATION_TRANCHE_VALUES),
        an=' -an '.join(SNP_RECALIBRATION_ANNOTATION_VALUES)
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
    output: "output/recalibrated_vcf_interval/" + CALLSET_NAME + ".cohort.filtered.{geno_interval}.vcf.gz"
    params: "output/recalibrated_vcf_interval/tmp.{geno_interval}.indel.recalibrated.vcf"
    benchmark: "benchmarks/recalibrated_vcf_interval/" + CALLSET_NAME + ".cohort.filtered.{geno_interval}.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        
        set -euo pipefail

        gatk --java-options "-Xms5000m -Xmx6500m" \
          ApplyVQSR \
          -O {params} \
          -V {input.vcf} \
          --recal-file {input.indels_recal} \
          --tranches-file {input.indels_tranches} \
          --truth-sensitivity-filter-level {INDEL_FILTER_LEVEL} \
          --create-output-variant-index true \
          -mode INDEL

        gatk --java-options "-Xms5000m -Xmx6500m" \
          ApplyVQSR \
          -O {output} \
          -V {params} \
          --recal-file {input.snps_recal} \
          --tranches-file {input.snps_tranches} \
          --truth-sensitivity-filter-level {SNP_FILTER_LEVEL} \
          --create-output-variant-index true \
          -mode SNP

        #rm output/recalibrated_vcf_interval/tmp.{wildcards.geno_interval}.indel.recalibrated.vcf
        #rm output/recalibrated_vcf_interval/tmp.{wildcards.geno_interval}.indel.recalibrated.vcf.idx
        rm {params}
	rm {params}.idx
	"""
#   x=expand("output/recalibrated_vcf_interval/" + CALLSET_NAME + ".cohort.filtered.{geno_interval}.vcf.gz",
def aggregate_recal_intervals(wildcards):
   checkpoint_output = checkpoints.split_interval_list.get(**wildcards).output[0]
   x=expand(rules.apply_recal.output,
       geno_interval=glob_wildcards(os.path.join(checkpoint_output,"{geno_interval}.interval_list")).geno_interval)
   print(x)
   return x

rule recal_intervals_aggregate:
    input: aggregate_recal_intervals
    output: "output/recalibrated_vcf_interval/" + CALLSET_NAME + "interval_recalibrated_vcf_list.txt"
    shell: "ls -v {input} > {output}"


# if is small callset ie <= 1000 gvcfs or num_gvcfs <= 1000 then can gather VCF intervals and then collect metrics
# if callset is > 1000 then collect metrics on interval gvcf and gather them later
rule final_gather_vcf:
    #input: expand(rules.apply_recal.output, interval=INTERVALS)
    input: rules.recal_intervals_aggregate.output
    output: "output/recal_final_vcf/" + CALLSET_NAME + ".cohort.vcf.gz"
    benchmark:"output/recal_final_vcf/benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3 htslib/1.9

        #input=$(ls output/recalibrated_vcf_interval/cohort.filtered.*.vcf.gz | sed 's/output/--input output/g')
        input=$(cat {input} | sed 's/output/--input output/g')
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
    output: "output/final_vcf_metrics/" + CALLSET_NAME + ".cohort.variant_calling_detail_metrics","output/final_vcf_metrics/" + CALLSET_NAME + ".cohort.variant_calling_summary_metrics"
    params: prefix="output/final_vcf_metrics/" + CALLSET_NAME + ".cohort",
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
    #input: expand(rules.apply_recal.output, interval=INTERVALS)
    input:
        vcf_list=rules.recal_intervals_aggregate.output,
        name_map=rules.check_unique.output.name_map
    output: "output/final_vcf_metrics/" + CALLSET_NAME + ".cohort.fingerprintcheck"
    #params: sample_name_map="output/sample_name_map"
    benchmark:"benchmarks/final_vcf_metrics/cohort.fingerprintcheck.benchmark.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3 htslib/1.9
        set -eu

        
        num_gvcfs=$(wc -l {input.name_map} | awk '{{print $1}}')
        cpu=$(( $num_gvcfs < 32 ? $num_gvcfs :32))
        memMb=$((1024*$cpu*375/100))
        java=$(($memMb-512))
        java_mem=$(printf "\\"-Xms%sm -Xmx%sm\\"" $java $java)
        echo $java_mem

        
        cat {input.name_map} | awk '{{print $2}}' > output/final_vcf_metrics/gvcf_inputs.list
        cp {input.vcf_list}  output/final_vcf_metrics/vcf_inputs.list

        gatk --java-options "-Xms8000m -Xmx8000m"\
          CrosscheckFingerprints \
          --INPUT output/final_vcf_metrics/gvcf_inputs.list \
          --SECOND_INPUT output/final_vcf_metrics/vcf_inputs.list \
          --HAPLOTYPE_MAP {HAPLOTYPE_DATABASE_FILE} \
          --INPUT_SAMPLE_FILE_MAP {input.name_map} \
          --CROSSCHECK_BY SAMPLE \
          --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
          --NUM_THREADS $cpu \
          --OUTPUT {output}
         """


