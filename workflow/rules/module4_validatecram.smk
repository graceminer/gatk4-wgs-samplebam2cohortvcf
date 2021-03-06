configfile: "config/config.yaml"
import glob
import os

REF_FASTA = config["REF_FASTA"]
COMPRESSION_LEVEL=2   

## cross check fingerprints (on the aligned.duplicates_marked.sorted.bam)
HAPLOTYPE_DATABASE_FILE="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.haplotype_database.txt"

## check contamination

# base_recalibration_report
#REF_DICT=config["REF_DICT"]
DBSNP_VCF=config["DBSNP_VCF"]

#known_indels_sites_vcfs
KNOWN_SITES_MILLS1000G="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_SITES_HOMOSAPIENS="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Homo_sapiens_assembly38.known_indels.vcf.gz"

# bqsr
GENOTYPES_FINGERPRINT_FILE="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/NA12878_NA12878.hg38.reference.fingerprint.vcf"
SAMPLE_NAME="NA12878"
# check pre validation
MAX_DUPLICATION_IN_REASONABLE_SAMPLE = 0.30
MAX_CHIMERISM_IN_REASONABLE_SAMPLE = 0.15

# wgs metrics
WGS_COVERAGE_INTERVAL_LIST="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/intervals_wgs_coverage_regions.hg38.interval_list"
READ_LENGTH=250




(SAMPLE,RG) = glob_wildcards("output/split_rg/{sample}/{rg}.bam")
INTERVALS, = glob_wildcards("output/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_{interval}.list")
    
    
localrules: main

#default combinatorial function is product, which can be replaced with other functions in this case "zip()" to create a tuple so that expand doesnt create all combos of sample and rg which ends up swapping sample dirs and rg and producing all possible cobinations as opposed to sampleA/sampleA_rg_bam

rule main4444:
    input:
        #expand("output/recalibrated_bam/{sample}.recal_data_{interval}.csv", sample=SAMPLE, interval=INTERVALS),
        #expand("output/recalibrated_bam/{sample}.recal_data.csv", sample=SAMPLE),
        #expand("output/recalibrated_bam/{sample}.aligned.duplicates_marked.recalibrated_{interval}.bam", sample=SAMPLE, interval=INTERVALS),
        #expand("output/final_bam/{sample}.final.bam", sample=SAMPLE),
        expand("output/final_bam/{sample}.final_bam.collect_aggregated_metrics.{metrics}", sample=SAMPLE, metrics=["alignment_summary_metrics","bait_bias_detail_metrics","bait_bias_summary_metrics","gc_bias.detail_metrics","gc_bias.pdf","gc_bias.summary_metrics","insert_size_histogram.pdf","insert_size_metrics","pre_adapter_detail_metrics", "pre_adapter_summary_metrics","quality_distribution.pdf","quality_distribution_metrics", "error_summary_metrics"]),
        expand("output/final_bam/{sample}.final_bam.collect_rg_metrics.{metrics}", sample=SAMPLE, metrics=["alignment_summary_metrics","gc_bias.detail_metrics","gc_bias.pdf","gc_bias.summary_metrics"]),
        expand("output/final_bam/{sample}.final_bam.{ext}", sample=SAMPLE,ext=["duplication.csv","chimerism.csv","fingerprinting_detail_metrics","fingerprinting_summary_metrics","read_group_md5","wgs_metrics","raw_wgs_metrics"]),
        expand("output/final_cram/{sample}.final_cram.{ext}", sample=SAMPLE, ext=["validation_report"]),

        
# note, intervals must have unmapped interval at the end and sorted by chromosome otherwise the gather steps error out
rule base_recalibrator_report:
   input:
       bam=directory("output/sorted_bam/{sample}/"),
   output: "output/recalibrated_bam/{sample}.recal_data_{intervals}.csv"
   params:
       interval="output/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_{intervals}.list"
   log:
       "logs/recalibrated_bam/{sample}.base_recal_report.gc_log.{intervals}.log"
   benchmark:
       "benchmarks/recalibrated_bam/{sample}.bqsr_report.{intervals}.txt"
   shell:
       """
       module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
       bam=$(ls {input.bam}*.bam)
       interval=$(cat {params.interval} | sed 's#^#-L #' | sed 's#$# \\\n#')
       ### BaseRecalibrator
       gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
         -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
         -Xloggc:{log} -Xms5g" \
         BaseRecalibrator \
         -R {REF_FASTA} \
         -I $bam \
         --use-original-qualities \
         -O {output} \
         --known-sites {DBSNP_VCF} \
         --known-sites {KNOWN_SITES_MILLS1000G} \
         --known-sites {KNOWN_SITES_HOMOSAPIENS} \
         $interval
	"""
       
rule gather_bqsr_reports:
    input: expand(rules.base_recalibrator_report.output, sample=SAMPLE, intervals=INTERVALS)
    output: "output/recalibrated_bam/{sample}.recal_data.csv"
    params: "output/recalibrated_bam/{sample}.recal_data_*.csv"
    benchmark:"benchmarks/recalibrated_bam/{sample}.gather_bqsr_reports.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        #input=$(cat {input} | grep {wildcards.sample} |sed 's/output/-I output/g')
        input=$(ls {params} | sed 's/output/-I output/g' )
        gatk --java-options "-Xms3000m" \
          GatherBQSRReports \
          $input \
          -O {output}
        """
        
rule apply_bqsr:
    input: bam=directory("output/sorted_bam/{sample}/"),recalibration_report=rules.gather_bqsr_reports.output
    output: "output/recalibrated_bam/{sample}.aligned.duplicates_marked.recalibrated_{interval}.bam"
    params: interval="output/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_{interval}.list"        
    benchmark: "benchmarks/recalibrated_bam/{sample}.apply_bqsr.{interval}.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        bam=$(ls {input.bam}*.bam)
        interval=$(cat {params.interval} | sed 's#^#-L #' | sed 's#$# \\\n#')
        gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
          -XX:+PrintGCDetails -Xloggc:{log} \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level={COMPRESSION_LEVEL} -Xms3000m" \
          ApplyBQSR \
          --create-output-bam-md5 \
          --add-output-sam-program-record \
          -R {REF_FASTA} \
          -I $bam \
          --use-original-qualities \
          -O {output} \
          --bqsr-recal-file {input.recalibration_report} \
          $interval
        """
 
rule gather_recal_bams:
     input: expand(rules.apply_bqsr.output, sample=SAMPLE, interval=INTERVALS)
     #input: directory("output/recalibrated_bam/")
     output: "output/final_bam/{sample}.final.bam"
     params: "output/recalibrated_bam/{sample}.aligned.duplicates_marked.recalibrated_*.bam"
     log: "logs/final_bam/{sample}.final.bam.log"
     benchmark: "benchmarks/final_bam/{sample}.final.bam.txt"
     shell:
         """
         module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

         #test=$(ls output/recalibrated_bam/{wildcards.sample}/*.bam | sed 's/output/INPUT=output/g')
         input=$(ls {params} | sed 's/output/INPUT=output/g' )
        
         java -Dsamjdk.compression_level={COMPRESSION_LEVEL} -Xms2000m -jar $PICARD \
           GatherBamFiles \
           $input \
           OUTPUT={output} \
           CREATE_INDEX=true \
           CREATE_MD5_FILE=true
         """

rule collect_aggregated_metrics:
    input: rules.gather_recal_bams.output           
    output:
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.alignment_summary_metrics",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.bait_bias_detail_metrics",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.bait_bias_summary_metrics",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.gc_bias.detail_metrics",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.gc_bias.pdf",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.gc_bias.summary_metrics",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.insert_size_histogram.pdf",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.insert_size_metrics",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.pre_adapter_detail_metrics",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.pre_adapter_summary_metrics",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.quality_distribution.pdf",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.quality_distribution_metrics",
        "output/final_bam/{sample}.final_bam.collect_aggregated_metrics.error_summary_metrics"
    params:"output/final_bam/{sample}.final_bam.collect_aggregated_metrics"
    benchmark: "benchmarks/final_bam/{sample}.collect_aggregated_metrics.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        java -Xms5000m -jar $PICARD \
          CollectMultipleMetrics \
          INPUT={input} \
          REFERENCE_SEQUENCE={REF_FASTA} \
          OUTPUT={params} \
          ASSUME_SORTED=true \
          PROGRAM=null \
          PROGRAM=CollectAlignmentSummaryMetrics \
          PROGRAM=CollectInsertSizeMetrics \
          PROGRAM=CollectSequencingArtifactMetrics \
          PROGRAM=QualityScoreDistribution \
          PROGRAM=CollectGcBiasMetrics \
          METRIC_ACCUMULATION_LEVEL=null \
          METRIC_ACCUMULATION_LEVEL=SAMPLE \
          METRIC_ACCUMULATION_LEVEL=LIBRARY
        """
rule check_prevalidation:
    #input: duplication_metrics="output/duplicates_marked/{sample}/{sample}.aligned.unsorted.duplicates_marked.duplicate_metrics",chimerism_metrics=rules.collect_aggregated_metrics.output[0]
    input: chimerism_metrics=rules.collect_aggregated_metrics.output[0]  
    output: duplication="output/final_bam/{sample}.final_bam.duplication.csv", chimerism="output/final_bam/{sample}.final_bam.chimerism.csv"
    benchmark: "benchmarks/final_bam/{sample}.check_prevalidation.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        set -o pipefail
        set -e

        grep -A 1 PERCENT_DUPLICATION output/duplicates_marked/{wildcards.sample}/{wildcards.sample}.aligned.unsorted.duplicates_marked.duplicate_metrics > {output.duplication}
        grep -A 3 PCT_CHIMERAS {input.chimerism_metrics} | grep -v OF_PAIR > {output.chimerism}

        python <<CODE

import csv
with open('{output.duplication}') as dupfile:
  reader = csv.DictReader(dupfile, delimiter='\t')
  for row in reader:
    with open("output/final_bam/{wildcards.sample}.duplication_value.txt","w") as file:
      file.write(row['PERCENT_DUPLICATION'])
      file.close()

with open('{output.chimerism}') as chimfile:
  reader = csv.DictReader(chimfile, delimiter='\t')
  for row in reader:
    with open("output/final_bam/{wildcards.sample}.chimerism_value.txt","w") as file:
      file.write(row['PCT_CHIMERAS'])
      file.close()
CODE

        """
#Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
rule collect_rg_bam_metrics:
    input: rules.gather_recal_bams.output  
    output: "output/final_bam/{sample}.final_bam.collect_rg_metrics.alignment_summary_metrics", "output/final_bam/{sample}.final_bam.collect_rg_metrics.gc_bias.detail_metrics", "output/final_bam/{sample}.final_bam.collect_rg_metrics.gc_bias.pdf", "output/final_bam/{sample}.final_bam.collect_rg_metrics.gc_bias.summary_metrics"
    params: "output/final_bam/{sample}.final_bam.collect_rg_metrics"
    benchmark: "benchmarks/final_bam/{sample}.collect_rg_metrics.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        java -Xms5000m -jar $PICARD \
          CollectMultipleMetrics \
          INPUT={input} \
          REFERENCE_SEQUENCE={REF_FASTA} \
          OUTPUT={params} \
          ASSUME_SORTED=true \
          PROGRAM=null \
          PROGRAM=CollectAlignmentSummaryMetrics \
          PROGRAM=CollectGcBiasMetrics \
          METRIC_ACCUMULATION_LEVEL=null \
          METRIC_ACCUMULATION_LEVEL=READ_GROUP
        """

rule check_fingerprint:
    input: rules.gather_recal_bams.output  
    output: detail="output/final_bam/{sample}.final_bam.fingerprinting_detail_metrics",summary="output/final_bam/{sample}.final_bam.fingerprinting_summary_metrics"
    benchmark: "benchmarks/final_bam/{sample}.check_fingerprint.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        java -Dsamjdk.buffer_size=131072 \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2g  \
          -jar $PICARD \
          CheckFingerprint \
          INPUT={input} \
          SUMMARY_OUTPUT={output.summary} \
          DETAIL_OUTPUT={output.detail} \
          GENOTYPES={GENOTYPES_FINGERPRINT_FILE} \
          HAPLOTYPE_MAP={HAPLOTYPE_DATABASE_FILE} \
          SAMPLE_ALIAS="{SAMPLE_NAME}" \
          IGNORE_READ_GROUPS=true
        """
rule calc_rg_checksum:
    input: rules.gather_recal_bams.output  
    output: "output/final_bam/{sample}.final_bam.read_group_md5"
    benchmark: "benchmarks/final_bam/{sample}.calculate_rg_checksum.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        java -Xms1000m -jar $PICARD \
            CalculateReadGroupChecksum \
            INPUT={input} \
            OUTPUT={output}
        """
# comment thread for the use fast algorith : https://github.com/broadinstitute/picard/issues/1402
rule collect_wgs_metrics:
    input: rules.gather_recal_bams.output 
    output: "output/final_bam/{sample}.final_bam.wgs_metrics"
    benchmark: "benchmarks/final_bam/{sample}.final_bam.wgs_metrics.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        # "USE_FAST_ALGORITHM=true" errors out so I changed it to false.  comment thread onthis error suggest that the quality is the same, it just takes longer. 
        java -Xms2000m -jar $PICARD \
          CollectWgsMetrics \
          INPUT={input} \
          VALIDATION_STRINGENCY=SILENT \
          REFERENCE_SEQUENCE={REF_FASTA} \
          INCLUDE_BQ_HISTOGRAM=true \
          INTERVALS={WGS_COVERAGE_INTERVAL_LIST} \
          OUTPUT={output} \
          USE_FAST_ALGORITHM=false \
          READ_LENGTH={READ_LENGTH}
        """

rule collect_raw_wgs_metrics:
    input: rules.gather_recal_bams.output  
    output: "output/final_bam/{sample}.final_bam.raw_wgs_metrics"
    benchmark: "benchmarks/final_bam/{sample}.final_bam.raw_wgs_metrics.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        # "USE_FAST_ALGORITHM=true" errors out so I changed it to false.  comment thread onthis error suggest that the quality is the same, it just takes longer. 
        java -Xms2000m -jar $PICARD \
          CollectRawWgsMetrics \
          INPUT={input} \
          VALIDATION_STRINGENCY=SILENT \
          REFERENCE_SEQUENCE={REF_FASTA} \
          INCLUDE_BQ_HISTOGRAM=true \
          INTERVALS={WGS_COVERAGE_INTERVAL_LIST} \
          OUTPUT={output} \
          USE_FAST_ALGORITHM=false \
          READ_LENGTH={READ_LENGTH}
        """
rule convert_to_cram:
    input: rules.gather_recal_bams.output  
    output: cram="output/final_cram/{sample}.final.cram",md5="output/final_cram/{sample}.final.cram.md5"
    benchmark: "benchmarks/final_cram/{sample}.final.cram.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3 samtools/1.9
        set -e
        set -o pipefail

        samtools view -C -T {REF_FASTA} {input} | \
        tee {output.cram} | \
        md5sum | awk '{{print $1}}' > {output.md5}

        # Create REF_CACHE. Used when indexing a CRAM
        seq_cache_populate.pl -root ./ref/cache {REF_FASTA}
        export REF_PATH=:
        export REF_CACHE=./ref/cache/%2s/%2s/%s

        samtools index {output.cram}
        
        """

rule validate_cram:
    input: rules.convert_to_cram.output.cram  
    output: "output/final_cram/{sample}.final_cram.validation_report"
    benchmark: "benchmarks/final_cram/{sample}.validate_cram.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        java -Xms2g -jar $PICARD \
          ValidateSamFile \
          INPUT={input} \
          OUTPUT={output} \
          REFERENCE_SEQUENCE={REF_FASTA} \
          MAX_OUTPUT=1000000000 \
          IGNORE=MISSING_TAG_NM \
          MODE=VERBOSE \
          SKIP_MATE_VALIDATION=false \
          IS_BISULFITE_SEQUENCED=false
        """




######################################################################################################                                                                                                                                                                                                                                                                                                                                         
