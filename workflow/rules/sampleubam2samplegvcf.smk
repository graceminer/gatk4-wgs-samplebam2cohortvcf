configfile: "config/config.yaml"

import glob
import os

# rule create_fastq:
#     input: "output/unmapped_bam/{sample}/{rg}.unmapped.bam"
#     output: "output/unmapped_bam/{sample}/{rg}.interleaved.fastq"
#     benchmark: "benchmarks/create_fq.{sample}.{rg}.txt"
#     #conda: "../envs/test6.yaml"
#     shell:
#         """
# 	module load python/3.8.2  java/1.8.0_211 bwa/0.7.8  picard/2.22.3  R/3.5.3 #python/3.7.3
# 	java -Xms1000m -Xmx1000m -jar $PICARD \
#           SamToFastq \
#           INPUT={input} \
#           FASTQ={output} \
#           INTERLEAVE=true \
#           NON_PF=true
#         """
# rule bwa_mem_alignment:
#     input: rules.create_fastq.output
#     output:
#         bam="output/aligned_bam/{sample}/{rg}.aligned.unsorted.bam",
#         bwalog="output/aligned_bam/{sample}/{rg}.aligned.unsorted.bwa.stderr.log"
#     shell:
#         """
#         module load java/1.8.0_211 bwa/0.7.8 python/3.8.2  picard/2.22.3 R/3.5.3 # python/3.5.0
#         set -o pipefail
#         set -e

#         bwa mem -p -v 3 -t 16 {REF_FASTA} /dev/stdin - 2> >(tee {output.bwalog} >&2) | \
#         java -Dsamjdk.compression_level=2 -Xms1000m -Xmx1000m -jar $PICARD \
#           MergeBamAlignment \
#           VALIDATION_STRINGENCY=SILENT \
#           EXPECTED_ORIENTATIONS=FR \
#           ATTRIBUTES_TO_RETAIN=X0 \
#           ATTRIBUTES_TO_REMOVE=NM \
#           ATTRIBUTES_TO_REMOVE=MD \
#           ALIGNED_BAM=/dev/stdin \
#           UNMAPPED_BAM={input} \
#           OUTPUT={output.bam} \
#           REFERENCE_SEQUENCE={REF_FASTA} \
#           PAIRED_RUN=true \
#           SORT_ORDER="unsorted" \
#           IS_BISULFITE_SEQUENCE=false \
#           ALIGNED_READS_ONLY=false \
#           CLIP_ADAPTERS=false \
#           MAX_RECORDS_IN_RAM=2000000 \
#           ADD_MATE_CIGAR=true \
#           MAX_INSERTIONS_OR_DELETIONS=-1 \
#           PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
#           PROGRAM_RECORD_ID="bwamem" \
#           PROGRAM_GROUP_VERSION={BWA_VERSION} \
#           PROGRAM_GROUP_COMMAND_LINE="bwa mem -p -v 3 -t 16 {REF_FASTA}" \
#           PROGRAM_GROUP_NAME="bwamem" \
#           UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
#           ALIGNER_PROPER_PAIR_FLAGS=true \
#           UNMAP_CONTAMINANT_READS=true \
#           ADD_PG_TAG_TO_READS=false        
#         """
            
rule align_bam_bwa:
    input: "output/unmapped_bam/{sample}/{rg}.unmapped.bam"
    output: bam="output/aligned_bam/{sample}/{rg}.aligned.unsorted.bam",bwalog="output/aligned_bam/{sample}/{rg}.aligned.unsorted.bwa.stderr.log"
    shell:
        """
	module load R/3.5.3 java/1.8.0_211 bwa/0.7.8  python/3.7.3 picard/2.22.3
	set -o pipefail
      set -e
	
	java -Xms1000m -Xmx1000m -jar $PICARD \
          SamToFastq \
          INPUT={input} \
          FASTQ=/dev/stdout \
          INTERLEAVE=true NON_PF=true | \
        bwa mem -p -v 3 -t 16 {REF_FASTA} /dev/stdin - 2> >(tee {output.bwalog} >&2) | \
        java -Dsamjdk.compression_level=2 -Xms1000m -Xmx1000m -jar $PICARD \
          MergeBamAlignment \
          VALIDATION_STRINGENCY=SILENT \
          EXPECTED_ORIENTATIONS=FR \
          ATTRIBUTES_TO_RETAIN=X0 \
          ATTRIBUTES_TO_REMOVE=NM \
          ATTRIBUTES_TO_REMOVE=MD \
          ALIGNED_BAM=/dev/stdin \
          UNMAPPED_BAM={input} \
          OUTPUT={output.bam} \
          REFERENCE_SEQUENCE={REF_FASTA} \
          PAIRED_RUN=true \
          SORT_ORDER="unsorted" \
          IS_BISULFITE_SEQUENCE=false \
          ALIGNED_READS_ONLY=false \
          CLIP_ADAPTERS=false \
          MAX_RECORDS_IN_RAM=2000000 \
          ADD_MATE_CIGAR=true \
          MAX_INSERTIONS_OR_DELETIONS=-1 \
          PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
          PROGRAM_RECORD_ID="bwamem" \
          PROGRAM_GROUP_VERSION={BWA_VERSION} \
          PROGRAM_GROUP_COMMAND_LINE="bwa mem -p -v 3 -t 16 {REF_FASTA}" \
          PROGRAM_GROUP_NAME="bwamem" \
          UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
          ALIGNER_PROPER_PAIR_FLAGS=true \
          UNMAP_CONTAMINANT_READS=true \
          ADD_PG_TAG_TO_READS=false

        java -Xms5000m -jar $PICARD \
          CollectMultipleMetrics \
          INPUT={output.bam} \
          OUTPUT={output.bam} \
          ASSUME_SORTED=true \
          PROGRAM=null \
          PROGRAM=CollectBaseDistributionByCycle \
          PROGRAM=CollectInsertSizeMetrics \
          PROGRAM=MeanQualityByCycle \
          PROGRAM=QualityScoreDistribution \
          METRIC_ACCUMULATION_LEVEL=null \
          METRIC_ACCUMULATION_LEVEL=ALL_READS
	"""

def aggregate_readgroups(wildcards):
   checkpoint_output = checkpoints.revertsam.get(**wildcards).output[0]
   x=expand("output/aligned_bam/{sample}/{rg}.aligned.unsorted.bam",
       sample=wildcards.sample,
       rg=glob_wildcards(os.path.join(checkpoint_output,"{rg}.bam")).rg)
   print(x)
   return x

rule readgroup_aggregate:
    input: aggregate_readgroups
    output: "output/aligned_bam/{sample}.aligned.unsorted.bam.readgroup_list.txt"
    shell: "ls {input} > {output}"    

rule mark_duplicates:
    input: rules.readgroup_aggregate.output
    #output: bam="output/duplicates_marked/{sample}/{sample}.aligned.unsorted.duplicates_marked.bam",metrics="output/duplicates_marked/{sample}/{sample}.aligned.unsorted.duplicates_marked.duplicate_metrics"
    output: bam="output/duplicates_marked/{sample}/{sample}.aligned.unsorted.duplicates_marked.bam",metrics="output/duplicates_marked/{sample}/{sample}.aligned.unsorted.duplicates_marked.duplicate_metrics"
    benchmark: "benchmarks/{sample}.mark_duplicates.txt"
    shell:
        """
        module load java/1.8.0_211 python/3.7.3 picard/2.22.3

        #mkdir -p output/duplicates_marked/{wildcards.sample}/
        #bamlist=$(ls {input}*.bam |  sed 's/output/INPUT=output/g')
        bamlist=$(cat {input} |  sed 's/output/INPUT=output/g')
        echo $bamlist
        
        java -Dsamjdk.compression_level={COMPRESSION_LEVEL} -Xms10g -jar $PICARD \
          MarkDuplicates \
          $bamlist \
          OUTPUT={output.bam} \
          METRICS_FILE={output.metrics} \
          VALIDATION_STRINGENCY=SILENT \
          OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
          ASSUME_SORT_ORDER="queryname" \
          CLEAR_DT="false" \
          ADD_PG_TAG_TO_READS=false
        """

rule sort_bam:
    input: rules.mark_duplicates.output.bam
    output: bam="output/sorted_bam/{sample}/{sample}.aligned.duplicates_marked.sorted.bam",bai="output/sorted_bam/{sample}/{sample}.aligned.duplicates_marked.sorted.bai"
    benchmark: "benchmarks/sorted_bam/{sample}.sortbam.txt"
    shell:
        """
	module load java/1.8.0_211 python/3.7.3 picard/2.22.3

        java -Dsamjdk.compression_level={COMPRESSION_LEVEL} -Xms20g -jar $PICARD \
          SortSam \
          INPUT={input} \
          OUTPUT={output.bam} \
          SORT_ORDER=coordinate \
          CREATE_INDEX=true \
          CREATE_MD5_FILE=true #\
          #MAX_RECORDS_IN_RAM=300000
        """
# rule sort_bam:
#     input: rules.mark_duplicates.output.bam
#     output: "output/sorted_bam/{sample}/{sample}.aligned.duplicates_marked.sorted.bam",
#     params:
#         sort_order="coordinate",
#         extra="CREATE_INDEX=true CREATE_MD5_FILE=true MAX_RECORDS_IN_RAM=300000"
#     log: "logs/picard/sortsam/{sample}.log"
#     wrapper:
#         "v1.1.0/bio/picard/sortsam"

rule check_contamination:
    input: rules.sort_bam.output.bam
    output: "output/sorted_bam/{sample}/{sample}.aligned.duplicates_marked.sorted.preBqsr.selfSM"
    shell:
        """
	module load R/3.5.3 java/1.8.0_211 python/3.7.3 picard/2.22.3 verifyBamID/2014-02-13
        
        # UDPath, MeanPath, and BedPath are functional but depreciated abnd have been replaced with --SVDPrefix where you provide the prefix for all three files
        # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
        # used to read from the selfSM file and calculate contamination, which gets printed out
        
	set -e
        # creates a "output_prefix".selfSM file, a TSV file with 2 rows, 19 columns.
        # need to download and compile the new version of VerifyBamID verifyBamID2 from https://github.com/Griffan/VerifyBamID
        {VerifyBamID} \
        --Verbose \
        --NumPC 4 \
        --Output output/sorted_bam/{wildcards.sample}/{wildcards.sample}.aligned.duplicates_marked.sorted.preBqsr \
        --BamFile {input} \
        --Reference {REF_FASTA} \
        --SVDPrefix {CONTAMINATION_SITES_PREFIX} \
        --DisableSanityCheck \
        
        1>/dev/null
python3 <<CODE
import csv
import sys
with open('{output}') as selfSM:
    reader = csv.DictReader(selfSM, delimiter='\t')
    i = 0
    for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
            # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
            # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
            # vcf and bam.
           sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
           sys.exit(1)
        print(float(row["FREEMIX"])/0.75)
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
            sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
            sys.exit(2)
CODE     
	"""
checkpoint create_sequence_grouping_bqsr:
    input: {REF_DICT}
    #input: "/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict"
    #output: sg="output/intervals/bqsr/sequence_grouping.tsv",sg_unmapped="output/intervals/bqsr/sequence_grouping.unmapped.tsv"
    output: directory("output/intervals/bqsr")
    #params: sg="output/intervals/bqsr/sequence_grouping.tsv",sg_unmapped="output/intervals/bqsr/sequence_grouping.unmapped.tsv"
    shell:
        """
        module load R/3.5.3 java/1.8.0_211 python/3.7.3 picard/2.22.3 verifyBamID/2014-02-13
        mkdir -p {output}
        ### CreateSequenceGroupingTSV (file used to greate sequence groupings uses in BQSR and PrintReads Scatter)
        python3 <<CODE
with open("{REF_DICT}", "r") as ref_dict_file:
    sequence_tuple_list = []
    longest_sequence = 0
    for line in ref_dict_file:
        if line.startswith("@SQ"):
            line_split = line.split("\t")
            # (Sequence_Name, Sequence_Length)
            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
    longest_sequence=sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
# We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
# the last element after a :, so we add this as a sacrificial element.
hg38_protection_tag = ":1+"
# initialize the tsv string with the first sequence
tsv_string = "-L " + sequence_tuple_list[0][0] + hg38_protection_tag
temp_size = sequence_tuple_list[0][1]
for sequence_tuple in sequence_tuple_list[1:]:
    if temp_size + sequence_tuple[1] <= longest_sequence:
        temp_size += sequence_tuple[1]
        tsv_string += "\\t" + "-L " + sequence_tuple[0] + hg38_protection_tag #need to escape backslash in forr snakemake
    else:
        tsv_string += "\\n" + "-L " + sequence_tuple[0] + hg38_protection_tag
        temp_size = sequence_tuple[1]
# add the unmapped sequences as a separate line to ensure that they are recalibrated as well
# with open("sequence_grouping.txt","w") as tsv_file:
with open("{output}/sequence_grouping.txt","w") as tsv_file:
    tsv_file.write(tsv_string)
    tsv_file.close()

tsv_string += '\\n' + "-L " + "unmapped"

# with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
with open("{output}/sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
    tsv_file_with_unmapped.write(tsv_string)
    tsv_file_with_unmapped.close()

# with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
with open("{output}/sequence_grouping_with_unmapped.txt") as f:
    for i, line in enumerate(f):
        with open("output/intervals/bqsr/{{}}.bqsr_interval.list".format(i+1), "w") as g:
            g.write(line)
CODE
           """
#        mkdir -p input/intervals/sequence_grouping_newline/
#        mkdir -p input/intervals/sequence_grouping_unmapped_newline/

#        head -25 {output.sg} | split --lines=1 --numeric-suffixes --additional-suffix=".list" - input/intervals/sequence_grouping_newline/sequence_grouping_
#        sed -n '/chrM:1+/,$p' {output.sg} | sed '1d'  > input/intervals/sequence_grouping_newline/sequence_grouping_25.list

#        head -25 {output.sg_unmapped} | split --lines=1 --numeric-suffixes --additional-suffix=".list" - input/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_
#        sed -n '/chrM:1+/,$p' {output.sg_unmapped} | sed '1d' > input/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_25.list
        

rule base_recalibrator_report:
   input:
       bam=rules.sort_bam.output.bam,
       interval="output/intervals/bqsr/{bqsr_interval}.bqsr_interval.list",
   output: "output/recalibrated_bam/{sample}.recal_data_{bqsr_interval}.csv"
   #params: interval="output/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_{bqsr_interval}.list"
   log:
       "logs/recalibrated_bam/{sample}.base_recal_report.gc_log.{bqsr_interval}.log"
   benchmark:
       "benchmarks/recalibrated_bam/{sample}.bqsr_report.{bqsr_interval}.txt"
   shell:
       """
       module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
       #bam=$(ls {input.bam}*.bam)
       #interval=$(cat {params} | sed 's#^#-L #' | sed 's#$# \\\n#')
       interval=$(cat {input.interval})
       ### BaseRecalibrator
         gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
         -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
         -Xloggc:{log} -Xms5g" \
         BaseRecalibrator \
         -R {REF_FASTA} \
         -I {input.bam} \
         --use-original-qualities \
         -O {output} \
         --known-sites {DBSNP_VCF} \
         --known-sites {KNOWN_SITES_MILLS1000G} \
         --known-sites {KNOWN_SITES_HOMOSAPIENS} \
         $interval
	"""
# def aggregate_readgroups(wildcards):
#    checkpoint_output = checkpoints.revertsam.get(**wildcards).output[0]
#    x=expand("output/aligned_bam/{sample}/{rg}.aligned.unsorted.bam",
#        sample=wildcards.sample,
#        rg=glob_wildcards(os.path.join(checkpoint_output,"{rg}.bam")).rg)
#    print(x)
#    return x

# rule readgroup_aggregate:
#     input: aggregate_readgroups
#     output: "output/aligned_bam/{sample}.aligned.unsorted.bam.readgroup_list.txt"
#     shell: "ls {input} > {output}"    

# rule mark_duplicates:
#     input: rules.readgroup_aggregate.output
           
def aggregate_bqsr_reports(wildcards):
   checkpoint_output = checkpoints.create_sequence_grouping_bqsr.get(**wildcards).output[0]
   x=expand("output/recalibrated_bam/{sample}.recal_data_{bqsr_interval}.csv",
       sample=wildcards.sample,
       bqsr_interval=glob_wildcards(os.path.join(checkpoint_output,"{bqsr_interval}.bqsr_interval.list")).bqsr_interval)
   print(x)
   return x

rule bqsr_intervalreport_aggregate:
    input: aggregate_bqsr_reports
    output: "output/recalibrated_bam/{sample}.bqsr.report_list.txt"
    shell: "ls -v {input}  > {output}"
#     #shell: "ls -v {input} | grep *.list > {output}"

rule gather_bqsr_reports:
    input: rules.bqsr_intervalreport_aggregate.output
    #input:  aggregate_bqsr_reports
    output: "output/recalibrated_bam/{sample}.recal_data.csv"
    #params: "output/recalibrated_bam/{sample}.recal_data_*.csv"
    benchmark:"benchmarks/recalibrated_bam/{sample}.gather_bqsr_reports.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        #input=$(cat {input} | grep {wildcards.sample} |sed 's/output/-I output/g')
        input=$(cat {input} | sed 's/output/-I output/g' )
        #input=$(sed 's/output/-I output/g' {input})
        gatk --java-options "-Xms3000m" \
          GatherBQSRReports \
          $input \
          -O {output}
        """

rule apply_bqsr:
    input:
        bam=rules.sort_bam.output.bam,
        recalibration_report=rules.gather_bqsr_reports.output,
        interval="output/intervals/bqsr/{bqsr_interval}.bqsr_interval.list",
    output:
        bam="output/recalibrated_bam/{sample}.aligned.duplicates_marked.recalibrated_{bqsr_interval}.bam",
        bai="output/recalibrated_bam/{sample}.aligned.duplicates_marked.recalibrated_{bqsr_interval}.bai",
        md5="output/recalibrated_bam/{sample}.aligned.duplicates_marked.recalibrated_{bqsr_interval}.bam.md5",
    params: interval="output/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_{bqsr_interval}.list"
    benchmark: "benchmarks/recalibrated_bam/{sample}.apply_bqsr.{bqsr_interval}.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        #bam=$(ls {input.bam}*.bam)
        #interval=$(cat {params.interval} | sed 's#^#-L #' | sed 's#$# \\\n#')
        interval=$(cat {input.interval})
        gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
          -XX:+PrintGCDetails -Xloggc:{log} \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level={COMPRESSION_LEVEL} -Xms3000m" \
          ApplyBQSR \
          --create-output-bam-md5 \
          --add-output-sam-program-record \
          -R {REF_FASTA} \
          -I {input.bam} \
          --use-original-qualities \
          -O {output.bam} \
          --bqsr-recal-file {input.recalibration_report} \
          $interval
        """

def aggregate_bqsr_bam_intervals(wildcards):
   checkpoint_output = checkpoints.create_sequence_grouping_bqsr.get(**wildcards).output[0]
   x=expand("output/recalibrated_bam/{sample}.aligned.duplicates_marked.recalibrated_{bqsr_interval}.bam",
       sample=wildcards.sample,
       bqsr_interval=glob_wildcards(os.path.join(checkpoint_output,"{bqsr_interval}.bqsr_interval.list")).bqsr_interval)
   print(x)
   return x

rule bam_intervals_aggregate:
    input: aggregate_bqsr_bam_intervals
    output: "output/recalibrated_bam/{sample}.interval_bam_list.txt"
    shell: "ls -v {input} > {output}"
           
rule gather_recal_bams:
     input: rules.bam_intervals_aggregate.output
     #input: expand(rules.apply_bqsr.output, sample=SAMPLE, interval=INTERVALS)
     #input: directory("output/recalibrated_bam/")
     output: "output/final_bam/{sample}.final.bam"
     #params: "output/recalibrated_bam/{sample}.aligned.duplicates_marked.recalibrated_*.bam"
     log: "logs/final_bam/{sample}.final.bam.log"
     benchmark: "benchmarks/final_bam/{sample}.final.bam.txt"
     shell:
         """
         module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

         #test=$(ls output/recalibrated_bam/{wildcards.sample}/*.bam | sed 's/output/INPUT=output/g')
         #input=$(ls {params} | sed 's/output/INPUT=output/g' )
         input=$(cat {input} | sed 's/output/INPUT=output/g' )
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
        seq_cache_populate.pl -root ./output/ref/cache {REF_FASTA}
        export REF_PATH=:
        export REF_CACHE=./output/ref/cache/%2s/%2s/%s

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
# java -Xms1g -jar $PICARD \
#           IntervalListTools \
#           SCATTER_COUNT=50 \
#           SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
#           UNIQUE=true \
#           SORT=true \
#           BREAK_BANDS_AT_MULTIPLES_OF=1000000 \
#           INPUT=/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-wgs_calling_regions.hg38.interval_lis \
#           OUTPUT="test/hc_interval"


checkpoint haplotypecaller_interval:
    input: {CALLING_INTERVAL_LIST}
    #input: "/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-wgs_calling_regions.hg38.interval_list"
    #output: "output/intervals/hc_interval_list/{hc_interval}/scattered.interval_list"
    output: directory("output/intervals/hc_interval_list")
    benchmark: "benchmarks/scatter_interval_list/scatter_interval_list.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        mkdir -p {output}
        set -e
        java -Xms1g -jar $PICARD \
          IntervalListTools \
          SCATTER_COUNT={HAPLOTYPE_SCATTER_COUNT} \
          SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
          UNIQUE=true \
          SORT=true \
          BREAK_BANDS_AT_MULTIPLES_OF={BREAK_BANDS_AT_MULTIPLES_OF} \
          INPUT={input} \
          OUTPUT={output}

        python3 <<CODE
import glob, os
# Works around a JES limitation where multiples files with the same name overwrite each other when globbed
intervals = sorted(glob.glob("{output}/*/scatter.interval_list"))
for i, interval in enumerate(intervals):
  (directory, filename) = os.path.split(interval)
  newName = os.path.join(directory, str(i + 1) + filename)
  os.rename(interval, newName)
print(len(intervals))
CODE
        """

rule gatk4_haplotype_caller:
    input:
        bam=rules.gather_recal_bams.output,  
        interval="output/intervals/hc_interval_list/{hc_interval}/scattered.interval_list"
    output: "output/gvcf/{sample}.{hc_interval}.g.vcf.gz"
    #params: interval_file="output/intervals/scatter_interval_list/*/{hc_interval}scattered.interval_list"        
    benchmark: "benchmarks/gvcf/{sample}.gvcf.{hc_interval}.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        #this also works and uses $interval_file as opposed to params
        #interval_file=$(ls input/scatter_interval_list/{wildcards.hc_interval}/scattered.interval_list)
        set -e
        gatk --java-options "-Xms5g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
          HaplotypeCaller \
          -R {REF_FASTA} \
          -I {input.bam} \
          -L {input.interval} \
          -O {output} \
          -contamination 0.0 \
          -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
          -new-qual \
          -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
           -ERC GVCF 
        """

def aggregate_haplotypecaller_intervals(wildcards):
   checkpoint_output = checkpoints.haplotypecaller_interval.get(**wildcards).output[0]
   x=expand("output/gvcf/{sample}.{hc_interval}.g.vcf.gz",
       sample=wildcards.sample,
       hc_interval=glob_wildcards(os.path.join(checkpoint_output,"{hc_interval}","scattered.interval_list")).hc_interval)
   print(x)
   return x

rule haplotypecaller_intervals_aggregate:
    input: aggregate_haplotypecaller_intervals
    output: "output/gvcf/{sample}.interval.list"
    shell: "ls -v {input} > {output}"
# def aggregate_bqsr_bam_intervals(wildcards):
#    checkpoint_output = checkpoints.create_sequence_grouping_bqsr.get(**wildcards).output[0]
#    x=expand("output/recalibrated_bam/{sample}.aligned.duplicates_marked.recalibrated_{bqsr_interval}.bam",
#        sample=wildcards.sample,
#        bqsr_interval=glob_wildcards(os.path.join(checkpoint_output,"{bqsr_interval}.bqsr_interval.list")).bqsr_interval)
#    print(x)
#    return x

# rule bam_intervals_aggregate:
#     input: aggregate_bqsr_bam_intervals
#     output: "output/recalibrated_bam/{sample}.interval_bam_list.txt"
#     shell: "ls -v {input} > {output}"    
           
rule merge_interval_gvcfs:
    input: rules.haplotypecaller_intervals_aggregate.output
    output: "output/merge_gvcf/{sample}.g.vcf.gz"
    benchmark: "benchmarks/merge_gvcf/{sample}.merge_interval_gvcf.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        #input=$(cat {input} | sed 's/output/INPUT=output/g')
        java -Xms2000m -jar $PICARD \
          MergeVcfs \
          I={input} \
          O={output}
        """
rule validate_gvcf:
    input: rules.merge_interval_gvcfs.output  
    output: "output/merge_gvcf_metrics/{sample}.validate_gvcf.metrics"
    benchmark: "benchmarks/merge_gvcf/{sample}.validate_gvcf.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        
        gatk --java-options -Xms6000m \
          ValidateVariants \
          -V {input} \
          -R {REF_FASTA} \
          -L {CALLING_INTERVAL_LIST} \
          --validate-GVCF \
          --validation-type-to-exclude ALLELES \
          --dbsnp {DBSNP_VCF}  \
        | tee {output}
        """
rule collect_variant_calling_metrics:
    input: rules.merge_interval_gvcfs.output 
    output: "output/merge_gvcf_metrics/{sample}.g.vcf.variant_calling_summary_metrics","output/merge_gvcf_metrics/{sample}.g.vcf.variant_calling_detail_metrics"
    params: "output/merge_gvcf_metrics/{sample}.g.vcf"
    benchmark: "benchmarks/merge_gvcf_metrics/{sample}.variant_calling_summary_metrics.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        java -Xms2000m -jar $PICARD \
          CollectVariantCallingMetrics \
          INPUT={input} \
          OUTPUT={params} \
          DBSNP={DBSNP_VCF} \
          SEQUENCE_DICTIONARY={REF_DICT} \
          TARGET_INTERVALS={EVALUATION_INTERVAL_LIST} \
          GVCF_INPUT=true
        """
