configfile: "config/config.yaml"

#REF_DICT=config["REF_DICT"]

# scatter intervals
#CALLING_INTERVAL_LIST=config["CALLING_INTERVAL_LIST"]
#HAPLOTYPE_SCATTER_COUNT=config["HAPLOTYPE_SCATTER_COUNT"]
#BREAK_BANDS_AT_MULTIPLES_OF=config["BREAK_BANDS_AT_MULTIPLES_OF"]

#SNPEFF = config["SNPEFF"]
#SNPSIFT = config["SNPSIFT"]
#GNOMAD = config["GNOMAD"]
#MAGMA = config["MAGMA"]
#GENOME = config["GENOME"]
#NAME = "1096013095.final.bam"

#CHROMOSOMES = range(1, 23) # ["all"]
#PATH1 = os.getcwd()

SAMPLES, = glob_wildcards("input/{sample}.final.bam")
#(SAMPLES,RGS) = glob_wildcards("output/split_rg/{sample}/{rg}.bam")


wildcard_constraints:
    chr="\d+"
localrules: #bam2ubam

rule main1:
    input:
        #expand("input/{sample}.final.bam", sample=SAMPLES),
        expand("output/split_rg/{sample}/", sample=SAMPLES),
        #expand("output/unmapped_bam/{sample}/{rg}.quality_yield_metrics",zip, sample=SAMPLES,rg=RGS),
        expand("output/unmapped_bam/{sample}.unmapped_bam.readgroup_list.txt", sample=SAMPLES),
        
#         expand("input/intervals/sequence_grouping.tsv"),
#         expand("input/intervals/sequence_grouping.unmapped.tsv"),
#         "input/scatter_interval_list/",
# rule test1:
#     input: "../input/{sample}.final.bam"
#     output: "input/{sample}.input.bam"
#     shell:
#         """ ln -s {input} {output}"""
# rule revertsam:
#     input: "input/{sample}.final.bam"
#     output: directory("output/split_rg/{sample}/")
#     shell:
#         """
#         module load java/1.8.0_211  python/3.7.3 gatk/4.2.0.0
        

#         gatk --java-options "-Xmx1g" \
#   	    RevertSam \
# 	    --INPUT {input} \
# 	    --OUTPUT {output} \
#             --OUTPUT_BY_READGROUP true \
#             --VALIDATION_STRINGENCY LENIENT \
#             --ATTRIBUTE_TO_CLEAR FT \
#             --ATTRIBUTE_TO_CLEAR CO \
#             --SORT_ORDER coordinate
# 	touch {output}
#         module unload java/1.8.0_211  python/3.7.3 singularity/3.2.1
#         """
# #this revert sam rule works, but does not allow exectution of sort sam rule bc the readgroup info is unknown and the glob wildcards cannot access it until the output files of revert sam exist
# rule revert_bam:
#     input:
#         "input/{sample}.final.bam",
#     output:
#         directory("output/split_rg/{sample}/")
#     log:
#         "logs/picard/revert_sam/{sample}.log",
#     params:
#         extra="--OUTPUT_BY_READGROUP true --VALIDATION_STRINGENCY LENIENT --ATTRIBUTE_TO_CLEAR FT --ATTRIBUTE_TO_CLEAR CO --SORT_ORDER coordinate",
#         # optional: Extra arguments for picard.
#     # optional specification of memory usage of the JVM that snakemake will respect with global
#     # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
#     # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
#     # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
#     resources:
#         mem_mb=1024,
#     wrapper:
#         "v1.1.0/bio/picard/revertsam"

# checkpoint revert_bam:
#     input:
#         "input/{sample}.final.bam",
#     output:
#         directory("output/split_rg/{sample}")
#     log:
#         "logs/picard/revert_sam/{sample}.log",
#     params:
#         extra="--OUTPUT_BY_READGROUP true --VALIDATION_STRINGENCY LENIENT --ATTRIBUTE_TO_CLEAR FT --ATTRIBUTE_TO_CLEAR CO --SORT_ORDER coordinate",
#         # optional: Extra arguments for picard.
#     # optional specification of memory usage of the JVM that snakemake will respect with global
#     # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
#     # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
#     # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
#     resources:
#         mem_mb=1024,
#     wrapper:
#         "v1.1.0/bio/picard/revertsam"


checkpoint revertsam:
    input: "input/{sample}.final.bam"
    output: readgroups=directory("output/split_rg/{sample}")
    shell:
        """
        module load java/1.8.0_211  python/3.5.0  gatk/4.2.0.0
        mkdir -p output/split_rg/{wildcards.sample}
        gatk --java-options "-Xmx1g" \
  	    RevertSam \
	    --INPUT {input} \
	    --OUTPUT {output} \
            --OUTPUT_BY_READGROUP true \
            --VALIDATION_STRINGENCY LENIENT \
            --ATTRIBUTE_TO_CLEAR FT \
            --ATTRIBUTE_TO_CLEAR CO \
            --SORT_ORDER coordinate
	#touch {output}
        module unload java/1.8.0_211  python/3.7.3 singularity/3.2.1
        """
# rule create_sequence_grouping_newline:
#     input: "/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict"
#     output: sg="output/intervals/sequence_grouping.tsv",sg_unmapped="output/intervals/sequence_grouping.unmapped.tsv"
#     shell:
#         """
#         module load R/3.5.3 java/1.8.0_211 python/3.7.3 picard/2.22.3 verifyBamID/2014-02-13
#         ### CreateSequenceGroupingTSV (file used to greate sequence groupings uses in BQSR and PrintReads Scatter)
#         python3 <<CODE
# with open("{REF_DICT}", "r") as ref_dict_file:
#     sequence_tuple_list = []
#     longest_sequence = 0
#     for line in ref_dict_file:
#         if line.startswith("@SQ"):
#             line_split = line.split("\t")
#             # (Sequence_Name, Sequence_Length)
#             sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
#     longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
# # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
# # the last element after a :, so we add this as a sacrificial element.
# hg38_protection_tag = ":1+"
# # initialize the tsv string with the first sequence
# tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
# temp_size = sequence_tuple_list[0][1]
# for sequence_tuple in sequence_tuple_list[1:]:
#     if temp_size + sequence_tuple[1] <= longest_sequence:
#         temp_size += sequence_tuple[1]
#         tsv_string += '''\n''' + sequence_tuple[0] + hg38_protection_tag
#     else:
#         tsv_string += '''\n''' + sequence_tuple[0] + hg38_protection_tag
#         temp_size = sequence_tuple[1]
# # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
# # with open("sequence_grouping.txt","w") as tsv_file:
# with open("{output.sg}","w") as tsv_file:
#     tsv_file.write(tsv_string)
#     tsv_file.close()
                                   
# tsv_string += '''\n''' + "unmapped"

# # with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
# with open("{output.sg_unmapped}","w") as tsv_file_with_unmapped:
#     tsv_file_with_unmapped.write(tsv_string)
#     tsv_file_with_unmapped.close()
# CODE
#         mkdir -p input/intervals/sequence_grouping_newline/
#         mkdir -p input/intervals/sequence_grouping_unmapped_newline/
        
#         head -25 {output.sg} | split --lines=1 --numeric-suffixes --additional-suffix=".list" - input/intervals/sequence_grouping_newline/sequence_grouping_
#         sed -n '/chrM:1+/,$p' {output.sg} | sed '1d'  > input/intervals/sequence_grouping_newline/sequence_grouping_25.list
        
#         head -25 {output.sg_unmapped} | split --lines=1 --numeric-suffixes --additional-suffix=".list" - input/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_
#         sed -n '/chrM:1+/,$p' {output.sg_unmapped} | sed '1d' > input/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_25.list

# 	#head -25 input/intervals/sequence_grouping.tsv | split --lines=1 --numeric-suffixes=1 --additional-suffix=".list" - input/intervals/sequence_grouping_newline/sequence_grouping_
#         #sed -n '/chrM:1+/,$p'  input/intervals/sequence_grouping.tsv | sed '1d' > input/intervals/sequence_grouping_newline/sequence_grouping_00.list	

#         #split --lines=1 --numeric-suffixes=1 --additional-suffix=".list" input/intervals/temp_chr25.tsv input/intervals/sequence_grouping_newline/sequence_grouping_
#         #split --lines=1 --numeric-suffixes=1 --additional-suffix=".list" input/intervals/temp_chr25_unmapped.tsv input/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_

#         #split --lines=1 -d --additional-suffix=".list" {output.sg} input/intervals/sequence_grouping_newline/sequence_grouping_ 
#         #split --lines=1 -d --additional-suffix=".list" {output.sg_unmapped}  input/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_
#         """        
# rule scatter_interval_list:
#     #input: {CALLING_INTERVAL_LIST}
#     input: "/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-wgs_calling_regions.hg38.interval_list"
#     output: directory("output/intervals/scatter_interval_list/")
#     benchmark: "benchmarks/scatter_interval_list/scatter_interval_list.txt"
#     shell:
#         """
#         module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

#         set -e
#         #mkdir out
#         java -Xms1g -jar $PICARD \
#           IntervalListTools \
#           SCATTER_COUNT={HAPLOTYPE_SCATTER_COUNT} \
#           SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
#           UNIQUE=true \
#           SORT=true \
#           BREAK_BANDS_AT_MULTIPLES_OF={BREAK_BANDS_AT_MULTIPLES_OF} \
#           INPUT={input} \
#           OUTPUT={output}

#         python3 <<CODE
# import glob, os
# # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
# intervals = sorted(glob.glob("{output}*/*.interval_list"))
# for i, interval in enumerate(intervals):
#   (directory, filename) = os.path.split(interval)
#   newName = os.path.join(directory, str(i + 1) + filename)
#   os.rename(interval, newName)
# print(len(intervals))
# CODE
#         """

# rule sortsam:
#     input: lambda wildcards: [os.path.join('output/','split_rg/',SAMPLES[i], x + '.bam') for i,x in enumerate(RGS) if x == wildcards.rg]
#     output: "output/unmapped_bam/{sample}/{rg}.unmapped.bam"
#     shell:
#         """
#         module load java/1.8.0_211  python/3.7.3  gatk/4.2.0.0

#         gatk --java-options "-Xmx3g" \
#         SortSam \
#             --INPUT {input} \
#             --OUTPUT {output} \
#             --SORT_ORDER queryname \
#             --MAX_RECORDS_IN_RAM 1000000

#         module unload java/1.8.0_211  python/3.7.3 gatk/4.2.0.0
#         """
# def aggregate_revertsam(wildcards):
#     checkpoint_output = checkpoints.revert_bam.get(**wildcards).output
#     file_names = expand("output/split_rg/{sample}/{rg}.bam",RGS = glob_wildcards(os.path.join(checkpoint_output, "{rg}.bam")).rg)

rule sortsam:
    input: "output/split_rg/{sample}/{rg}.bam"
    output: "output/unmapped_bam/{sample}/{rg}.unmapped.bam"
    shell:
        """
        module load java/1.8.0_211  python/3.5.0  gatk/4.2.0.0

        gatk --java-options "-Xmx3g" \
        SortSam \
            --INPUT {input} \
            --OUTPUT {output} \
            --SORT_ORDER queryname \
            --MAX_RECORDS_IN_RAM 1000000

        module unload java/1.8.0_211  python/3.7.3 gatk/4.2.0.0
        """

    
    #file_names = expand("output/split_rg/{sample}/{rg}.bam",RGS = glob_wildcards(os.path.join(checkpoint_output, "{rg}.bam")).rg)


    
rule collect_quality_yield_metrics:
    input: rules.sortsam.output
    output: "output/unmapped_bam/{sample}/{rg}.quality_yield_metrics"
    #output: "output/{sample}/{rg}.unmapped.bam.quality_yield_metrics"
    shell:
        """
        module load java/1.8.0_211  python/3.7.3 picard/2.22.3

        java -Xms2000m -jar $PICARD \
             CollectQualityYieldMetrics \
             INPUT={input} \
             OQ=true  \
             OUTPUT={output}
         module unload java/1.8.0_211  python/3.7.3 picard/2.22.3
        """
# def aggregate_input(wildcards):
#     checkpoint_output = checkpoints.revertsam.get(**wildcards).output[0]
#     x=expand("output/unmapped_bam/{sample}/{rg}.unmapped.bam",
#         sample=wildcards.sample,
#         rg=glob_wildcards(os.path.join(checkpoint_output,"{rg}.bam")).rg)
#     print(x)
#     return x
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.revertsam.get(**wildcards).output[0]
    x=expand("output/unmapped_bam/{sample}/{rg}.quality_yield_metrics",
        sample=wildcards.sample,
        rg=glob_wildcards(os.path.join(checkpoint_output,"{rg}.bam")).rg)
    print(x)
    return x
rule test_aggregate:
    input: aggregate_input
    output: "output/unmapped_bam/{sample}.unmapped_bam.readgroup_list.txt"
    shell: "ls {input} > {output}"
#REDO THE align_bam_bwa rule so that the log files are part of the output (they arent erased on job fail)
#separate the lines of the sam2fq and merge alignment steps to clean up

rule align_bam_bwa:
    input: rules.sortsam.output
    output: bam="output/aligned_bam/{sample}/{rg}.aligned.unsorted.bam",bwalog="output/aligned_bam/{sample}/{rg}.aligned.unsorted.bwa.stderr.log"
    shell:
        """
        module load R/3.5.3 java/1.8.0_211 bwa/0.7.8  python/3.7.3 picard/2.22.3
        echo "modules loaded"
        #mkdir -p ./output/{wildcards.sample}/aligned_bam

        #BWA_VERSION=$(bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')
        #echo $BWA_VERSION

        #set -o pipefail
        #set -e

        #if [-z $BWA_VERSION]; then
        #    exit 1;
        #fi

        java -Xms1000m -Xmx1000m -jar $PICARD SamToFastq INPUT={input} FASTQ=/dev/stdout INTERLEAVE=true NON_PF=true | bwa mem -p -v 3 -t 16 /{REF_FASTA} /dev/stdin - 2> >(tee {output.bwalog} >&2) | java -Dsamjdk.compression_level=2 -Xms1000m -Xmx1000m -jar $PICARD MergeBamAlignment VALIDATION_STRINGENCY=SILENT EXPECTED_ORIENTATIONS=FR ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_REMOVE=NM ATTRIBUTES_TO_REMOVE=MD ALIGNED_BAM=/dev/stdin UNMAPPED_BAM={input} OUTPUT={output.bam} REFERENCE_SEQUENCE={REF_FASTA} PAIRED_RUN=true SORT_ORDER="unsorted" IS_BISULFITE_SEQUENCE=false ALIGNED_READS_ONLY=false CLIP_ADAPTERS=false MAX_RECORDS_IN_RAM=2000000 ADD_MATE_CIGAR=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant PROGRAM_RECORD_ID="bwamem" PROGRAM_GROUP_VERSION={BWA_VERSION} PROGRAM_GROUP_COMMAND_LINE="bwa mem -p -v 3 -t 16 {REF_FASTA}" PROGRAM_GROUP_NAME="bwamem" UNMAPPED_READ_STRATEGY=COPY_TO_TAG ALIGNER_PROPER_PAIR_FLAGS=true UNMAP_CONTAMINANT_READS=true ADD_PG_TAG_TO_READS=false

        java -Xms5000m -jar $PICARD \
          CollectMultipleMetrics \
          INPUT={output.bam} \
          OUTPUT="output/aligned_bam/{wildcards.sample}/{wildcards.rg}.aligned.unsorted" \
          ASSUME_SORTED=true \
          PROGRAM=null \
          PROGRAM=CollectBaseDistributionByCycle \
          PROGRAM=CollectInsertSizeMetrics \
          PROGRAM=MeanQualityByCycle \
          PROGRAM=QualityScoreDistribution \
          METRIC_ACCUMULATION_LEVEL=null \
          METRIC_ACCUMULATION_LEVEL=ALL_READS
        """





   
