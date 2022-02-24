configfile: "config/config.yaml"

#REF_DICT=config["REF_DICT"]

# scatter intervals
CALLING_INTERVAL_LIST=config["CALLING_INTERVAL_LIST"]
#HAPLOTYPE_SCATTER_COUNT=config["HAPLOTYPE_SCATTER_COUNT"]
#BREAK_BANDS_AT_MULTIPLES_OF=config["BREAK_BANDS_AT_MULTIPLES_OF"]



SAMPLES, = glob_wildcards("input/{samples}.final.bam")


wildcard_constraints:
    chr="\d+"
localrules: #bam2ubam

# rule main1:
#     input:
#         expand("input/{sample}.final.bam", sample=SAMPLES),
# 	expand("output/split_rg/{sample}/", sample=SAMPLES),
#         expand("input/intervals/sequence_grouping.tsv"),
#         expand("input/intervals/sequence_grouping.unmapped.tsv"),
#         "input/scatter_interval_list/",


rule revertsam:
    input: "input/{sample}.final.bam"
    output: directory("output/split_rg/{sample}/")
    shell:
        """
        module load java/1.8.0_211  python/3.7.3 gatk/4.2.0.0
        

        gatk --java-options "-Xmx1g" \
  	    RevertSam \
	    --INPUT {input} \
	    --OUTPUT {output} \
            --OUTPUT_BY_READGROUP true \
            --VALIDATION_STRINGENCY LENIENT \
            --ATTRIBUTE_TO_CLEAR FT \
            --ATTRIBUTE_TO_CLEAR CO \
            --SORT_ORDER coordinate
	touch {output}
        module unload java/1.8.0_211  python/3.7.3 singularity/3.2.1
        """
rule create_sequence_grouping_newline:
    input: "/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict"
    output: sg="output/intervals/sequence_grouping.tsv",sg_unmapped="output/intervals/sequence_grouping.unmapped.tsv"
    shell:
        """
        module load R/3.5.3 java/1.8.0_211 python/3.7.3 picard/2.22.3 verifyBamID/2014-02-13
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
    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
# We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
# the last element after a :, so we add this as a sacrificial element.
hg38_protection_tag = ":1+"
# initialize the tsv string with the first sequence
tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
temp_size = sequence_tuple_list[0][1]
for sequence_tuple in sequence_tuple_list[1:]:
    if temp_size + sequence_tuple[1] <= longest_sequence:
        temp_size += sequence_tuple[1]
        tsv_string += '''\n''' + sequence_tuple[0] + hg38_protection_tag
    else:
        tsv_string += '''\n''' + sequence_tuple[0] + hg38_protection_tag
        temp_size = sequence_tuple[1]
# add the unmapped sequences as a separate line to ensure that they are recalibrated as well
# with open("sequence_grouping.txt","w") as tsv_file:
with open("{output.sg}","w") as tsv_file:
    tsv_file.write(tsv_string)
    tsv_file.close()
                                   
tsv_string += '''\n''' + "unmapped"

# with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
with open("{output.sg_unmapped}","w") as tsv_file_with_unmapped:
    tsv_file_with_unmapped.write(tsv_string)
    tsv_file_with_unmapped.close()
CODE
        mkdir -p input/intervals/sequence_grouping_newline/
        mkdir -p input/intervals/sequence_grouping_unmapped_newline/
        
        head -25 {output.sg} | split --lines=1 --numeric-suffixes --additional-suffix=".list" - input/intervals/sequence_grouping_newline/sequence_grouping_
        sed -n '/chrM:1+/,$p' {output.sg} | sed '1d'  > input/intervals/sequence_grouping_newline/sequence_grouping_25.list
        
        head -25 {output.sg_unmapped} | split --lines=1 --numeric-suffixes --additional-suffix=".list" - input/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_
        sed -n '/chrM:1+/,$p' {output.sg_unmapped} | sed '1d' > input/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_25.list

        """        
rule scatter_interval_list:
    #input: {CALLING_INTERVAL_LIST}
    input: "/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-wgs_calling_regions.hg38.interval_list"
    output: directory("output/intervals/scatter_interval_list/")
    benchmark: "benchmarks/scatter_interval_list/scatter_interval_list.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

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
intervals = sorted(glob.glob("{output}*/*.interval_list"))
for i, interval in enumerate(intervals):
  (directory, filename) = os.path.split(interval)
  newName = os.path.join(directory, str(i + 1) + filename)
  os.rename(interval, newName)
print(len(intervals))
CODE
        """


   
