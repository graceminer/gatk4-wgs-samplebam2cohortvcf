configfile: "config/config.yaml"

import glob
import os

REF_FASTA = config["REF_FASTA"]

BWA_VERSION = "0.7.8-r455"
COMPRESSION_LEVEL=2   

# cross check fingerprints (on the aligned.duplicates_marked.sorted.bam)
LOD_THRESHOLD=-20.0
CROSSCHECK_BY="READGROUP"
HAPLOTYPE_DATABASE_FILE="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.haplotype_database.txt"

# check contamination
CONTAMINATION_SITES_UD="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.contam.UD"
CONTAMINATION_SITES_BED="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.contam.bed"
CONTAMINATION_SITES_MU="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.contam.mu"
CONTAMINATION_UNDERESTIMATION_FACTOR=0.75
CONTAMINATION_SITES_PREFIX="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.contam"
VerifyBamID = config["VerifyBamID"]

# base_recalibration_report
REF_DICT=config["REF_DICT"]
DBSNP_VCF=config["DBSNP_VCF"]

INTERVALS =  glob.glob("output/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_*.list")
#known_indels_sites_vcfs
KNOWN_SITES_MILLS1000G="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_SITES_HOMOSAPIENS="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Homo_sapiens_assembly38.known_indels.vcf.gz"




(SAMPLE,RG) = glob_wildcards("output/split_rg/{sample}/{rg}.bam")



    
localrules: main

#default combinatorial function is product, which can be replaced with other functions in this case "zip()" to create a tuple so that expand doesnt create all combos of sample and rg which ends up swapping sample dirs and rg and producing all possible cobinations as opposed to sampleA/sampleA_rg_bam

rule main3:
    input:
        expand("output/aligned_bam/{sample}/", sample=SAMPLE),
	expand("output/duplicates_marked/{sample}/{sample}.aligned.unsorted.duplicates_marked.bam", zip, sample=SAMPLE),
        expand("output/duplicates_marked/{sample}/{sample}.aligned.unsorted.duplicates_marked.duplicate_metrics",zip, sample=SAMPLE),
	expand("output/sorted_bam/{sample}/{sample}.aligned.duplicates_marked.sorted.bam", zip, sample=SAMPLE),
	expand("output/sorted_bam/{sample}/{sample}.aligned.duplicates_marked.sorted.crosscheck", zip, sample=SAMPLE), 	
	expand("output/sorted_bam/{sample}/{sample}.aligned.duplicates_marked.sorted.preBqsr.selfSM", zip, sample=SAMPLE),
        #expand("output/recalibrated_bam/{sample}/{rg}.aligned.duplicates_marked.recalibrated.bam", zip, sample=SAMPLE, rg=RG),
        #expand("output/recalibrated_bam/{sample}/{rg}.recal_data.csv", zip, sample=SAMPLE, rg=RG),
        

rule mark_duplicates:
    input: directory("output/aligned_bam/{sample}/")
    #output: bam="output/duplicates_marked/{sample}/{sample}.aligned.unsorted.duplicates_marked.bam",metrics="output/duplicates_marked/{sample}/{sample}.aligned.unsorted.duplicates_marked.duplicate_metrics"
    output: bam="output/duplicates_marked/{sample}/{sample}.aligned.unsorted.duplicates_marked.bam",metrics="output/duplicates_marked/{sample}/{sample}.aligned.unsorted.duplicates_marked.duplicate_metrics"
    shell:
        """
        module load java/1.8.0_211 python/3.7.3 picard/2.22.3

        #mkdir -p output/duplicates_marked/{wildcards.sample}/
        bamlist=$(ls {input}*.bam |  sed 's/output/INPUT=output/g')

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
    output: "output/sorted_bam/{sample}/{sample}.aligned.duplicates_marked.sorted.bam"
    shell:
        """
	module load R/3.5.3 java/1.8.0_211 python/3.7.3 picard/2.22.3

        java -Dsamjdk.compression_level={COMPRESSION_LEVEL} -Xms4000m -jar $PICARD \
          SortSam \
          INPUT={input} \
          OUTPUT={output} \
          SORT_ORDER="coordinate" \
          CREATE_INDEX=true \
          CREATE_MD5_FILE=true \
          MAX_RECORDS_IN_RAM=300000
        """
rule crosscheck_fingerprints:
    input: rules.sort_bam.output
    output: "output/sorted_bam/{sample}/{sample}.aligned.duplicates_marked.sorted.crosscheck" 
    shell:
        """
	module load R/3.5.3 java/1.8.0_211 python/3.7.3 picard/2.22.3

        java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2000m -jar $PICARD \
          CrosscheckFingerprints \
          OUTPUT={output} \
          HAPLOTYPE_MAP={HAPLOTYPE_DATABASE_FILE} \
          EXPECT_ALL_GROUPS_TO_MATCH=true \
          INPUT={input} \
          LOD_THRESHOLD={LOD_THRESHOLD} \
          CROSSCHECK_BY={CROSSCHECK_BY}
	"""


rule check_contamination:
    input: rules.sort_bam.output
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

    



                                                                                                                                                                                                                                                                                                                                         
