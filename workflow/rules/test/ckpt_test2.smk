configfile: "config/config.yaml"

import glob
import os

# REF_FASTA = config["REF_FASTA"]

# #[minerg01@li03c03 ~]$ ml bwa/0.7.8
# #[minerg01@li03c03 ~]$ BWA_VERSION=$(bwa 2>&1 |     grep -e '^Version' |     sed 's/Version: //')
# #[minerg01@li03c03 ~]$ echo $BWA_VERSION
# #0.7.8-r455
# # BWA_VERSION = "0.7.8-r455"
# # COMPRESSION_LEVEL=2   

# # # cross check fingerprints (on the aligned.duplicates_marked.sorted.bam)
# # LOD_THRESHOLD=-20.0
# # CROSSCHECK_BY="READGROUP"
# # HAPLOTYPE_DATABASE_FILE="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.haplotype_database.txt"

# # # check contamination
# # CONTAMINATION_SITES_UD="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.contam.UD"
# # CONTAMINATION_SITES_BED="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.contam.bed"
# # CONTAMINATION_SITES_MU="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.contam.mu"
# # CONTAMINATION_UNDERESTIMATION_FACTOR=0.75
# # CONTAMINATION_SITES_PREFIX="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/references-hg38-v0-Homo_sapiens_assembly38.contam"
# # VerifyBamID = config["VerifyBamID"]

# # # base_recalibration_report
# # REF_DICT=config["REF_DICT"]
# # DBSNP_VCF=config["DBSNP_VCF"]

# # #INTERVALS = glob.glob("input/sequence_grouping/sequence_grouping_*.list")
# # INTERVALS =  glob.glob("input/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_*.list")
# # #known_indels_sites_vcfs
# # KNOWN_SITES_MILLS1000G="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
# # KNOWN_SITES_HOMOSAPIENS="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Homo_sapiens_assembly38.known_indels.vcf.gz"



##SNPEFF = config["SNPEFF"]
#SNPSIFT = config["SNPSIFT"]
##GNOMAD = config["GNOMAD"]
#MAGMA = config["MAGMA"]
##SNPSIFT_FILTER_IMPACT = config["SNPSIFT_FILTER_IMPACT"]
#SNPSIFT_FILTER_IMPACT_TWO="(((ANN[0].IMPACT = 'HIGH') | (ANN[0].IMPACT = 'MODERATE')) & !(ANN[0].EFFECT = 'sequence_feature'))"
#SNPSIFT_FILTER_BIOTYPE_PROTEIN_CODING = "(ANN[0].BIOTYPE = 'protein_coding')"
#DATE=20201204
#PHENO_COVARS_ALL = "output/modules/pca_final/pheno_covars_all"
#GENOME = config["GENOME"]
#NAME = "1096013095.final.bam"
#UBAM = ["A","B","C"] #,"D", "E" ,"F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U"]
#CHROMOSOMES = range(1, 23) # ["all"]
#PATH1 = os.getcwd()



#(SAMPLE,RG) = glob_wildcards("output/split_rg/{sample}/{rg}.bam")


# scattergather:
#     split=25
    
wildcard_constraints:
    chr="\d+"
    
#localrules: main
#ruleorder: revertsam > sortsam    

#rule main2:
#    input:
#        expand("output/unmapped_bam/{sample}/{rg}.unmapped.bam",zip, sample=SAMPLE,rg=RG), #default combinatorial function is product, which can be replaced with other functions in this case "zip()" to create a tuple so that expand doesnt create all combos of sample and rg which ends up swapping sample dirs and rg and producing all possible cobinations as opposed to sampleA/sampleA_rg_bam
#	expand("output/unmapped_bam/{sample}/{rg}.quality_yield_metrics",zip, sample=SAMPLE,rg=RG),
#	expand("output/aligned_bam/{sample}/{rg}.aligned.unsorted.bam",zip, sample=SAMPLE,rg=RG),
#	expand("output/aligned_bam/{sample}/{rg}.aligned.unsorted.bwa.stderr.log",zip, sample=SAMPLE,rg=RG),


        
# rule sortsam2:
#     input: lambda wildcards: [os.path.join('output/','split_rg/',SAMPLE[i], x + '.bam') for i,x in enumerate(RG) if x == wildcards.rg]
#     output: "output/unmapped_bam/{sample}/{rg}.unmapped.bam"
#     shell:
#         """
#         module load java/1.8.0_211  python/3.7.3  gatk/4.2.0.0
        
# 	gatk --java-options "-Xmx3g" \
#         SortSam \
#             --INPUT {input} \
#             --OUTPUT {output} \
#             --SORT_ORDER queryname \
#             --MAX_RECORDS_IN_RAM 1000000

#         module unload java/1.8.0_211  python/3.7.3 gatk/4.2.0.0
#         """

# rule collect_quality_yield_metrics:
#     input: rules.sortsam.output
#     output: "output/unmapped_bam/{sample}/{rg}.quality_yield_metrics"
#     #output: "output/{sample}/{rg}.unmapped.bam.quality_yield_metrics"
#     shell:
#         """
#         module load java/1.8.0_211  python/3.7.3 picard/2.22.3
		
#      	java -Xms2000m -jar $PICARD \
#              CollectQualityYieldMetrics \
#              INPUT={input} \
#              OQ=true  \
#              OUTPUT={output}
# 	 module unload java/1.8.0_211  python/3.7.3 picard/2.22.3
#         """
#REDO THE align_bam_bwa rule so that the log files are part of the output (they arent erased on job fail)
#separatte the lines of the sam2fq and merge alignment steps to clean up
rule align_bam_bwa:
    input: "output/unmapped_bam/{sample}/{rg}.unmapped.bam"
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

def aggregate_input(wildcards):
   checkpoint_output = checkpoints.revertsam.get(**wildcards).output[0]
   x=expand("output/aligned_bam/{sample}/{rg}.aligned.unsorted.bam",
       sample=wildcards.sample,
       rg=glob_wildcards(os.path.join(checkpoint_output,"{rg}.bam")).rg)
   print(x)
   return x

rule readgroup_aggregate:
    input: aggregate_input
    output: "output/aligned_bam/{sample}.aligned.unsorted.bam.readgroup_list.txt"
    shell: "ls {input} > {output}"    



                                                                                                                                                                                                                                                                                                                                         
