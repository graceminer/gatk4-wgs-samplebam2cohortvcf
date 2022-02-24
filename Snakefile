configfile: "config.yaml"

import glob
import os

REF_FASTA = config["REF_FASTA"]

#[minerg01@li03c03 ~]$ ml bwa/0.7.8
#[minerg01@li03c03 ~]$ BWA_VERSION=$(bwa 2>&1 |     grep -e '^Version' |     sed 's/Version: //')
#[minerg01@li03c03 ~]$ echo $BWA_VERSION
#0.7.8-r455
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

#INTERVALS = glob.glob("input/sequence_grouping/sequence_grouping_*.list")
INTERVALS =  glob.glob("input/intervals/sequence_grouping_unmapped_newline/sequence_grouping.unmapped_*.list")
#known_indels_sites_vcfs
KNOWN_SITES_MILLS1000G="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_SITES_HOMOSAPIENS="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-Homo_sapiens_assembly38.known_indels.vcf.gz"



SNPEFF = config["SNPEFF"]
SNPSIFT = config["SNPSIFT"]
#GNOMAD = config["GNOMAD"]
MAGMA = config["MAGMA"]
#SNPSIFT_FILTER_IMPACT = config["SNPSIFT_FILTER_IMPACT"]
SNPSIFT_FILTER_IMPACT_TWO="(((ANN[0].IMPACT = 'HIGH') | (ANN[0].IMPACT = 'MODERATE')) & !(ANN[0].EFFECT = 'sequence_feature'))"
SNPSIFT_FILTER_BIOTYPE_PROTEIN_CODING = "(ANN[0].BIOTYPE = 'protein_coding')"
DATE=20201204
PHENO_COVARS_ALL = "output/modules/pca_final/pheno_covars_all"
GENOME = config["GENOME"]
NAME = "1096013095.final.bam"
UBAM = ["A","B","C"] #,"D", "E" ,"F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U"]
CHROMOSOMES = range(1, 23) # ["all"]
#PATH1 = os.getcwd()

#SAMPLES, = glob_wildcards("/sc/arion/projects/MMAAAS/ngs_res/gatk20200105/input/{samples}.final.bam")
#SAMPLES, = glob_wildcards("input/{samples}.final.bam")

#RGS, = glob_wildcards("output/{wildcards.samples}/{rgs}.bam")
#RGS, = glob_wildcards("output/{wildcards.samples}/{wildcards.samples}_{rgs}.bam")

#(SAMPLE,RG) = glob_wildcards("output/split_rg/{sample}/{rg}.bam")

#(SAMPLE,RG) = glob_wildcards("output/{sample}/{rg}.bam")

include: "snakefile1_revertsam"
#include: "smk/snakefilecurrent"
