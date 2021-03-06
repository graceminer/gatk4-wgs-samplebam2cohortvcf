configfile: "config/config.yaml"
import glob
import os

REF_FASTA = config["REF_FASTA"]


REF_DICT=config["REF_DICT"]
DBSNP_VCF=config["DBSNP_VCF"]



CALLING_INTERVAL_LIST="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/resources-broad-hg38-v0-wgs_calling_regions.hg38.interval_list"
EVALUATION_INTERVAL_LIST="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/intervals_wgs_evaluation_regions.hg38.interval_list"
(SAMPLE,RG) = glob_wildcards("output/split_rg/{sample}/{rg}.bam")
INTERVAL_LIST = glob.glob("output/intervals/scatter_interval_list/*/*.interval_list")
INTERVALS=range(1,51)

# check unique
SAMPLE_NUM_THRESHOLD=2

# scatter interval list
SCATTER_MODE = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
REF_FASTA = config["REF_FASTA"]
UNPADDED_INTERVAL_FILE="/sc/arion/projects/MMAAAS/src/gatk-resources/hg38_20201116/hg38/v0/intervals-hg38.even.handcurated.20k.intervals"
UNBOUNDED_SCATTER_COUNT_SCALE_FACTOR=2.5

    
localrules: main
#default combinatorial function is product, which can be replaced with other functions in this case "zip()" to create a tuple so that expand doesnt create all combos of sample and rg which ends up swapping sample dirs and rg and producing all possible cobinations as opposed to sampleA/sampleA_rg_bam

rule main5:
    input:
        #expand("output/gvcf/{sample}.{intervals}.g.vcf.gz", sample=SAMPLE, intervals=INTERVALS),
        expand("output/merge_gvcf_metrics/{sample}.g.vcf.variant_calling_summary_metrics", sample=SAMPLE),
        expand("output/merge_gvcf_metrics/{sample}.validate_gvcf.metrics", sample=SAMPLE),
        expand("output/merge_gvcf_metrics/{sample}.g.vcf.variant_calling_detail_metrics", sample=SAMPLE),
	"output/merge_gvcf_metrics/check_unique.gvcf.txt",
	"output/variant_calling_scatter_interval/",


rule gatk4_haplotype_caller:
    input: bam="output/final_bam/{sample}.final.bam"
    output: "output/gvcf/{sample}.{intervals}.g.vcf.gz"
    params: interval_file="output/intervals/scatter_interval_list/*/{intervals}scattered.interval_list"        
    benchmark: "benchmarks/gvcf/{sample}.gvcf.{intervals}.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3
        #this also works and uses $interval_file as opposed to params
        #interval_file=$(ls input/scatter_interval_list/*/{wildcards.intervals}scattered.interval_list)
        set -e
        gatk --java-options "-Xms5g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
          HaplotypeCaller \
          -R {REF_FASTA} \
          -I {input.bam} \
          -L {params.interval_file} \
          -O {output} \
          -contamination 0.0 \
          -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
          -new-qual \
          -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
          -ERC GVCF 
        """

rule merge_interval_gvcfs:
    input: expand(rules.gatk4_haplotype_caller.output, sample=SAMPLE, intervals=INTERVALS) 
    output: "output/merge_gvcf/{sample}.g.vcf.gz"
    params: "output/gvcf/{sample}.*.g.vcf.gz"
    benchmark: "benchmarks/merge_gvcf/{sample}.merge_interval_gvcf.txt"
    shell:
        """
        module load R/3.5.3  gatk/4.2.0.0 java/1.8.0_211 python/3.7.3 picard/2.22.3

        input=$(ls {params} | sed 's/output/INPUT=output/g')
        java -Xms2000m -jar $PICARD \
          MergeVcfs \
          $input \
          OUTPUT={output}
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

rule check_unique:
    input: expand(rules.merge_interval_gvcfs.output, sample=SAMPLE)
    output: name_map="output/sample_name_map",stdout="output/merge_gvcf_metrics/check_unique.gvcf.txt"
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
    
rule split_interval_list:
    input: rules.check_unique.output.name_map
    output: directory("output/variant_calling_scatter_interval/")
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

#####################################################################################################                                                                                                                                                                                                                                                                                                                                         
