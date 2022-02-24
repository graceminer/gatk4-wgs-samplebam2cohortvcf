#!/bin/sh
# properties = {"type": "single", "rule": "gatk4_haplotype_caller", "local": false, "input": ["output/final_bam/1102777134.final.bam"], "output": ["output/gvcf/1102777134.17.g.vcf.gz"], "wildcards": {"sample": "1102777134", "intervals": "17"}, "params": {"interval_file": "input/scatter_interval_list/*/17scattered.interval_list"}, "log": [], "threads": 1, "resources": {}, "jobid": 67, "cluster": {"project": "acc_MMAAAS", "queue": "premium", "cores": 4, "mem": 8000, "time": "04:00"}}
 cd /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113 && \
PATH='/hpc/packages/minerva-centos7/python/3.7.3/bin':$PATH /hpc/packages/minerva-centos7/python/3.7.3/bin/python3.7 \
-m snakemake output/gvcf/1102777134.17.g.vcf.gz --snakefile /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113/snakefilecurrent \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113/.snakemake/tmp.l2jcm2hx output/final_bam/1102777134.final.bam --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules gatk4_haplotype_caller --nocolor --notemp --no-hooks --nolock \
--mode 2 

