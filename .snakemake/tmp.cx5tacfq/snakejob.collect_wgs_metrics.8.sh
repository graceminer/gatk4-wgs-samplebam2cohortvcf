#!/bin/sh
# properties = {"type": "single", "rule": "collect_wgs_metrics", "local": false, "input": ["output/final_bam/1102781231.final.bam"], "output": ["output/final_bam/1102781231.final_bam.wgs_metrics"], "wildcards": {"sample": "1102781231"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 8, "cluster": {"project": "acc_MMAAAS", "queue": "premium", "cores": 4, "mem": 8000, "time": "05:00"}}
 cd /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113 && \
PATH='/hpc/packages/minerva-centos7/python/3.7.3/bin':$PATH /hpc/packages/minerva-centos7/python/3.7.3/bin/python3.7 \
-m snakemake output/final_bam/1102781231.final_bam.wgs_metrics --snakefile /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113/snakefilecurrent \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113/.snakemake/tmp.cx5tacfq output/final_bam/1102781231.final.bam --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules collect_wgs_metrics --nocolor --notemp --no-hooks --nolock \
--mode 2 

