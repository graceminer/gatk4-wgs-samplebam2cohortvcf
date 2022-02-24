#!/bin/sh
# properties = {"type": "single", "rule": "convert_to_cram", "local": false, "input": ["output/final_bam/1102777134.final.bam"], "output": ["output/final_cram/1102777134.final.cram", "output/final_cram/1102777134.final.cram.md5"], "wildcards": {"sample": "1102777134"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 20, "cluster": {"project": "acc_MMAAAS", "queue": "premium", "cores": 4, "mem": 8000, "time": "08:00"}}
 cd /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113 && \
PATH='/hpc/packages/minerva-centos7/python/3.7.3/bin':$PATH /hpc/packages/minerva-centos7/python/3.7.3/bin/python3.7 \
-m snakemake output/final_cram/1102777134.final.cram --snakefile /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113/snakefilecurrent \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113/.snakemake/tmp.cx5tacfq output/final_bam/1102777134.final.bam --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules convert_to_cram --nocolor --notemp --no-hooks --nolock \
--mode 2 

