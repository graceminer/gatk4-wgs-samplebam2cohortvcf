#!/bin/sh
# properties = {"type": "single", "rule": "validate_cram", "local": false, "input": ["output/final_cram/1102777134.final.cram"], "output": ["output/final_cram/1102777134.final_cram.validation_report"], "wildcards": {"sample": "1102777134"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 16, "cluster": {"project": "acc_MMAAAS", "queue": "premium", "cores": 4, "mem": 8000, "time": "04:00"}}
 cd /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113 && \
PATH='/hpc/packages/minerva-centos7/python/3.7.3/bin':$PATH /hpc/packages/minerva-centos7/python/3.7.3/bin/python3.7 \
-m snakemake output/final_cram/1102777134.final_cram.validation_report --snakefile /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113/snakefile4_validate_cram \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113/.snakemake/tmp.tnnopj8s output/final_cram/1102777134.final.cram --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules validate_cram --nocolor --notemp --no-hooks --nolock \
--mode 2 

