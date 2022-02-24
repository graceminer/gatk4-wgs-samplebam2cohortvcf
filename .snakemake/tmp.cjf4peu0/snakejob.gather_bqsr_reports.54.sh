#!/bin/sh
# properties = {"type": "single", "rule": "gather_bqsr_reports", "local": false, "input": ["output/recalibrated_bam/1102781231.recal_data_23.csv", "output/recalibrated_bam/1102781231.recal_data_13.csv", "output/recalibrated_bam/1102781231.recal_data_04.csv", "output/recalibrated_bam/1102781231.recal_data_24.csv", "output/recalibrated_bam/1102781231.recal_data_25.csv", "output/recalibrated_bam/1102781231.recal_data_21.csv", "output/recalibrated_bam/1102781231.recal_data_11.csv", "output/recalibrated_bam/1102781231.recal_data_18.csv", "output/recalibrated_bam/1102781231.recal_data_15.csv", "output/recalibrated_bam/1102781231.recal_data_20.csv", "output/recalibrated_bam/1102781231.recal_data_06.csv", "output/recalibrated_bam/1102781231.recal_data_12.csv", "output/recalibrated_bam/1102781231.recal_data_07.csv", "output/recalibrated_bam/1102781231.recal_data_08.csv", "output/recalibrated_bam/1102781231.recal_data_02.csv", "output/recalibrated_bam/1102781231.recal_data_09.csv", "output/recalibrated_bam/1102781231.recal_data_00.csv", "output/recalibrated_bam/1102781231.recal_data_10.csv", "output/recalibrated_bam/1102781231.recal_data_03.csv", "output/recalibrated_bam/1102781231.recal_data_22.csv", "output/recalibrated_bam/1102781231.recal_data_05.csv", "output/recalibrated_bam/1102781231.recal_data_16.csv", "output/recalibrated_bam/1102781231.recal_data_14.csv", "output/recalibrated_bam/1102781231.recal_data_01.csv", "output/recalibrated_bam/1102781231.recal_data_19.csv", "output/recalibrated_bam/1102781231.recal_data_17.csv", "output/recalibrated_bam/1102781231.recal_data_23.csv", "output/recalibrated_bam/1102781231.recal_data_13.csv", "output/recalibrated_bam/1102781231.recal_data_04.csv", "output/recalibrated_bam/1102781231.recal_data_24.csv", "output/recalibrated_bam/1102781231.recal_data_25.csv", "output/recalibrated_bam/1102781231.recal_data_21.csv", "output/recalibrated_bam/1102781231.recal_data_11.csv", "output/recalibrated_bam/1102781231.recal_data_18.csv", "output/recalibrated_bam/1102781231.recal_data_15.csv", "output/recalibrated_bam/1102781231.recal_data_20.csv", "output/recalibrated_bam/1102781231.recal_data_06.csv", "output/recalibrated_bam/1102781231.recal_data_12.csv", "output/recalibrated_bam/1102781231.recal_data_07.csv", "output/recalibrated_bam/1102781231.recal_data_08.csv", "output/recalibrated_bam/1102781231.recal_data_02.csv", "output/recalibrated_bam/1102781231.recal_data_09.csv", "output/recalibrated_bam/1102781231.recal_data_00.csv", "output/recalibrated_bam/1102781231.recal_data_10.csv", "output/recalibrated_bam/1102781231.recal_data_03.csv", "output/recalibrated_bam/1102781231.recal_data_22.csv", "output/recalibrated_bam/1102781231.recal_data_05.csv", "output/recalibrated_bam/1102781231.recal_data_16.csv", "output/recalibrated_bam/1102781231.recal_data_14.csv", "output/recalibrated_bam/1102781231.recal_data_01.csv", "output/recalibrated_bam/1102781231.recal_data_19.csv", "output/recalibrated_bam/1102781231.recal_data_17.csv", "output/recalibrated_bam/1102781231.recal_data_23.csv", "output/recalibrated_bam/1102781231.recal_data_13.csv", "output/recalibrated_bam/1102781231.recal_data_04.csv", "output/recalibrated_bam/1102781231.recal_data_24.csv", "output/recalibrated_bam/1102781231.recal_data_25.csv", "output/recalibrated_bam/1102781231.recal_data_21.csv", "output/recalibrated_bam/1102781231.recal_data_11.csv", "output/recalibrated_bam/1102781231.recal_data_18.csv", "output/recalibrated_bam/1102781231.recal_data_15.csv", "output/recalibrated_bam/1102781231.recal_data_20.csv", "output/recalibrated_bam/1102781231.recal_data_06.csv", "output/recalibrated_bam/1102781231.recal_data_12.csv", "output/recalibrated_bam/1102781231.recal_data_07.csv", "output/recalibrated_bam/1102781231.recal_data_08.csv", "output/recalibrated_bam/1102781231.recal_data_02.csv", "output/recalibrated_bam/1102781231.recal_data_09.csv", "output/recalibrated_bam/1102781231.recal_data_00.csv", "output/recalibrated_bam/1102781231.recal_data_10.csv", "output/recalibrated_bam/1102781231.recal_data_03.csv", "output/recalibrated_bam/1102781231.recal_data_22.csv", "output/recalibrated_bam/1102781231.recal_data_05.csv", "output/recalibrated_bam/1102781231.recal_data_16.csv", "output/recalibrated_bam/1102781231.recal_data_14.csv", "output/recalibrated_bam/1102781231.recal_data_01.csv", "output/recalibrated_bam/1102781231.recal_data_19.csv", "output/recalibrated_bam/1102781231.recal_data_17.csv", "output/recalibrated_bam/1102781231.recal_data_23.csv", "output/recalibrated_bam/1102781231.recal_data_13.csv", "output/recalibrated_bam/1102781231.recal_data_04.csv", "output/recalibrated_bam/1102781231.recal_data_24.csv", "output/recalibrated_bam/1102781231.recal_data_25.csv", "output/recalibrated_bam/1102781231.recal_data_21.csv", "output/recalibrated_bam/1102781231.recal_data_11.csv", "output/recalibrated_bam/1102781231.recal_data_18.csv", "output/recalibrated_bam/1102781231.recal_data_15.csv", "output/recalibrated_bam/1102781231.recal_data_20.csv", "output/recalibrated_bam/1102781231.recal_data_06.csv", "output/recalibrated_bam/1102781231.recal_data_12.csv", "output/recalibrated_bam/1102781231.recal_data_07.csv", "output/recalibrated_bam/1102781231.recal_data_08.csv", "output/recalibrated_bam/1102781231.recal_data_02.csv", "output/recalibrated_bam/1102781231.recal_data_09.csv", "output/recalibrated_bam/1102781231.recal_data_00.csv", "output/recalibrated_bam/1102781231.recal_data_10.csv", "output/recalibrated_bam/1102781231.recal_data_03.csv", "output/recalibrated_bam/1102781231.recal_data_22.csv", "output/recalibrated_bam/1102781231.recal_data_05.csv", "output/recalibrated_bam/1102781231.recal_data_16.csv", "output/recalibrated_bam/1102781231.recal_data_14.csv", "output/recalibrated_bam/1102781231.recal_data_01.csv", "output/recalibrated_bam/1102781231.recal_data_19.csv", "output/recalibrated_bam/1102781231.recal_data_17.csv", "output/recalibrated_bam/1102781231.recal_data_23.csv", "output/recalibrated_bam/1102781231.recal_data_13.csv", "output/recalibrated_bam/1102781231.recal_data_04.csv", "output/recalibrated_bam/1102781231.recal_data_24.csv", "output/recalibrated_bam/1102781231.recal_data_25.csv", "output/recalibrated_bam/1102781231.recal_data_21.csv", "output/recalibrated_bam/1102781231.recal_data_11.csv", "output/recalibrated_bam/1102781231.recal_data_18.csv", "output/recalibrated_bam/1102781231.recal_data_15.csv", "output/recalibrated_bam/1102781231.recal_data_20.csv", "output/recalibrated_bam/1102781231.recal_data_06.csv", "output/recalibrated_bam/1102781231.recal_data_12.csv", "output/recalibrated_bam/1102781231.recal_data_07.csv", "output/recalibrated_bam/1102781231.recal_data_08.csv", "output/recalibrated_bam/1102781231.recal_data_02.csv", "output/recalibrated_bam/1102781231.recal_data_09.csv", "output/recalibrated_bam/1102781231.recal_data_00.csv", "output/recalibrated_bam/1102781231.recal_data_10.csv", "output/recalibrated_bam/1102781231.recal_data_03.csv", "output/recalibrated_bam/1102781231.recal_data_22.csv", "output/recalibrated_bam/1102781231.recal_data_05.csv", "output/recalibrated_bam/1102781231.recal_data_16.csv", "output/recalibrated_bam/1102781231.recal_data_14.csv", "output/recalibrated_bam/1102781231.recal_data_01.csv", "output/recalibrated_bam/1102781231.recal_data_19.csv", "output/recalibrated_bam/1102781231.recal_data_17.csv", "output/recalibrated_bam/1102777134.recal_data_23.csv", "output/recalibrated_bam/1102777134.recal_data_13.csv", "output/recalibrated_bam/1102777134.recal_data_04.csv", "output/recalibrated_bam/1102777134.recal_data_24.csv", "output/recalibrated_bam/1102777134.recal_data_25.csv", "output/recalibrated_bam/1102777134.recal_data_21.csv", "output/recalibrated_bam/1102777134.recal_data_11.csv", "output/recalibrated_bam/1102777134.recal_data_18.csv", "output/recalibrated_bam/1102777134.recal_data_15.csv", "output/recalibrated_bam/1102777134.recal_data_20.csv", "output/recalibrated_bam/1102777134.recal_data_06.csv", "output/recalibrated_bam/1102777134.recal_data_12.csv", "output/recalibrated_bam/1102777134.recal_data_07.csv", "output/recalibrated_bam/1102777134.recal_data_08.csv", "output/recalibrated_bam/1102777134.recal_data_02.csv", "output/recalibrated_bam/1102777134.recal_data_09.csv", "output/recalibrated_bam/1102777134.recal_data_00.csv", "output/recalibrated_bam/1102777134.recal_data_10.csv", "output/recalibrated_bam/1102777134.recal_data_03.csv", "output/recalibrated_bam/1102777134.recal_data_22.csv", "output/recalibrated_bam/1102777134.recal_data_05.csv", "output/recalibrated_bam/1102777134.recal_data_16.csv", "output/recalibrated_bam/1102777134.recal_data_14.csv", "output/recalibrated_bam/1102777134.recal_data_01.csv", "output/recalibrated_bam/1102777134.recal_data_19.csv", "output/recalibrated_bam/1102777134.recal_data_17.csv", "output/recalibrated_bam/1102777134.recal_data_23.csv", "output/recalibrated_bam/1102777134.recal_data_13.csv", "output/recalibrated_bam/1102777134.recal_data_04.csv", "output/recalibrated_bam/1102777134.recal_data_24.csv", "output/recalibrated_bam/1102777134.recal_data_25.csv", "output/recalibrated_bam/1102777134.recal_data_21.csv", "output/recalibrated_bam/1102777134.recal_data_11.csv", "output/recalibrated_bam/1102777134.recal_data_18.csv", "output/recalibrated_bam/1102777134.recal_data_15.csv", "output/recalibrated_bam/1102777134.recal_data_20.csv", "output/recalibrated_bam/1102777134.recal_data_06.csv", "output/recalibrated_bam/1102777134.recal_data_12.csv", "output/recalibrated_bam/1102777134.recal_data_07.csv", "output/recalibrated_bam/1102777134.recal_data_08.csv", "output/recalibrated_bam/1102777134.recal_data_02.csv", "output/recalibrated_bam/1102777134.recal_data_09.csv", "output/recalibrated_bam/1102777134.recal_data_00.csv", "output/recalibrated_bam/1102777134.recal_data_10.csv", "output/recalibrated_bam/1102777134.recal_data_03.csv", "output/recalibrated_bam/1102777134.recal_data_22.csv", "output/recalibrated_bam/1102777134.recal_data_05.csv", "output/recalibrated_bam/1102777134.recal_data_16.csv", "output/recalibrated_bam/1102777134.recal_data_14.csv", "output/recalibrated_bam/1102777134.recal_data_01.csv", "output/recalibrated_bam/1102777134.recal_data_19.csv", "output/recalibrated_bam/1102777134.recal_data_17.csv", "output/recalibrated_bam/1102777134.recal_data_23.csv", "output/recalibrated_bam/1102777134.recal_data_13.csv", "output/recalibrated_bam/1102777134.recal_data_04.csv", "output/recalibrated_bam/1102777134.recal_data_24.csv", "output/recalibrated_bam/1102777134.recal_data_25.csv", "output/recalibrated_bam/1102777134.recal_data_21.csv", "output/recalibrated_bam/1102777134.recal_data_11.csv", "output/recalibrated_bam/1102777134.recal_data_18.csv", "output/recalibrated_bam/1102777134.recal_data_15.csv", "output/recalibrated_bam/1102777134.recal_data_20.csv", "output/recalibrated_bam/1102777134.recal_data_06.csv", "output/recalibrated_bam/1102777134.recal_data_12.csv", "output/recalibrated_bam/1102777134.recal_data_07.csv", "output/recalibrated_bam/1102777134.recal_data_08.csv", "output/recalibrated_bam/1102777134.recal_data_02.csv", "output/recalibrated_bam/1102777134.recal_data_09.csv", "output/recalibrated_bam/1102777134.recal_data_00.csv", "output/recalibrated_bam/1102777134.recal_data_10.csv", "output/recalibrated_bam/1102777134.recal_data_03.csv", "output/recalibrated_bam/1102777134.recal_data_22.csv", "output/recalibrated_bam/1102777134.recal_data_05.csv", "output/recalibrated_bam/1102777134.recal_data_16.csv", "output/recalibrated_bam/1102777134.recal_data_14.csv", "output/recalibrated_bam/1102777134.recal_data_01.csv", "output/recalibrated_bam/1102777134.recal_data_19.csv", "output/recalibrated_bam/1102777134.recal_data_17.csv", "output/recalibrated_bam/1102777134.recal_data_23.csv", "output/recalibrated_bam/1102777134.recal_data_13.csv", "output/recalibrated_bam/1102777134.recal_data_04.csv", "output/recalibrated_bam/1102777134.recal_data_24.csv", "output/recalibrated_bam/1102777134.recal_data_25.csv", "output/recalibrated_bam/1102777134.recal_data_21.csv", "output/recalibrated_bam/1102777134.recal_data_11.csv", "output/recalibrated_bam/1102777134.recal_data_18.csv", "output/recalibrated_bam/1102777134.recal_data_15.csv", "output/recalibrated_bam/1102777134.recal_data_20.csv", "output/recalibrated_bam/1102777134.recal_data_06.csv", "output/recalibrated_bam/1102777134.recal_data_12.csv", "output/recalibrated_bam/1102777134.recal_data_07.csv", "output/recalibrated_bam/1102777134.recal_data_08.csv", "output/recalibrated_bam/1102777134.recal_data_02.csv", "output/recalibrated_bam/1102777134.recal_data_09.csv", "output/recalibrated_bam/1102777134.recal_data_00.csv", "output/recalibrated_bam/1102777134.recal_data_10.csv", "output/recalibrated_bam/1102777134.recal_data_03.csv", "output/recalibrated_bam/1102777134.recal_data_22.csv", "output/recalibrated_bam/1102777134.recal_data_05.csv", "output/recalibrated_bam/1102777134.recal_data_16.csv", "output/recalibrated_bam/1102777134.recal_data_14.csv", "output/recalibrated_bam/1102777134.recal_data_01.csv", "output/recalibrated_bam/1102777134.recal_data_19.csv", "output/recalibrated_bam/1102777134.recal_data_17.csv", "output/recalibrated_bam/1102777134.recal_data_23.csv", "output/recalibrated_bam/1102777134.recal_data_13.csv", "output/recalibrated_bam/1102777134.recal_data_04.csv", "output/recalibrated_bam/1102777134.recal_data_24.csv", "output/recalibrated_bam/1102777134.recal_data_25.csv", "output/recalibrated_bam/1102777134.recal_data_21.csv", "output/recalibrated_bam/1102777134.recal_data_11.csv", "output/recalibrated_bam/1102777134.recal_data_18.csv", "output/recalibrated_bam/1102777134.recal_data_15.csv", "output/recalibrated_bam/1102777134.recal_data_20.csv", "output/recalibrated_bam/1102777134.recal_data_06.csv", "output/recalibrated_bam/1102777134.recal_data_12.csv", "output/recalibrated_bam/1102777134.recal_data_07.csv", "output/recalibrated_bam/1102777134.recal_data_08.csv", "output/recalibrated_bam/1102777134.recal_data_02.csv", "output/recalibrated_bam/1102777134.recal_data_09.csv", "output/recalibrated_bam/1102777134.recal_data_00.csv", "output/recalibrated_bam/1102777134.recal_data_10.csv", "output/recalibrated_bam/1102777134.recal_data_03.csv", "output/recalibrated_bam/1102777134.recal_data_22.csv", "output/recalibrated_bam/1102777134.recal_data_05.csv", "output/recalibrated_bam/1102777134.recal_data_16.csv", "output/recalibrated_bam/1102777134.recal_data_14.csv", "output/recalibrated_bam/1102777134.recal_data_01.csv", "output/recalibrated_bam/1102777134.recal_data_19.csv", "output/recalibrated_bam/1102777134.recal_data_17.csv"], "output": ["output/recalibrated_bam/1102777134.recal_data.csv"], "wildcards": {"sample": "1102777134"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 54, "cluster": {"project": "acc_MMAAAS", "queue": "premium", "cores": 4, "mem": 8000, "time": "04:00"}}
 cd /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113 && \
PATH='/hpc/packages/minerva-centos7/python/3.7.3/bin':$PATH /hpc/packages/minerva-centos7/python/3.7.3/bin/python3.7 \
-m snakemake output/recalibrated_bam/1102777134.recal_data.csv --snakefile /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113/snakefilecurrent \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /sc/arion/projects/MMAAAS/ngs_res/aaa_pilot_20201113/gatk20220113/.snakemake/tmp.cjf4peu0 output/recalibrated_bam/1102781231.recal_data_23.csv output/recalibrated_bam/1102781231.recal_data_13.csv output/recalibrated_bam/1102781231.recal_data_04.csv output/recalibrated_bam/1102781231.recal_data_24.csv output/recalibrated_bam/1102781231.recal_data_25.csv output/recalibrated_bam/1102781231.recal_data_21.csv output/recalibrated_bam/1102781231.recal_data_11.csv output/recalibrated_bam/1102781231.recal_data_18.csv output/recalibrated_bam/1102781231.recal_data_15.csv output/recalibrated_bam/1102781231.recal_data_20.csv output/recalibrated_bam/1102781231.recal_data_06.csv output/recalibrated_bam/1102781231.recal_data_12.csv output/recalibrated_bam/1102781231.recal_data_07.csv output/recalibrated_bam/1102781231.recal_data_08.csv output/recalibrated_bam/1102781231.recal_data_02.csv output/recalibrated_bam/1102781231.recal_data_09.csv output/recalibrated_bam/1102781231.recal_data_00.csv output/recalibrated_bam/1102781231.recal_data_10.csv output/recalibrated_bam/1102781231.recal_data_03.csv output/recalibrated_bam/1102781231.recal_data_22.csv output/recalibrated_bam/1102781231.recal_data_05.csv output/recalibrated_bam/1102781231.recal_data_16.csv output/recalibrated_bam/1102781231.recal_data_14.csv output/recalibrated_bam/1102781231.recal_data_01.csv output/recalibrated_bam/1102781231.recal_data_19.csv output/recalibrated_bam/1102781231.recal_data_17.csv output/recalibrated_bam/1102781231.recal_data_23.csv output/recalibrated_bam/1102781231.recal_data_13.csv output/recalibrated_bam/1102781231.recal_data_04.csv output/recalibrated_bam/1102781231.recal_data_24.csv output/recalibrated_bam/1102781231.recal_data_25.csv output/recalibrated_bam/1102781231.recal_data_21.csv output/recalibrated_bam/1102781231.recal_data_11.csv output/recalibrated_bam/1102781231.recal_data_18.csv output/recalibrated_bam/1102781231.recal_data_15.csv output/recalibrated_bam/1102781231.recal_data_20.csv output/recalibrated_bam/1102781231.recal_data_06.csv output/recalibrated_bam/1102781231.recal_data_12.csv output/recalibrated_bam/1102781231.recal_data_07.csv output/recalibrated_bam/1102781231.recal_data_08.csv output/recalibrated_bam/1102781231.recal_data_02.csv output/recalibrated_bam/1102781231.recal_data_09.csv output/recalibrated_bam/1102781231.recal_data_00.csv output/recalibrated_bam/1102781231.recal_data_10.csv output/recalibrated_bam/1102781231.recal_data_03.csv output/recalibrated_bam/1102781231.recal_data_22.csv output/recalibrated_bam/1102781231.recal_data_05.csv output/recalibrated_bam/1102781231.recal_data_16.csv output/recalibrated_bam/1102781231.recal_data_14.csv output/recalibrated_bam/1102781231.recal_data_01.csv output/recalibrated_bam/1102781231.recal_data_19.csv output/recalibrated_bam/1102781231.recal_data_17.csv output/recalibrated_bam/1102781231.recal_data_23.csv output/recalibrated_bam/1102781231.recal_data_13.csv output/recalibrated_bam/1102781231.recal_data_04.csv output/recalibrated_bam/1102781231.recal_data_24.csv output/recalibrated_bam/1102781231.recal_data_25.csv output/recalibrated_bam/1102781231.recal_data_21.csv output/recalibrated_bam/1102781231.recal_data_11.csv output/recalibrated_bam/1102781231.recal_data_18.csv output/recalibrated_bam/1102781231.recal_data_15.csv output/recalibrated_bam/1102781231.recal_data_20.csv output/recalibrated_bam/1102781231.recal_data_06.csv output/recalibrated_bam/1102781231.recal_data_12.csv output/recalibrated_bam/1102781231.recal_data_07.csv output/recalibrated_bam/1102781231.recal_data_08.csv output/recalibrated_bam/1102781231.recal_data_02.csv output/recalibrated_bam/1102781231.recal_data_09.csv output/recalibrated_bam/1102781231.recal_data_00.csv output/recalibrated_bam/1102781231.recal_data_10.csv output/recalibrated_bam/1102781231.recal_data_03.csv output/recalibrated_bam/1102781231.recal_data_22.csv output/recalibrated_bam/1102781231.recal_data_05.csv output/recalibrated_bam/1102781231.recal_data_16.csv output/recalibrated_bam/1102781231.recal_data_14.csv output/recalibrated_bam/1102781231.recal_data_01.csv output/recalibrated_bam/1102781231.recal_data_19.csv output/recalibrated_bam/1102781231.recal_data_17.csv output/recalibrated_bam/1102781231.recal_data_23.csv output/recalibrated_bam/1102781231.recal_data_13.csv output/recalibrated_bam/1102781231.recal_data_04.csv output/recalibrated_bam/1102781231.recal_data_24.csv output/recalibrated_bam/1102781231.recal_data_25.csv output/recalibrated_bam/1102781231.recal_data_21.csv output/recalibrated_bam/1102781231.recal_data_11.csv output/recalibrated_bam/1102781231.recal_data_18.csv output/recalibrated_bam/1102781231.recal_data_15.csv output/recalibrated_bam/1102781231.recal_data_20.csv output/recalibrated_bam/1102781231.recal_data_06.csv output/recalibrated_bam/1102781231.recal_data_12.csv output/recalibrated_bam/1102781231.recal_data_07.csv output/recalibrated_bam/1102781231.recal_data_08.csv output/recalibrated_bam/1102781231.recal_data_02.csv output/recalibrated_bam/1102781231.recal_data_09.csv output/recalibrated_bam/1102781231.recal_data_00.csv output/recalibrated_bam/1102781231.recal_data_10.csv output/recalibrated_bam/1102781231.recal_data_03.csv output/recalibrated_bam/1102781231.recal_data_22.csv output/recalibrated_bam/1102781231.recal_data_05.csv output/recalibrated_bam/1102781231.recal_data_16.csv output/recalibrated_bam/1102781231.recal_data_14.csv output/recalibrated_bam/1102781231.recal_data_01.csv output/recalibrated_bam/1102781231.recal_data_19.csv output/recalibrated_bam/1102781231.recal_data_17.csv output/recalibrated_bam/1102781231.recal_data_23.csv output/recalibrated_bam/1102781231.recal_data_13.csv output/recalibrated_bam/1102781231.recal_data_04.csv output/recalibrated_bam/1102781231.recal_data_24.csv output/recalibrated_bam/1102781231.recal_data_25.csv output/recalibrated_bam/1102781231.recal_data_21.csv output/recalibrated_bam/1102781231.recal_data_11.csv output/recalibrated_bam/1102781231.recal_data_18.csv output/recalibrated_bam/1102781231.recal_data_15.csv output/recalibrated_bam/1102781231.recal_data_20.csv output/recalibrated_bam/1102781231.recal_data_06.csv output/recalibrated_bam/1102781231.recal_data_12.csv output/recalibrated_bam/1102781231.recal_data_07.csv output/recalibrated_bam/1102781231.recal_data_08.csv output/recalibrated_bam/1102781231.recal_data_02.csv output/recalibrated_bam/1102781231.recal_data_09.csv output/recalibrated_bam/1102781231.recal_data_00.csv output/recalibrated_bam/1102781231.recal_data_10.csv output/recalibrated_bam/1102781231.recal_data_03.csv output/recalibrated_bam/1102781231.recal_data_22.csv output/recalibrated_bam/1102781231.recal_data_05.csv output/recalibrated_bam/1102781231.recal_data_16.csv output/recalibrated_bam/1102781231.recal_data_14.csv output/recalibrated_bam/1102781231.recal_data_01.csv output/recalibrated_bam/1102781231.recal_data_19.csv output/recalibrated_bam/1102781231.recal_data_17.csv output/recalibrated_bam/1102777134.recal_data_23.csv output/recalibrated_bam/1102777134.recal_data_13.csv output/recalibrated_bam/1102777134.recal_data_04.csv output/recalibrated_bam/1102777134.recal_data_24.csv output/recalibrated_bam/1102777134.recal_data_25.csv output/recalibrated_bam/1102777134.recal_data_21.csv output/recalibrated_bam/1102777134.recal_data_11.csv output/recalibrated_bam/1102777134.recal_data_18.csv output/recalibrated_bam/1102777134.recal_data_15.csv output/recalibrated_bam/1102777134.recal_data_20.csv output/recalibrated_bam/1102777134.recal_data_06.csv output/recalibrated_bam/1102777134.recal_data_12.csv output/recalibrated_bam/1102777134.recal_data_07.csv output/recalibrated_bam/1102777134.recal_data_08.csv output/recalibrated_bam/1102777134.recal_data_02.csv output/recalibrated_bam/1102777134.recal_data_09.csv output/recalibrated_bam/1102777134.recal_data_00.csv output/recalibrated_bam/1102777134.recal_data_10.csv output/recalibrated_bam/1102777134.recal_data_03.csv output/recalibrated_bam/1102777134.recal_data_22.csv output/recalibrated_bam/1102777134.recal_data_05.csv output/recalibrated_bam/1102777134.recal_data_16.csv output/recalibrated_bam/1102777134.recal_data_14.csv output/recalibrated_bam/1102777134.recal_data_01.csv output/recalibrated_bam/1102777134.recal_data_19.csv output/recalibrated_bam/1102777134.recal_data_17.csv output/recalibrated_bam/1102777134.recal_data_23.csv output/recalibrated_bam/1102777134.recal_data_13.csv output/recalibrated_bam/1102777134.recal_data_04.csv output/recalibrated_bam/1102777134.recal_data_24.csv output/recalibrated_bam/1102777134.recal_data_25.csv output/recalibrated_bam/1102777134.recal_data_21.csv output/recalibrated_bam/1102777134.recal_data_11.csv output/recalibrated_bam/1102777134.recal_data_18.csv output/recalibrated_bam/1102777134.recal_data_15.csv output/recalibrated_bam/1102777134.recal_data_20.csv output/recalibrated_bam/1102777134.recal_data_06.csv output/recalibrated_bam/1102777134.recal_data_12.csv output/recalibrated_bam/1102777134.recal_data_07.csv output/recalibrated_bam/1102777134.recal_data_08.csv output/recalibrated_bam/1102777134.recal_data_02.csv output/recalibrated_bam/1102777134.recal_data_09.csv output/recalibrated_bam/1102777134.recal_data_00.csv output/recalibrated_bam/1102777134.recal_data_10.csv output/recalibrated_bam/1102777134.recal_data_03.csv output/recalibrated_bam/1102777134.recal_data_22.csv output/recalibrated_bam/1102777134.recal_data_05.csv output/recalibrated_bam/1102777134.recal_data_16.csv output/recalibrated_bam/1102777134.recal_data_14.csv output/recalibrated_bam/1102777134.recal_data_01.csv output/recalibrated_bam/1102777134.recal_data_19.csv output/recalibrated_bam/1102777134.recal_data_17.csv output/recalibrated_bam/1102777134.recal_data_23.csv output/recalibrated_bam/1102777134.recal_data_13.csv output/recalibrated_bam/1102777134.recal_data_04.csv output/recalibrated_bam/1102777134.recal_data_24.csv output/recalibrated_bam/1102777134.recal_data_25.csv output/recalibrated_bam/1102777134.recal_data_21.csv output/recalibrated_bam/1102777134.recal_data_11.csv output/recalibrated_bam/1102777134.recal_data_18.csv output/recalibrated_bam/1102777134.recal_data_15.csv output/recalibrated_bam/1102777134.recal_data_20.csv output/recalibrated_bam/1102777134.recal_data_06.csv output/recalibrated_bam/1102777134.recal_data_12.csv output/recalibrated_bam/1102777134.recal_data_07.csv output/recalibrated_bam/1102777134.recal_data_08.csv output/recalibrated_bam/1102777134.recal_data_02.csv output/recalibrated_bam/1102777134.recal_data_09.csv output/recalibrated_bam/1102777134.recal_data_00.csv output/recalibrated_bam/1102777134.recal_data_10.csv output/recalibrated_bam/1102777134.recal_data_03.csv output/recalibrated_bam/1102777134.recal_data_22.csv output/recalibrated_bam/1102777134.recal_data_05.csv output/recalibrated_bam/1102777134.recal_data_16.csv output/recalibrated_bam/1102777134.recal_data_14.csv output/recalibrated_bam/1102777134.recal_data_01.csv output/recalibrated_bam/1102777134.recal_data_19.csv output/recalibrated_bam/1102777134.recal_data_17.csv output/recalibrated_bam/1102777134.recal_data_23.csv output/recalibrated_bam/1102777134.recal_data_13.csv output/recalibrated_bam/1102777134.recal_data_04.csv output/recalibrated_bam/1102777134.recal_data_24.csv output/recalibrated_bam/1102777134.recal_data_25.csv output/recalibrated_bam/1102777134.recal_data_21.csv output/recalibrated_bam/1102777134.recal_data_11.csv output/recalibrated_bam/1102777134.recal_data_18.csv output/recalibrated_bam/1102777134.recal_data_15.csv output/recalibrated_bam/1102777134.recal_data_20.csv output/recalibrated_bam/1102777134.recal_data_06.csv output/recalibrated_bam/1102777134.recal_data_12.csv output/recalibrated_bam/1102777134.recal_data_07.csv output/recalibrated_bam/1102777134.recal_data_08.csv output/recalibrated_bam/1102777134.recal_data_02.csv output/recalibrated_bam/1102777134.recal_data_09.csv output/recalibrated_bam/1102777134.recal_data_00.csv output/recalibrated_bam/1102777134.recal_data_10.csv output/recalibrated_bam/1102777134.recal_data_03.csv output/recalibrated_bam/1102777134.recal_data_22.csv output/recalibrated_bam/1102777134.recal_data_05.csv output/recalibrated_bam/1102777134.recal_data_16.csv output/recalibrated_bam/1102777134.recal_data_14.csv output/recalibrated_bam/1102777134.recal_data_01.csv output/recalibrated_bam/1102777134.recal_data_19.csv output/recalibrated_bam/1102777134.recal_data_17.csv output/recalibrated_bam/1102777134.recal_data_23.csv output/recalibrated_bam/1102777134.recal_data_13.csv output/recalibrated_bam/1102777134.recal_data_04.csv output/recalibrated_bam/1102777134.recal_data_24.csv output/recalibrated_bam/1102777134.recal_data_25.csv output/recalibrated_bam/1102777134.recal_data_21.csv output/recalibrated_bam/1102777134.recal_data_11.csv output/recalibrated_bam/1102777134.recal_data_18.csv output/recalibrated_bam/1102777134.recal_data_15.csv output/recalibrated_bam/1102777134.recal_data_20.csv output/recalibrated_bam/1102777134.recal_data_06.csv output/recalibrated_bam/1102777134.recal_data_12.csv output/recalibrated_bam/1102777134.recal_data_07.csv output/recalibrated_bam/1102777134.recal_data_08.csv output/recalibrated_bam/1102777134.recal_data_02.csv output/recalibrated_bam/1102777134.recal_data_09.csv output/recalibrated_bam/1102777134.recal_data_00.csv output/recalibrated_bam/1102777134.recal_data_10.csv output/recalibrated_bam/1102777134.recal_data_03.csv output/recalibrated_bam/1102777134.recal_data_22.csv output/recalibrated_bam/1102777134.recal_data_05.csv output/recalibrated_bam/1102777134.recal_data_16.csv output/recalibrated_bam/1102777134.recal_data_14.csv output/recalibrated_bam/1102777134.recal_data_01.csv output/recalibrated_bam/1102777134.recal_data_19.csv output/recalibrated_bam/1102777134.recal_data_17.csv --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules gather_bqsr_reports --nocolor --notemp --no-hooks --nolock \
--mode 2 

