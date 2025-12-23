#!/usr/bin/env python3

from expTools import execute
from common import PATH_TO_EXEC, PATH_TO_PERF_FROM_EXEC, PATH_TO_DATA_FROM_EXEC

options = {}
options["--ex, --ey, --ez"] = [100], [100], [100]
options["--timemax"] = [1.536]
options["--insitu-slice"] = [""]
options["--insitu-slice-interval"] = [2**i for i in range(3,9)]  # 8,16,32,64,128,256
options["--insitu-slice-folder"] = [PATH_TO_DATA_FROM_EXEC + "insitu_slice/"]
options["--perf-file"] = [PATH_TO_PERF_FROM_EXEC + "perf_insitu_slice.csv"]


var_env = {}

nbruns = 1

execute(f'./semproxy', var_env, options, nbruns=nbruns, verbose=True, execPath=PATH_TO_EXEC)



