#!/usr/bin/env python3

from expTools import execute
from common import PATH_TO_EXEC, PATH_TO_PERF_FROM_EXEC, PATH_TO_DATA_FROM_EXEC

options = {}
options["--ex, --ey, --ez"] = [100], [100], [100]
options["--timemax"] = [1.536]
options["--insitu-histogram"] = [""]
options["--insitu-interval"] = [2**i for i in range(3,9)]  # 8,16,32,64,128,256
options["--insitu-folder"] = [PATH_TO_DATA_FROM_EXEC + "insitu_histo/"]
options["--perf-file"] = [PATH_TO_PERF_FROM_EXEC + "perf_insitu_histo.csv"]

var_env = {}

nbruns = 1

execute(f'./semproxy', var_env, options, nbruns=nbruns, verbose=True, execPath=PATH_TO_EXEC)


