#!/usr/bin/env python3

from expTools import execute
from common import PATH_TO_EXEC, PATH_TO_PERF_FROM_EXEC, PATH_TO_DATA_FROM_EXEC

options = {}
options["--ex, --ey, --ez"] = [10], [10], [10]
options["--timemax"] = [1.536]
options["--perf-file"] = [PATH_TO_PERF_FROM_EXEC + "perf_insitu_slice.csv"]


var_env = {}

nbruns = 1

execute(f'./semproxy', var_env, options, nbruns=nbruns, verbose=True, execPath=PATH_TO_EXEC)

###### with in-situ slice ######

options["--insitu-slice"] = [""]
options["--insitu-slice-interval"] = [32]
options["--insitu-slice-folder"] = [PATH_TO_DATA_FROM_EXEC + "insitu_slice/"]

execute(f'./semproxy', var_env, options, nbruns=nbruns, verbose=True, execPath=PATH_TO_EXEC)


