#!/usr/bin/env python3

from expTools import execute
from common import PATH_TO_EXEC, PATH_TO_PERF_FROM_EXEC, PATH_TO_DATA_FROM_EXEC

options = {}
options["--ex, --ey, --ez"] = [100], [100], [100]
options["--timemax"] = [1.536]
options["--save-snapshot"] = [""]
options["--snap-interval"] = [2**i for i in range(5,9)]  # 32,64,128,256
options["--snap-folder"] = [PATH_TO_DATA_FROM_EXEC + "snap_interval/"]
options["--perf-file"] = [PATH_TO_PERF_FROM_EXEC + "perf_snap-interval.csv"]


var_env = {}

nbruns = 1

execute(f'./semproxy', var_env, options, nbruns=nbruns, verbose=True, execPath=PATH_TO_EXEC)




