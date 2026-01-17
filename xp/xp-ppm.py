#!/usr/bin/env python3

from expTools import execute
from common import PATH_TO_EXEC, PATH_TO_PERF_FROM_EXEC, PATH_TO_DATA_FROM_EXEC

options = {}

#Domaine
options["--ex, --ey, --ez"] = [100], [100], [100]
options["--timemax"] = [1.536]

#Activation PPM
options["--save-ppm"] = [""]
options["--ppm-plane"] = ["xy"]
options["--ppm-colormap"] = ["coolwarm"]
options["--ppm-slice-index"] = [0]

#Intervalle PPM (fréquence d’écriture)
options["--ppm-interval"] = [2**i for i in range(3, 9)]  # 8,16,32,64,128,256

#Dossiers
options["--ppm-folder"] = [PATH_TO_DATA_FROM_EXEC + "ppm/"]
options["--perf-file"] = [PATH_TO_PERF_FROM_EXEC + "perf_ppm.csv"]

var_env = {}
nbruns = 1

execute(
    "./semproxy",
    var_env,
    options,
    nbruns=nbruns,
    verbose=True,
    execPath=PATH_TO_EXEC
)