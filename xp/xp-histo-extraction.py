#!/usr/bin/env python3

from expTools import execute
from common import PATH_TO_EXEC, PATH_TO_PERF_FROM_EXEC, PATH_TO_DATA_FROM_EXEC

ex = 10
ey = 10
ez = 10

options = {}
options["--ex, --ey, --ez"] = [ex], [ey], [ez]
options["--timemax"] = [1.536]
options["--save-snapshot"] = [""]
options["--snap-interval"] = [32]
options["--snap-folder"] = [PATH_TO_DATA_FROM_EXEC + "snap_ad-hoc_histo/"]
options["--perf-file"] = [PATH_TO_PERF_FROM_EXEC + "perf_ad-hoc_histo.csv"]


var_env = {}

nbruns = 1

execute(f'./semproxy', var_env, options, nbruns=nbruns, verbose=True, execPath=PATH_TO_EXEC)


######################  Ad-hoc histo extraction from snapshots ######################

import os
import sys
import glob
import time
PATH_TO_AD_HOC_SCRIPTS = "../ad-hoc_scripts/"
sys.path.append(os.path.abspath(PATH_TO_AD_HOC_SCRIPTS))

from histogram import histogram

snap_folder = "data/snap_ad-hoc_histo/"
histo_folder = "data/histo_ad-hoc_histo/"
if not os.path.exists(histo_folder):
    os.makedirs(histo_folder)

start_time = time.time()

snapshots = glob.glob( os.path.join(snap_folder, "snapshot_*.csv") )

for snap_file in snapshots:
    histo_filename = os.path.join(histo_folder, os.path.basename(snap_file).replace("snapshot", "histo"))
    histogram(snap_file, bins=100, output_file = histo_filename, verbose=False)

end_time = time.time()
elapsed_time = end_time - start_time

print(f"\nAd-hoc histo extraction completed in {elapsed_time:.4f} seconds.\n")

## save performance
perf_histo_file = "perf/perf_ad-hoc_histo_extraction.csv"

f = open(perf_histo_file, "a")

if f.tell() == 0:
    f.write("ex,ey,ez,nb_snapshots,histo_extraction_time\n")

f.write(f"{ex},{ey},{ez},{len(snapshots)},{elapsed_time:.4f}\n")









