#!/usr/bin/env python3

from expTools import execute
from common import PATH_TO_EXEC, PATH_TO_PERF_FROM_EXEC, PATH_TO_DATA_FROM_EXEC

ex = [10, 20, 50, 100, 150, 200]
ey = ex.copy()
ez = ex.copy()

snap_folders = [PATH_TO_DATA_FROM_EXEC + "snap_ad-hoc_slice_" + str(i) + "/" for i in ex]


options = {}
options["--ex, --ey, --ez, --snap-folder"] = ex, ey, ez, snap_folders
options["--timemax"] = [1.536]
options["--save-snapshot"] = [""]
options["--snap-interval"] = [256]
options["--perf-file"] = [PATH_TO_PERF_FROM_EXEC + "perf_ad-hoc_slice.csv"]


var_env = {}

nbruns = 1

execute(f'./semproxy', var_env, options, nbruns=nbruns, verbose=True, execPath=PATH_TO_EXEC)


######################  Ad-hoc slice extraction from snapshots ######################

import os
import sys
import glob
import time
PATH_TO_AD_HOC_SCRIPTS = "../ad-hoc_scripts/"
sys.path.append(os.path.abspath(PATH_TO_AD_HOC_SCRIPTS))

from slice import slice_snapshot

i=0
for snap_folder in snap_folders:

    #print(snap_folder)

    snap_folder = snap_folder.replace(PATH_TO_DATA_FROM_EXEC, "data/") # path from here /xp/ (not from exec)

    slice_folder = str(snap_folder).replace("snap_ad-hoc_slice", "slice_ad-hoc_slice")


    #print(snap_folder)
    #print(slice_folder)

    if not os.path.exists(slice_folder):
        os.makedirs(slice_folder)

    start_time = time.time()

    snapshots = glob.glob( os.path.join(snap_folder, "snapshot_*.csv") )

    for snap_file in snapshots:
        slice_filename = os.path.join(slice_folder, os.path.basename(snap_file).replace("snapshot", "slice_x0"))
        slice_snapshot(snap_file, axis='x', value=0, output_file = slice_filename, verbose=False)

    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"\nAd-hoc ex {ex[i]} slice extraction completed in {elapsed_time:.4f} seconds.\n")

    ## save performance
    perf_slice_file = "perf/perf_ad-hoc_slice_extraction.csv"

    f = open(perf_slice_file, "a")

    if f.tell() == 0:
        f.write("ex,ey,ez,nb_snapshots,slice_extraction_time\n")
        f.flush()

    f.write(f"{ex[i]},{ey[i]},{ez[i]},{len(snapshots)},{elapsed_time:.4f}\n")
    f.flush()
    i+=1









