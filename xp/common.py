
# check if the user, exec this script from here /xp
import os
if os.path.basename(os.getcwd()) != "xp":
    print("Please, execute this script from the 'stencil/xp' folder")
    exit(1)


PATH_TO_EXEC = "../build/bin"

PATH_TO_PERF_FROM_EXEC = "../../xp/perf/"
PATH_TO_DATA_FROM_EXEC = "../../xp/data/"