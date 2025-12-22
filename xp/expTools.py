import os
from itertools import *
import subprocess
import time


def iterateur_option(dicopt, sep=' '):
    options = []
    for opt, listval in dicopt.items():

        #multiple keys
        if "," in opt:
            opt = opt.replace(' ', '').split(',')

            # check listval is list of lists
            if not isinstance(listval, (list, tuple) ) or not all(isinstance(i, (list, tuple)) for i in listval):
                print(f"Error: {str(opt)}, {len(opt)} keys given, then list of lists expected.", flush=True)
                exit(1)

            # check number of lists matches number of keys
            if len(listval) != len(opt):
                print(f"Error: {str(opt)}, {len(opt)} key given, then {len(opt)} listval expected, but {len(listval)} given.", flush=True)
                exit(1)

            # check all lists have same length
            if not all(len(i) == len(listval[0]) for i in listval):
                print(f"Error: {str(opt)}, all lists must have the same length.", flush=True)
                exit(1)

            optlist = []
            for i in range(len(listval[0])):
                v = ""
                for j in range(len(opt)):
                    v += opt[j] + sep + str(listval[j][i]) + " "
                # remove last space
                v = v[:-1]
                optlist.append(v)
            options.append(optlist)



        # single key
        else:
            if not isinstance(listval, (list, tuple)):
                print(f"Error: {opt}, listval type (list or tuple) expected, but {type(listval)} given.", flush=True)
                exit(1)

            optlist = []
            for val in listval:
                optlist.append(opt + sep + str(val))
            options.append(optlist)

    return [' '.join(value) for value in product(*options)]


def execute(commande, var_env, option, nbruns=1, verbose=True, execPath='.'):
    path = os.getcwd()
    os.chdir(execPath)

    var_env_it = iterateur_option(var_env, "=")
    opt_it = iterateur_option(option)

    total_xp = nbruns * len(list(var_env_it)) * len(list(opt_it))
    current_xp = 1

    start_time = time.time()

    runs_average = -1

    for i in range(nbruns):
        for var in var_env_it:
            for opt in opt_it:

                if (verbose):
                    print("\n#######################\n", flush=True)
                    print("[" + str(current_xp) + "/" + str(total_xp) + "]", end=" ", flush=True)

                    elapsed_time = time.time() - start_time
                    print("Elapsed time: " + str(int(elapsed_time / 60)) + "m" + str(int(elapsed_time % 60)) + "s", flush=True)

                    average_time = elapsed_time / current_xp if runs_average == -1 else runs_average
                    estimated_time = average_time * ((total_xp - current_xp) + 1)
                    print("Estimated time remaining: " + str(int(estimated_time / 60)) + "m" + str(int(estimated_time % 60)) + "s", flush=True)

                    print(var + " " + commande + " " + opt + "\n", flush=True)

                if (subprocess.call([var + " " + commande + " " + opt], shell=True) == 1):
                    os.chdir(path)
                    return ("Error on the command used")

                current_xp += 1

        elapsed_time = time.time() - start_time
        runs_average = elapsed_time / current_xp

    print("\nExperiences done", flush=True)
    elapsed_time = time.time() - start_time
    print("Total elapsed time: " + str(int(elapsed_time / 60)) + "m" + str(int(elapsed_time % 60)) + "s", flush=True)

    os.chdir(path)
