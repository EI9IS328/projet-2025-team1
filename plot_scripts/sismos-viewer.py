#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import pandas as pd
import argparse
import argcomplete
import matplotlib.pyplot as plt




def plot_sismos(sismos, output_file=None):
    """
    Plot the seismogram data.

    :param sismos: dataframe of the seismogram
    :param output_file: path to save the output plot image (optional)
    """

    # draw the seismogram
    plt.figure(figsize=(10, 6))
    plt.plot(sismos['time'], sismos['p'], linestyle='-', color='blue')
    plt.xlabel('Time')
    plt.ylabel('Pressure (p)')
    plt.title('Seismogram')
    plt.grid(True)

    if output_file is not None:
        plt.savefig(output_file)
    else:
        plt.show()




if __name__ == "__main__":
    ########## parse

    parser = argparse.ArgumentParser(description="Plot the seismogram data.")
    parser.add_argument('--input', type=str, required=True, help="Path to the input sismos CSV file")
    parser.add_argument('--output', type=str, required=False, help="Path to save the output plot image")

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    sismos = pd.read_csv(args.input)

    # check that sismos has columns time, p
    if not ('time' in sismos.columns and 'p' in sismos.columns):
        raise ValueError("Input sismos CSV must have 'time', 'p' columns.")

    plot_sismos(sismos, args.output)
