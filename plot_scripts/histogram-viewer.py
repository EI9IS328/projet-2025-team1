#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import pandas as pd
import argparse
import argcomplete
import matplotlib.pyplot as plt




def plot_histo(histo, output_file=None):
    """
    Plot the histogram.

    :param histo: dataframe of the histogram
    """

    # draw the histogram
    plt.bar((histo['left'] + histo['right']) / 2, histo['count'], width=(histo['right'] - histo['left']), align='center', edgecolor='black')
    plt.xlabel('Value')
    plt.ylabel('Count')
    plt.title('Histogram')

    if output_file is not None:
        plt.savefig(output_file)
    else :
        plt.show()














if __name__ == "__main__":
    ########## parse

    parser = argparse.ArgumentParser(description="Plot the histogram.")
    parser.add_argument('--input', type=str, required=True, help="Path to the input histogram CSV file")
    parser.add_argument('--output', type=str, required=False, help="Path to save the output plot image")

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    histo = pd.read_csv(args.input)

    # check that histo has columns left, right, count
    if not ('left' in histo.columns and 'right' in histo.columns and 'count' in histo.columns):
        raise ValueError("Input histogram CSV must have 'left', 'right', 'count' columns.")

    plot_histo(histo, args.output)













