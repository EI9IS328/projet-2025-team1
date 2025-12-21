#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import pandas as pd
import argparse
import argcomplete
import matplotlib.pyplot as plt



def plot_slice(slice, output_file=None):
    """
    Plot the 2D slice.

    :param slice: dataframe of the 2D slice
    """

    # get what slice dimensions are
    dim1 = slice.columns[0]
    dim2 = slice.columns[1]

    # draw the slice like a heatmap without the value in cells
    # draw like a image (dim1.size x dim2.size) with color as p
    pivot = slice.pivot(index=dim2, columns=dim1, values='p')
    plt.imshow(pivot, origin='lower', aspect='auto', extent=[pivot.columns.min(), pivot.columns.max(), pivot.index.min(), pivot.index.max()])
    plt.colorbar(label='Pressure (p)')
    plt.xlabel(dim1)
    plt.ylabel(dim2)
    plt.title(f'2D Slice')

    if output_file is not None:
        plt.savefig(output_file)
    else :
        plt.show()




if __name__ == "__main__":
    ########## parse

    parser = argparse.ArgumentParser(description="Plot the 2D slice.")
    parser.add_argument('--input', type=str, required=True, help="Path to the input snapshot slice CSV file")
    parser.add_argument('--output', type=str, required=False, help="Path to save the output plot image")

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    slice = pd.read_csv(args.input)

    # check that slice has x,y,p or y,z,p or x,z,p
    if not (('x' in slice.columns and 'y' in slice.columns and 'p' in slice.columns) or
            ('y' in slice.columns and 'z' in slice.columns and 'p' in slice.columns) or
            ('x' in slice.columns and 'z' in slice.columns and 'p' in slice.columns)):
        raise ValueError("Input slice CSV must have either 'x','y','p' or 'y','z','p' or 'x','z','p' columns.")

    plot_slice(slice, args.output)

















