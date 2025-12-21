#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
"""
## Histogram

The script `histogram.py` generates histograms from snapshot or slice files

inputs:
- snapshot or slice file
- number of bins

outputs:
- histogram file


# File Format Documentation

## Snapshot File Format
CSV file with the following columns:
- `x`: X-coordinate of the data point
- `y`: Y-coordinate of the data point
- `z`: Z-coordinate of the data point
- `p`: Pressure value at the data point

## Histogram File Format
CSV file with the following columns:
- `left`: Left edge of the histogram bin
- `right`: Right edge of the histogram bin
- `count`: Number of data points in the histogram bin
"""

import pandas as pd
import argparse
import argcomplete
import numpy as np
import time


def histogram(input_file, bins=100, output_file=None, verbose=False):
    """
    if output_file is not None, save the histogram to output_file

    :return: time taken, dataframe of the histogram
    """

    start_time = time.time()

    ######## load data
    data = pd.read_csv(input_file)

    if 'p' not in data.columns:
        raise ValueError(f"'p' column not found in data.")

    ######## compute histogram
    counts, bin_edges = np.histogram(data['p'], bins=bins)

    histogram_df = pd.DataFrame({
        'left': bin_edges[:-1],
        'right': bin_edges[1:],
        'count': counts
    })

    if output_file is not None:
        histogram_df.to_csv(output_file, index=False)
        if verbose:
            print(f"Histogram saved to {output_file}")

    end_time = time.time()
    if verbose:
        print(f"Histogram computed in {end_time - start_time:.4f} seconds.")



    return end_time - start_time, histogram_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True)
    parser.add_argument('--bins', type=int, default=100)
    parser.add_argument('--output', type=str, required=False)
    parser.add_argument('--verbose', action='store_true', default=True)
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    histogram(args.input, bins=args.bins, output_file=args.output, verbose=args.verbose)







































































