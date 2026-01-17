#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
"""
## Desctriptive stats

The script `desc.py` generates descriptives stats from snapshot or slice files

inputs:
- snapshot or slice file

outputs:
- descriptive stats file


# File Format Documentation

## Snapshot File Format
CSV file with the following columns:
- `x`: X-coordinate of the data point
- `y`: Y-coordinate of the data point
- `z`: Z-coordinate of the data point
- `p`: Pressure value at the data point

## Descriptive stats File Format
CSV file with the following columns:
        --TBD
    
- `left`: Left edge of the histogram bin
- `right`: Right edge of the histogram bin
- `count`: Number of data points in the histogram bin
"""

import pandas as pd
import argparse
import argcomplete
import numpy as np
import time


def descriptive_stats(input_file, bins=100, output_file=None, verbose=False):
    """
    if output_file is not None, save the descriptive stats to output_file

    :return: time taken, dataframe of the descriptive files
    """

    start_time = time.time()

    ######## load data
    data = pd.read_csv(input_file)

    if 'p' not in data.columns:
        raise ValueError(f"'p' column not found in data.")

    ######## compute histogram
    ##counts, bin_edges = np.histogram(data['p'], bins=bins)
    
    _mean = [np.mean(data['p'])]
    _min = [min(data['p'])]
    _max = [max(data['p'])]
    _variance = [np.var(data['p'])]
    _mediane = [np.median(data['p'])]
    desc_df = pd.DataFrame({
        'mean': _mean,
        'min': _min,
        'max': _max,
        'variance': _variance,
        'mediane': _mediane
    })

    if output_file is not None:
        desc_df.to_csv(output_file, index=False)
        if verbose:
            print(f"Descriptive stats saved to {output_file}")

    end_time = time.time()
    if verbose:
        print(f"Descriptive stats computed in {end_time - start_time:.4f} seconds.")



    return end_time - start_time, desc_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True)
    parser.add_argument('--bins', type=int, default=100)
    parser.add_argument('--output', type=str, required=False)
    parser.add_argument('--verbose', action='store_true', default=True)
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    descriptive_stats(args.input, bins=args.bins, output_file=args.output, verbose=args.verbose)







































































