#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
"""
## Slices

The script `slice.py` extracts 2D slices from 3D simulation (entire snapshot)

inputs:
- snapshot file
- axis (x, y, z)
- axis value

outputs:
- 2D slice file

# File Format Documentation

## Snapshot File Format
CSV file with the following columns:
- `x`: X-coordinate of the data point
- `y`: Y-coordinate of the data point
- `z`: Z-coordinate of the data point
- `p`: Pressure value at the data point

## 2D Slice Snapshot File Format
CSV file with the following columns:
- `x|y`: X-coordinate or Y-coordinate of the data point in the slice
- `y|z`: Y-coordinate or Z-coordinate of the data point in the slice
- `p`: Pressure value at the data point in the slice
"""

import pandas as pd
import argparse
import argcomplete
import time


def slice_snapshot(input_file, axis, value, output_file = None, verbose=False):
    """
    if output_file is not None, save the 2D slice to output_file

    :return: time taken, dataframe of the 2D slice
    """

    start_time = time.time()

    ######## load snapshot
    snapshot = pd.read_csv(input_file)

    ######## extract slice
    if axis not in snapshot.columns:
        raise ValueError(f"Axis '{axis}' not found in snapshot data.")

    slice_data = snapshot[snapshot[axis] == value]
    if slice_data.empty:
        raise ValueError(f"No value found for {axis} = {value}.")

    # Select relevant columns for 2D slice
    if axis == 'x':
        slice_data = slice_data[['y', 'z', 'p']]
    elif axis == 'y':
        slice_data = slice_data[['x', 'z', 'p']]
    elif axis == 'z':
        slice_data = slice_data[['x', 'y', 'p']]

    ######## save slice
    if output_file is not None:
        slice_data.to_csv(output_file, index=False)
        if verbose:
            print(f"2D slice saved to {output_file}")

    end_time = time.time()
    if verbose:
        print(f"Slice extraction completed in {end_time - start_time:.4f} seconds.")

    return end_time - start_time, slice_data







if __name__ == "__main__":
    ########## parse

    parser = argparse.ArgumentParser(description="Extract 2D slices from 3D simulation snapshot.")
    parser.add_argument('--input', type=str, required=True, help="Path to the input snapshot CSV file")
    parser.add_argument('--output', type=str, required=True, help="Path to save the 2D slice CSV file")
    parser.add_argument('--axis', type=str, choices=['x', 'y', 'z'], required=True, help="Axis along which to take the slice (x, y, or z)")
    parser.add_argument('--value', type=float, required=True, help="Value along the specified axis to take the slice")

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    slice_snapshot(args.input, args.axis, args.value, args.output, verbose=True)


















































