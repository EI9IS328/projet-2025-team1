#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

"""
## Slices movie
The script `slices-movie.py` creates a movie from a series of 2D slice files.
The visual style is consistent with `slice-viewer.py` (matplotlib-based).
Use alpha-numeric sorting for input files to ensure correct frame order.
### Inputs:
- 2D slice files (CSV format see FORMAT.md)
- Output movie file (MP4 format)
- Optional: Frame rate (default: 30 fps)
### Outputs:
- Movie file (MP4 format)
### Example usage:
```bash
python slices-movie.py --input slices/slice_*.csv --output_file output_movie.mp4 --frame_rate 24
```
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import argcomplete
import glob
import re
import os
import numpy as np
import cv2




def create_movie(inputs, output_file, frame_rate):

    ## get format from first file (x,y,p or y,z,p or x,z,p)
    first_file = inputs[0]
    first_slice = pd.read_csv(first_file)
    dim1 = None
    dim2 = None
    if ('x' in first_slice.columns and 'y' in first_slice.columns and 'p' in first_slice.columns):
        dim1, dim2 = 'x', 'y'
    elif ('y' in first_slice.columns and 'z' in first_slice.columns and 'p' in first_slice.columns):
        dim1, dim2 = 'y', 'z'
    elif ('x' in first_slice.columns and 'z' in first_slice.columns and 'p' in first_slice.columns):
        dim1, dim2 = 'x', 'z'
    else:
        raise ValueError("Input slice CSV must have either 'x','y','p' or 'y','z','p' or 'x','z','p' columns.")

    print(f"Generating movie with {len(inputs)} frames at {frame_rate} fps...")

    fig = plt.figure(figsize=(10, 8))
    video_writer = None

    try:
        for i, input_file in enumerate(inputs):
            print(f"Processing frame {i+1}/{len(inputs)}: {input_file}")

            try:
                slice_df = pd.read_csv(input_file)
            except Exception as e:
                print(f"Error reading {input_file}: {e}")
                continue

            # Clear figure to reset layout and colorbar
            fig.clf()
            ax = fig.add_subplot(111)

            # Plot logic consistent with slice-viewer.py
            pivot = slice_df.pivot(index=dim2, columns=dim1, values='p')
            im = ax.imshow(pivot, origin='lower', aspect='auto',
                           extent=[pivot.columns.min(), pivot.columns.max(), pivot.index.min(), pivot.index.max()])

            fig.colorbar(im, ax=ax, label='Pressure (p)')
            ax.set_xlabel(dim1)
            ax.set_ylabel(dim2)
            ax.set_title(f'2D Slice - {os.path.basename(input_file)}')

            # Render to buffer
            fig.canvas.draw()

            # Get buffer (Matplotlib 3.9+ compatible)
            buf = fig.canvas.buffer_rgba()
            img_array = np.frombuffer(buf, dtype=np.uint8)
            w, h = fig.canvas.get_width_height()
            img_array = img_array.reshape((h, w, 4))

            # Convert RGBA to BGR for OpenCV
            img_bgr = cv2.cvtColor(img_array, cv2.COLOR_RGBA2BGR)

            # Initialize video writer on first frame
            if video_writer is None:
                height, width, _ = img_bgr.shape
                fourcc = cv2.VideoWriter_fourcc(*'mp4v')
                video_writer = cv2.VideoWriter(output_file, fourcc, frame_rate, (width, height))

            video_writer.write(img_bgr)

    except KeyboardInterrupt:
        print("\nInterrupted by user.")
    finally:
        if video_writer:
            video_writer.release()
            print(f"Movie saved to {output_file}")
        plt.close(fig)


if __name__ == "__main__":
    ########## parse

    parser = argparse.ArgumentParser(description="Create a movie from 2D slice files. Ex: python slices-movie.py --input slices/slice_*.csv --output_file output_movie.mp4 --frame_rate 24")
    parser.add_argument('--input', type=str, nargs='+', required=True, help="List of input 2D slice files (CSV format)")
    parser.add_argument('--output', type=str, required=True, help="Path to save the output movie file (MP4 format)")
    parser.add_argument('--frame_rate', type=int, default=30, help="Frame rate of the output movie (default: 30 fps)")
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # check fps is positive
    if args.frame_rate <= 0:
        raise ValueError("Frame rate must be a positive integer.")

    # check output file ends with .mp4
    if not args.output.endswith('.mp4'):
        raise ValueError("Output file must be in MP4 format (end with .mp4).")

    # check inputs number is at least 1
    if len(args.input) < 1:
        raise ValueError("At least one input file must be provided.")

    # sort input files alphanumerically
    input_files = sorted(args.input, key=lambda f: [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', f)])

    # check all input files exist + *.csv
    for input_file in input_files:
        if not os.path.isfile(input_file):
            raise FileNotFoundError(f"Input file {input_file} does not exist.")
        if not input_file.endswith('.csv'):
            raise ValueError(f"Input file {input_file} must be in CSV format (end with .csv).")


    create_movie(input_files, args.output, args.frame_rate)
