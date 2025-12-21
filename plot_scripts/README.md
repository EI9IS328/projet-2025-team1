# PLOT SCRIPTS

## Slice viewer

The script `slice-viewer.py` visualizes 2D slices extracted from 3D simulation snapshots.
### Inputs:
- 2D slice file (CSV format see FORMAT.md)

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

## Histogram viewer
The script `histogram-viewer.py` visualizes histograms from histogram file.
### Inputs:
- Histogram file (CSV format see FORMAT.md)

## Sismos viewer
The script `sismos-viewer.py` visualizes seismogram data from sismos file.
### Inputs:
- Sismos file (CSV format see FORMAT.md)








