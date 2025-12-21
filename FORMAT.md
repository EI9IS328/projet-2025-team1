# File Format Documentation

## Sismos File Format
CSV file with the following columns:
- `time`: Time step of the simulation
- `p`: Pressure value at the given time step

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

## Histogram File Format
CSV file with the following columns:
- `left`: Left edge of the histogram bin
- `right`: Right edge of the histogram bin
- `count`: Number of data points in the histogram bin

## Min/Max File Format
??




