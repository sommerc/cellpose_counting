## Cut and crop ImageJ macro
see [crop_rois_czi.py](./imagej/README.md)

## CellPose filter

## Install
Anaconda Python distribution with Python >= 3.6 recommended

```pip install -U git+https://git.ist.ac.at/csommer/cellpose_filter.git```

## Usage
From Anaconda prompt

* Get help information
```
$ cellpose_filter -h
```
will output:

```
usage: cellpose_filter [-h] [--no-display] [--min-area MIN_AREA]
                       [--max-area MAX_AREA] [--min-circularity MIN_CIRC]
                       [--min-roundness MIN_ROUND] [--min-solidity MIN_SOLID]
                       input [input ...]

Filter cellpose segmentations and create filtered tif and csv output of the
segmentation

positional arguments:
  input                 Files or directory to process

optional arguments:
  -h, --help            show this help message and exit
  --no-display
  --min-area MIN_AREA   Minimum area in pixels
  --max-area MAX_AREA   Maximum area in pixels
  --min-circularity MIN_CIRC
                        Minimum circlularity (0-1;1 being perfect circle)
  --min-roundness MIN_ROUND
                        Minimum roundness (0-1;1 being perfect circle)
  --min-solidity MIN_SOLID
                        Minimum solidity (0-1; ratio of region and its convex hull
```

## Example processing:

### Single file: *<my_seg.npy>*
`$ cellpose_filter --min-area 100 --max-area 4000 --min-circularity 0.8 <my_seg.npy>`

### Entire folder containing *_seg.npy* files
`$ cellpose_filter --min-area 100 --max-area 4000 --min-circularity 0.8 <folder containing files>`

## Output:
In the same folder as input file
* Excel table with filtered results

