# Cell counting pipeline for brain organoids
---
This repository contains a pipeline for semi-automatic cell and cell positive
counting in multi-channel z-stacks.

The pipeline is based on cell segmentation with pretrained deep learning networks
 [CellPose](http://www.cellpose.org/) and supports images from Zeiss microscopes (.czi, .lsm).


 ## Steps
 ### 1. Extract counting regions
 * Define region of interests in ImageJ/Fiji and add to *ROIManager*
 * Use Fiji/jython script in [imagej](./imagej/README.md) to extract cropped counting regions
 * Select channel name for segmentation (e.g. *blue*)
 * Choose number of z-planes used as basis for later count (e.g. *+-1*)
 * Enter approximate cell diameter

 Die command to run CellPose on extracted regions of interest will automatically be copied to clipboard

 ### 2. Run CellPose
 * Open command shell
 * Paste CellPose command (`RIGHT-click` or `CTRL+V`)

### 3. CellPose filter

#### Example: single file: *<my_seg.npy>*
`$ cellpose_filter --min-area 100 --max-area 4000 --min-circularity 0.8 <my_seg.npy>`

#### Example: batch for folder containing *_seg.npy* files
`$ cellpose_filter --min-area 100 --max-area 4000 --min-circularity 0.8 <folder containing files>`

#### Full usage and options
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


#### Output:
In the same folder as input file
* Excel table with filtered results

---
## Installation
Anaconda Python distribution with Python >= 3.6 recommended

1. git clone this repository
2. `cd cellpose_filter`
3. `pip install -e .`