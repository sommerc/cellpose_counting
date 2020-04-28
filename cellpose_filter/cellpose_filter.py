import os
import glob
import numpy
import pandas
import argparse
import warnings
#import read_roi
import tifffile

from skimage import measure, segmentation, morphology

import matplotlib
from matplotlib import pyplot as plt

description = """Filter cellpose segmentations and create filtered tif and csv output of the segmentation"""

def seg_show(img, seg, seg2, ax=None, alpha=0.3):
    border  = segmentation.find_boundaries(seg, mode='inner')
    border2 = segmentation.find_boundaries(seg2, mode='inner')

    vmax = max(seg.max(), seg2.max())

    cmap = numpy.random.rand(vmax, 4)
    cmap[0, :] = [0,0,0, 1]
    cmap[:, 3] = alpha
    cmap = matplotlib.colors.ListedColormap(cmap)

    bcmap = numpy.random.rand(2, 4)
    bcmap[0, :] = [0,0,0, 0]
    bcmap[1, :] = [1,1,1, 0.5]

    bcmap = matplotlib.colors.ListedColormap(bcmap)


    f, ax = plt.subplots(1,2, sharex=True, sharey=True, figsize=(20,10))


    ax[0].imshow(img, 'gray')
    ax[1].imshow(img, 'gray')
    ax[0].imshow(seg, cmap=cmap, vmin=0, vmax=vmax)
    ax[1].imshow(seg2, cmap=cmap, vmin=0, vmax=vmax)
    ax[0].imshow(border, cmap=bcmap)
    ax[1].imshow(border2, cmap=bcmap)
    for a in ax: a.set_axis_off()
    ax[0].set_title("CellPose original")
    ax[1].set_title("Filtered")
    plt.tight_layout()

def apply_filter(seg, img, min_area, max_area, min_circularity):
    rp = measure.regionprops(seg, intensity_image=img)
    for r in rp:
        r.circularity = 4*numpy.pi*(r.area/r.perimeter**2)

    seg = seg.copy()
    for r in rp:
        remove = False
        if r.area < min_area: remove=True
        if r.area > max_area: remove=True
        if r.circularity < min_circularity: remove=True

        if remove: seg[r.coords[:, 0], r.coords[:, 1]] = 0

    return seg

def get_args():
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input',  type=str, nargs="+", help="Files or directory to process")
    parser.add_argument('--min-area', dest="min_area", type=int, help="Minimum area in pixels", default=1)
    parser.add_argument('--max-area', dest="max_area", type=int, help="Maximum area in pixels", default=999999)
    parser.add_argument('--min-circularity', dest="min_circ", type=float, help="Minimum circlularity (0-1;1 being perfect circle) ", default=0)

    return parser.parse_args()

def run_file(fn, args, display):
    print(f" * Process {fn}")
    args = get_args()

    npy = numpy.load(fn, allow_pickle=True)[()]
    img = npy["img"]
    seg = npy["masks"]
    seg = measure.label(seg)


    seg_new = apply_filter(seg, img, min_area=args.min_area, max_area=args.max_area, min_circularity=args.min_circ)

    if display:
        seg_show(img, seg, seg_new)

    seg_new = measure.label(seg_new).astype(numpy.uint16)
    rp_tab = measure.regionprops_table(seg_new, intensity_image=img, properties=['label', 'area', 'mean_intensity', 'centroid', ], cache=True, separator='-')
    tab = pandas.DataFrame(rp_tab)


    fn_base = os.path.splitext(fn)[0]
    print(fn, fn_base)

    tifffile.imsave(fn_base + "_flt.tif", seg_new, imagej=True)

    with open(fn_base + "_flt.tab", 'w', newline='') as f:
        f.write(f'# Count: {len(tab)}\n')
        f.write(f'#  - Min-area: {args.min_area}\n')
        f.write(f'#  - Max-area: {args.max_area}\n')
        f.write(f'#  - Min-circularity: {args.min_circ}\n')
        f.write("\n")

        tab.to_csv(f, sep="\t")

    print(f" -> Done")

def main():
    warnings.simplefilter("ignore")

    args = get_args()
    print("Running CellPose filter...")

    if os.path.isfile(args.input[0]):
        run_file(args.input[0], args, display=True)
        plt.show()
    elif os.path.isdir(args.input[0]):
        fn_list = glob.glob(os.path.join(args.input[0], "*seg.npy"))
        for fn in fn_list:
            run_file(fn, args, display=False)

    print("Finished")






if __name__ == "__main__":
    args = main()


