import os
import glob
import numpy
import pandas
import argparse
import warnings
#import read_roi
import tifffile

from skimage import measure, segmentation, morphology
from scipy.ndimage.morphology import binary_fill_holes

import matplotlib
from matplotlib import pyplot as plt

from mpldatacursor import datacursor

description = """Filter cellpose segmentations and create filtered tif and csv output of the segmentation"""

# uugly
to_delete = []

def seg_show(img, rp, seg, seg2, ax=None, alpha=0.3):
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

    class MyFormatter(object):
        def __init__(self, seg, rp):
            self.seg = seg
            self.rp = rp

        def __call__(self, x, y, **kwargs):
            try:
                label = self.seg[int(y+0.5), int(x+0.5)]
            except:
                return

            if label == 0 or (label > len(self.rp)-1):
                return

            res = ""
            for feat in ["area", "circularity", "roundness", "solidity"]:
                res += f"{feat}: {getattr(self.rp[label-1], feat):0.3f}\n"

            return res

    f, ax = plt.subplots(1,2, sharex=True, sharey=True, figsize=(20,10))

    ax[0].imshow(img, 'gray')
    ax[1].imshow(img, 'gray')
    ax[0].imshow(seg, cmap=cmap, vmin=0, vmax=vmax)
    ax[1].imshow(seg2, cmap=cmap, vmin=0, vmax=vmax)
    ax[0].imshow(border, cmap=bcmap)
    ax[1].imshow(border2, cmap=bcmap)
    for a in ax: a.set_axis_off()
    ax[0].set_title("CellPose original")
    ax[1].set_title("Filtered (double click on cell to delete)")
    plt.tight_layout()

    datacursor(artists=ax[0], formatter=MyFormatter(seg, rp), hover=True)

    def onclick(event):

        if event.xdata != None and event.ydata != None and event.inaxes == ax[1] and event.dblclick:
            label_to_delete = seg2[int(event.ydata+0.5), int(event.xdata+0.5)]
            if label_to_delete > 0:
                global to_delete
                to_delete.append(label_to_delete)
                ax[1].plot([event.xdata], [event.ydata], "rx", linewidth=2)

    cid = f.canvas.mpl_connect('button_press_event', onclick)


def apply_filter(seg, img, min_area, max_area, min_circularity, min_roundness, min_solidity):
    rp = measure.regionprops(seg, intensity_image=img)
    for r in rp:
        r.circularity = numpy.clip(4*numpy.pi*(r.area/r.perimeter**2), 0, 1)
        r.roundness   = numpy.clip(4 * (r.area) / (numpy.pi * r.major_axis_length**2), 0, 1)

    seg = seg.copy()
    for r in rp:
        remove = False
        if r.area < min_area: remove=True
        if r.area > max_area: remove=True
        if r.circularity < min_circularity: remove=True
        if r.roundness < min_roundness: remove=True
        if r.solidity < min_solidity: remove=True

        if remove: seg[r.coords[:, 0], r.coords[:, 1]] = 0



    return seg, rp

def get_args():
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input',  type=str, nargs="+", help="Files or directory to process")
    parser.add_argument('--min-area', dest="min_area", type=int, help="Minimum area in pixels", default=1)
    parser.add_argument('--max-area', dest="max_area", type=int, help="Maximum area in pixels", default=9999999)
    parser.add_argument('--min-circularity', dest="min_circ", type=float, help="Minimum circlularity (0-1;1 being perfect circle) ", default=0)
    parser.add_argument('--min-roundness', dest="min_round", type=float, help="Minimum roundness (0-1;1 being perfect circle) ", default=0)
    parser.add_argument('--min-solidity', dest="min_solid", type=float, help="Minimum solidity (0-1; ratio of region and its convex hull ", default=0)

    return parser.parse_args()

def run_file(fn, args, display):
    print(f" * Process {fn}")
    args = get_args()

    fn_base = os.path.splitext(fn)[0]

    npy = numpy.load(fn, allow_pickle=True)[()]
    img = npy["img"]
    seg = npy["masks"]

    # relabel
    seg = measure.label(seg, background=-1)

    # remove holes
    rp = measure.regionprops(seg)
    for r in rp:
        if r.euler_number < 1:
            seg[r.coords[:, 0], r.coords[:, 1]] = 0
            no_holes = binary_fill_holes(r.image)
            y_coords, x_coords = numpy.nonzero(no_holes)
            seg[r.bbox[0] + y_coords, r.bbox[1] + x_coords] = r.label

    seg_new, rp_tab = apply_filter(seg, img, min_area=args.min_area,
                                     max_area=args.max_area,
                                     min_circularity=args.min_circ,
                                     min_roundness=args.min_round,
                                     min_solidity=args.min_solid
                                     )
    if display:
        seg_show(img, rp_tab, seg, seg_new)
        plt.show()

    global to_delete
    print(f"  -- #{len(to_delete)} manual deletions")
    seg_new[numpy.isin(seg_new, to_delete)] = 0

    seg_new = measure.label(seg_new)

    rp_tab = measure.regionprops_table(seg_new, intensity_image=img, properties=['label', 'area', 'mean_intensity', 'centroid', ], cache=True, separator='-')
    tab = pandas.DataFrame(rp_tab)

    if fn_base.endswith("_seg"):
        fn_base = fn_base[:-4]

    other_tabs = []
    for other_img_fn in glob.glob(f"{fn_base}*.tif"):
        col_name = other_img_fn[len(fn_base)+1:-4]
        print("  -- reading intensity from", col_name)
        other_img = tifffile.imread(other_img_fn)

        other_tab = measure.regionprops_table(seg_new, intensity_image=other_img, properties=['mean_intensity' ], cache=True, separator='-')
        other_vals = other_tab['mean_intensity']

        tab[col_name] = other_vals




    with open(fn_base + "_filtered.tab", 'w', newline='') as f:
        f.write(f'# Count: {len(tab)}\n')
        f.write(f'#  - Min-area: {args.min_area}\n')
        f.write(f'#  - Max-area: {args.max_area}\n')
        f.write(f'#  - Min-circularity: {args.min_circ}\n')
        f.write(f'#  - Min-roundness: {args.min_round}\n')
        f.write(f'#  - Min-solidity: {args.min_solid}\n')
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
    else:
        print("  -- File or directory does not exist... aborting")

    print("Finished")






if __name__ == "__main__":
    main()


