"""
make_skeletons.py — Generate 1px skeletons from ilastik Multicut label TIFs.
Outputs <stem>_skeleton.tif for every .tif/.tiff in the current directory.

Usage:
  python make_skeletons.py                 # all .tif in cwd
  python make_skeletons.py image.tif       # specific file(s)
"""

import numpy as np
from skimage import io, segmentation, morphology, measure
from scipy.ndimage import label as ndlabel, binary_dilation, convolve
from pathlib import Path
import warnings, sys, glob

warnings.filterwarnings('ignore')

# ── Parameters ──────────────────────────────────────────────────────────
MIN_CELL_AREA_PX = 50
MIN_REAL_CELL_NEIGHBOURS = 3
EXPAND_DISTANCE = 3       # Gap closing expansion in pixels
PRUNE_ITERATIONS = 30     # Skeleton endpoint pruning iterations


def make_skeleton(label_image):
    """Correct segmentation → gap close → boundaries → 1px skeleton."""

    corrected = label_image.copy().astype(np.int32)
    struct3 = np.ones((3, 3), dtype=bool)

    # Stage 1: Split disconnected components
    next_new = int(corrected.max()) + 1
    for prop in measure.regionprops(label_image):
        if prop.label == 0:
            continue
        mask = label_image == prop.label
        cc_map, n = ndlabel(mask)
        if n <= 1:
            continue
        cc_sizes = sorted([(int(np.sum(cc_map == i)), i)
                           for i in range(1, n + 1)], reverse=True)
        for rank, (sz, cc_id) in enumerate(cc_sizes):
            if rank == 0:
                continue
            corrected[cc_map == cc_id] = next_new
            next_new += 1

    # Stage 2: Remove slivers
    for prop in measure.regionprops(corrected):
        if prop.label == 0 or prop.area >= MIN_CELL_AREA_PX:
            continue
        mask = corrected == prop.label
        dil = binary_dilation(mask, struct3)
        nbrs = set(int(v) for v in np.unique(corrected[dil & ~mask])
                   if v not in (0, prop.label))
        if len(nbrs) < MIN_REAL_CELL_NEIGHBOURS:
            corrected[corrected == prop.label] = 0

    # Gap closing
    corrected = segmentation.expand_labels(corrected, distance=EXPAND_DISTANCE)

    # Create skeleton
    boundaries = segmentation.find_boundaries(corrected, mode='thick')
    skeleton = morphology.skeletonize(boundaries)

    # Prune dangling branches
    kernel = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])
    for _ in range(PRUNE_ITERATIONS):
        nc = convolve(skeleton.astype(int), kernel, mode='constant')
        endpoints = skeleton & (nc == 1)
        if not np.any(endpoints):
            break
        skeleton[endpoints] = False

    return skeleton


def process(image_path):
    stem = Path(image_path).stem
    img = io.imread(image_path)
    if img.ndim > 2:
        img = img[..., 0]
    label_img = img.astype(np.int32)

    skeleton = make_skeleton(label_img)

    out_path = f"{stem}_skeleton.tif"
    io.imsave(out_path, (skeleton.astype(np.uint8) * 255))
    print(f"{Path(image_path).name} -> {out_path}  ({int(skeleton.sum())} px)")


if __name__ == "__main__":
    tif_files = [a for a in sys.argv[1:] if a.endswith(('.tif', '.tiff'))]
    if not tif_files:
        tif_files = sorted(glob.glob("*.tif") + glob.glob("*.tiff"))
    if not tif_files:
        print("No .tif files found.")
        sys.exit(1)

    for path in tif_files:
        try:
            process(path)
        except Exception as e:
            print(f"ERROR ({path}): {e}")

    print("Done.")