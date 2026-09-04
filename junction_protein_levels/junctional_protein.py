"""
Quantification of junctional protein intensity from a binary junction mask and a matched fluorescence image.

Inputs:
SKELETON_PATH   
- junction mask, TIFF. Must be a closed network
PROTEIN_PATH    
- fluorescence image, TIFF.
BACKGROUND_INTENSITY
- constant subtracted from every junction mean.
                   
1.  Loads a binary junction mask and the corresponding IF image (Cldn6 or ZO1 stained
2. Calculates:
2.1 Per junction: mean intensity in a band of BAND_HALFWIDTH_PX around the line 
(nearest-junction assignment by Euclidean distance transform),
background-subtracted mean, band area, skeleton length, centroid.
2.2 Per cell: length-weighted mean of its junctions (raw and background-subtracted), 
and the unweighted CV across its junctions, on raw intensity.
4.  Saves junction_intensities.csv, cell_intensities.csv, 
three heatmaps anda QC overlay into a folder named after the mask file.

Requirements:
Python 3.11.15 (expected to work on >= 3.9)

Contact:
Tuba Majid / tuba.majid@ki.se
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import ndimage as ndi
from skimage.morphology import skeletonize

# ----------------------------------------------------------------------
SKELETON_PATH = "veh_R_no1.tif"        # binary junction mask
PROTEIN_PATH = "veh_R_no1_SV.tif"      # stained image

BACKGROUND_INTENSITY = 1.672

# ----------------------------------------------------------------------
# Measurement settings
# ----------------------------------------------------------------------
PROTEIN_CHANNEL = 0          # which channel to read if the stained image is RGB
BAND_HALFWIDTH_PX = 1        # band half-width used to MEASURE intensity
DISPLAY_HALFWIDTH_PX = 5     # band half-width used only in the figures
MIN_SEG_LEN_PX = 2           # ignore junction segments shorter than this
EXCLUDE_BORDER = True        # ignore cells touching the image frame
CV_DDOF = 1                  # 1 = sample SD (N-1), 0 = population SD (N)
MIN_JUNCTIONS_FOR_CV = 2     # a cell needs this many junctions to get a CV
DPI = 600

# ----------------------------------------------------------------------
# Figure settings
# ----------------------------------------------------------------------
CMAP_NAME = "magma"
BG_COLOR = "#E5E5E0"         # light grey; try "#EDE4D3" for light brown/tan

VMIN, VMAX = 0, 60          # intensity heatmaps
INT_STEP = 10                # colourbar tick spacing for intensity

CV_VMIN, CV_VMAX = 0, 1.0   # CV heatmap
CV_STEP = 0.2                # colourbar tick spacing for CV


# ======================================================================
# Loading
# ======================================================================
def load_skeleton(path):
    """Read the junction mask and thin it to lines exactly 1 pixel wide."""
    mask = tifffile.imread(path)
    if mask.ndim == 3:                    # colour image: any coloured pixel is a line
        mask = mask.any(axis=-1)
    return skeletonize(mask > 0)


def load_protein(path, channel=0):
    """Read the stained image as floating-point numbers."""
    image = tifffile.imread(path)
    if image.ndim == 3:
        image = image[..., channel]
    return image.astype(np.float32)


# ======================================================================
# Cutting the mask into junctions and cells
# ======================================================================
def split_skeleton_into_segments(skeleton, min_seg_len):
    """Cut the skeleton at its vertices. Each remaining piece is one junction.

    Returns an image where every junction has its own number (1, 2, 3, ...)
    and everything else is 0, plus how many junctions there are.
    """
    skel = skeleton.astype(np.uint8)

    # How many skeleton neighbours does each pixel have? (3x3 sum minus itself)
    neighbours = ndi.convolve(skel, np.ones((3, 3), np.uint8), mode="constant") - skel

    # 3 or more neighbours = a vertex where junctions meet. Cut there.
    vertices = skeleton & (neighbours >= 3)
    pieces, _ = ndi.label(skeleton & ~vertices, structure=np.ones((3, 3)))

    # Throw away pieces that are too short, then renumber 1, 2, 3, ...
    sizes = np.bincount(pieces.ravel())
    keep = sizes >= min_seg_len
    keep[0] = False                       # 0 is the background, never keep it

    renumber = np.zeros_like(sizes, dtype=np.int32)
    renumber[keep] = np.arange(1, keep.sum() + 1)
    return renumber[pieces], int(keep.sum())


def label_cells(skeleton, exclude_border=True):
    """Number the areas enclosed by the junction lines (the cells).

    The space outside the tissue is the biggest area that reaches the edge
    of the image; it is always dropped. With exclude_border=True, every
    cell touching the edge of the image is dropped as well.
    """
    # 4-connectivity here, so cells cannot leak into each other diagonally
    areas, n_areas = ndi.label(~skeleton, structure=ndi.generate_binary_structure(2, 1))

    frame = np.zeros(areas.shape, dtype=bool)
    frame[0, :] = frame[-1, :] = frame[:, 0] = frame[:, -1] = True
    at_frame = set(np.unique(areas[frame])) - {0}

    sizes = np.bincount(areas.ravel(), minlength=n_areas + 1)
    sizes[0] = 0
    outside = max(at_frame, key=lambda i: sizes[i])

    drop = at_frame if exclude_border else {outside}

    renumber = np.zeros(n_areas + 1, dtype=np.int32)
    new_id = 0
    for old_id in range(1, n_areas + 1):
        if old_id not in drop:
            new_id += 1
            renumber[old_id] = new_id
    return renumber[areas], new_id


def map_cells_to_junctions(junction_label, cell_label, n_cells):
    """For each cell, list the junctions sitting on its border.

    Returns a dictionary: cell number -> [junction numbers].
    """
    height, width = junction_label.shape
    rows, cols = np.where(junction_label > 0)
    junction_ids = junction_label[rows, cols]

    pairs = set()
    for d_row in (-1, 0, 1):                  # look at all 8 neighbours
        for d_col in (-1, 0, 1):
            if d_row == 0 and d_col == 0:
                continue
            r, c = rows + d_row, cols + d_col
            inside = (r >= 0) & (r < height) & (c >= 0) & (c < width)
            neighbour_cell = cell_label[r[inside], c[inside]]
            touches_a_cell = neighbour_cell > 0
            pairs.update(zip(junction_ids[inside][touches_a_cell].tolist(),
                             neighbour_cell[touches_a_cell].tolist()))

    cell_to_junctions = {cell: [] for cell in range(1, n_cells + 1)}
    for junction, cell in pairs:
        cell_to_junctions[cell].append(junction)
    return cell_to_junctions


# ======================================================================
# Measuring
# ======================================================================
def assign_territories(junction_label, band_halfwidth_px):
    """Widen each junction line into a band of the given half-width.

    Every pixel goes to its nearest junction line; pixels further away than
    band_halfwidth_px stay as 0 and are not measured.
    """
    distance, (nearest_row, nearest_col) = ndi.distance_transform_edt(
        junction_label == 0, return_indices=True)
    nearest_junction = junction_label[nearest_row, nearest_col]
    return np.where(distance <= band_halfwidth_px, nearest_junction, 0)


def per_junction_stats(territory, junction_label, protein, n_junctions):
    """Mean brightness in each band, plus band area, line length and centre."""
    ids = np.arange(1, n_junctions + 1)
    ones = np.ones_like(protein)
    return dict(
        mean_int=np.asarray(ndi.mean(protein, labels=territory, index=ids)),
        area_px=np.asarray(ndi.sum(ones, labels=territory, index=ids)),
        length_px=np.asarray(ndi.sum(ones, labels=junction_label, index=ids)),
        centroids=ndi.center_of_mass(ones, labels=junction_label, index=ids),
    )


def per_cell_stats(cell_to_junctions, junction_mean, junction_mean_bgsub,
                   junction_length, n_cells):
    """For each cell: the length-weighted mean of its junctions (raw and
    background-subtracted), and the coefficient of variation of its
    junctions (unweighted, on raw intensity)."""
    rows = []
    for cell in range(1, n_cells + 1):
        junctions = cell_to_junctions[cell]
        row = dict(cell_id=cell, n_junctions=len(junctions),
                   mean_intensity=np.nan, mean_intensity_bgsub=np.nan,
                   cv_unweighted=np.nan)

        if junctions:
            index = np.asarray(junctions) - 1      # junction 1 is at position 0
            lengths = junction_length[index]
            raw = junction_mean[index]
            bgsub = junction_mean_bgsub[index]

            row["mean_intensity"] = float(np.average(raw, weights=lengths))
            row["mean_intensity_bgsub"] = float(np.average(bgsub, weights=lengths))

            if len(junctions) >= MIN_JUNCTIONS_FOR_CV:
                row["cv_unweighted"] = float(np.std(raw, ddof=CV_DDOF) / np.mean(raw))

        rows.append(row)
    return pd.DataFrame(rows)


def build_junction_df(stats, n_junctions, background, mean_int_bgsub):
    """One row per junction."""
    return pd.DataFrame(dict(
        junction_id=np.arange(1, n_junctions + 1),
        length_px=stats["length_px"].astype(int),
        area_px=stats["area_px"].astype(int),
        mean_intensity=stats["mean_int"],
        background_intensity=background,
        mean_intensity_bgsub=mean_int_bgsub,
    ))


# ======================================================================
# Figures
# ======================================================================
def colorize_by_label(label_image, values, n_labels):
    """Paint each numbered region with its measured value; 0 becomes NaN."""
    lookup = np.full(n_labels + 1, np.nan, dtype=np.float32)
    lookup[1:] = values
    return np.ma.masked_invalid(lookup[label_image])

def make_colormap(bad_color=BG_COLOR):
    cmap = plt.get_cmap(CMAP_NAME).copy()
    cmap.set_bad(bad_color)
    return cmap


def add_colorbar(fig, ax, image, label, vmin, vmax, step):
    cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.05)
    bar = fig.colorbar(image, cax=cax)
    bar.set_label(label)
    bar.set_ticks(np.arange(vmin, vmax + step / 2, step))


def new_figure(size=10):
    fig, ax = plt.subplots(figsize=(size, size))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor(BG_COLOR)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    return fig, ax


def fig_junction_heatmap(display_territory, values, n_junctions, out_path,
                         label="Background-subtracted intensity (a.u.)"):
    fig, ax = new_figure()
    image = ax.imshow(colorize_by_label(display_territory, values, n_junctions),
                      cmap=make_colormap(), norm=Normalize(VMIN, VMAX),
                      interpolation="nearest")
    add_colorbar(fig, ax, image, label, VMIN, VMAX, INT_STEP)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor=BG_COLOR)
    plt.close(fig)


def fig_cell_heatmap(cell_label, values, skeleton, n_cells, out_path,
                     label, vmin, vmax, step):
    fig, ax = new_figure()
    image = ax.imshow(colorize_by_label(cell_label, values, n_cells),
                      cmap=make_colormap(), norm=Normalize(vmin, vmax),
                      alpha=0.9, interpolation="nearest")
    # draw the junction lines on top so the cell outlines stay visible
    lines = np.ma.masked_where(~skeleton, np.ones(skeleton.shape))
    ax.imshow(lines, cmap="gray_r", alpha=0.6)
    add_colorbar(fig, ax, image, label, vmin, vmax, step)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor=BG_COLOR)
    plt.close(fig)


def fig_junction_labeled(protein, junction_label, values, junction_centroids,
                         cell_centroids, n_junctions, n_cells, out_path):
    """Stained image with junctions coloured and every ID written on it."""
    thicker = ndi.grey_dilation(junction_label, size=(3, 3))   # easier to see

    fig, ax = plt.subplots(figsize=(14, 14))
    ax.imshow(protein, cmap="gray", vmax=float(np.percentile(protein, 99.5)))
    ax.imshow(colorize_by_label(thicker, values, n_junctions),
              cmap=make_colormap(bad_color=(0, 0, 0, 0)), norm=Normalize(VMIN, VMAX))

    for i, (row, col) in enumerate(cell_centroids):            # cell IDs
        ax.text(col, row, str(i + 1), color="black", fontsize=7,
                ha="center", va="center",
                bbox=dict(facecolor="yellow", alpha=0.6, pad=0.4, edgecolor="none"))
    for i, (row, col) in enumerate(junction_centroids):        # junction IDs
        ax.text(col, row, str(i + 1), color="white", fontsize=6,
                ha="center", va="center",
                bbox=dict(facecolor="black", alpha=0.6, pad=0.3, edgecolor="none"))

    ax.set_title(f"Junction IDs (n={n_junctions}) and cell IDs (n={n_cells})")
    ax.axis("off")
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


# ======================================================================
def main():
    out_dir = os.path.splitext(os.path.basename(SKELETON_PATH))[0]
    os.makedirs(out_dir, exist_ok=True)

    skeleton = load_skeleton(SKELETON_PATH)
    protein = load_protein(PROTEIN_PATH, PROTEIN_CHANNEL)

    # Cut the mask into junctions and cells, and link them together
    junction_label, n_junctions = split_skeleton_into_segments(skeleton, MIN_SEG_LEN_PX)
    cell_label, n_cells = label_cells(skeleton, EXCLUDE_BORDER)
    cell_to_junctions = map_cells_to_junctions(junction_label, cell_label, n_cells)
    print(f"Found {n_junctions} junctions and {n_cells} cells.")

    # Measure the stain in a band around each junction
    territory = assign_territories(junction_label, BAND_HALFWIDTH_PX)
    stats = per_junction_stats(territory, junction_label, protein, n_junctions)
    mean_int_bgsub = stats["mean_int"] - BACKGROUND_INTENSITY

    # Tables
    junction_df = build_junction_df(stats, n_junctions, BACKGROUND_INTENSITY,
                                    mean_int_bgsub)
    cell_df = per_cell_stats(cell_to_junctions, stats["mean_int"], mean_int_bgsub,
                             stats["length_px"], n_cells)
    junction_df.to_csv(os.path.join(out_dir, "junction_intensities.csv"), index=False)
    cell_df.to_csv(os.path.join(out_dir, "cell_intensities.csv"), index=False)

    # Figures
    display_territory = assign_territories(junction_label, DISPLAY_HALFWIDTH_PX)
    cell_centroids = ndi.center_of_mass(np.ones(cell_label.shape), labels=cell_label,
                                        index=np.arange(1, n_cells + 1))

    fig_junction_heatmap(display_territory, mean_int_bgsub, n_junctions,
                         os.path.join(out_dir, "junction_heatmap_bgsub.png"))

    fig_cell_heatmap(cell_label, cell_df["mean_intensity_bgsub"].values, skeleton,
                     n_cells, os.path.join(out_dir, "cell_heatmap_bgsub.png"),
                     label="Mean intensity per cell (a.u.)",
                     vmin=VMIN, vmax=VMAX, step=INT_STEP)

    fig_cell_heatmap(cell_label, cell_df["cv_unweighted"].values, skeleton,
                     n_cells, os.path.join(out_dir, "cv.png"),
                     label="CoV within-cell",
                     vmin=CV_VMIN, vmax=CV_VMAX, step=CV_STEP)

    fig_junction_labeled(protein, junction_label, mean_int_bgsub, stats["centroids"],
                         cell_centroids, n_junctions, n_cells,
                         os.path.join(out_dir, "QC.png"))

    print(f"Done. All outputs are in: {out_dir}")


if __name__ == "__main__":
    main()