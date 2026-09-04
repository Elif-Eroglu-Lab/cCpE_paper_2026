"""
Cell shape and orientation analysis.

1. Loads a segmented label image (ilastik multicolor output).
2. Relabels so every connected region gets a unique ID.
3. Measures for each cell: area, perimeter, orienatation, aspect ratio and other shape descriptors.
4. Save a spreadsheet (cell_metrics.xlsx) and plots, next to the image.

Requirements:
Python 3.11.15 (expected to work on >= 3.9)

Usage: 
python shape_morphometrics.py ccpe_R_no3.tif 0.375
python shape_morphometrics.py 0.375          #all .tif in this folder

Contact:
Tuba Majid / tuba.majid@ki.se
"""

import os, sys, logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from scipy.ndimage import binary_dilation
from skimage import io, measure, segmentation

# import warnings
# warnings.filterwarnings("ignore")

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "figure.facecolor": "white",      
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.transparent": False,
})

# --- Settings ----------------------------------------------------------------
PIXEL_SIZE_UM = 0.375        

SMALL_AREA_UM2  = 250.0     # to visualize all cells under this value (in um) for plots 03 and 04
SCALE_AREA_MAX  = 12000.0   # colour scale max for cell area plot
CC_CONNECTIVITY = 1         # 1 = orthogonal neighbours only when relabelling
DPI             = 600

_AREA_COLORS = ['#0077BB', '#EE7733', '#332288', '#CCBB44', '#CC3311', '#33BBEE',
                '#882255', '#44BB99', '#EE3377', '#DABFFF', '#FFDD44', '#AA4499']
CMAP_AREA = mcolors.ListedColormap(_AREA_COLORS)

# Orientation-bar geometry 
ORIENT_LINE_SCALE = 0.65    # multiplies sqrt(area/pi) to get half-length
ORIENT_LINE_FLOOR = 10.0    # minimum half-length in px, for tiny cells
ORIENT_LINE_WIDTH = 2.0     # matplotlib linewidth
ORIENT_CMAP  = "twilight"
ORIENT_TICKS = [-90, -60, -30, 0, 30, 60, 90]   # colorbar ticks 
ROSE_BAR_COLOR = "#4A90D9"  # single colour for every bar in plot 06 orienattion_rose

# --- Loading and relabelling -------------------------------------------------
def load_label_image(path):
    """Read a 2D label image. Only takes channel 0 if the file is really RGB."""
    raw = io.imread(path)
    if raw.ndim == 3 and raw.shape[-1] in (3, 4):
        raw = raw[..., 0]
    return raw.astype(np.int32)


def detect_background(raw):
    """Background = the most common value around the image border."""
    border = np.concatenate([raw[0, :], raw[-1, :], raw[:, 0], raw[:, -1]])
    vals, counts = np.unique(border, return_counts=True)
    return int(vals[np.argmax(counts)])


def relabel_connected_components(raw, background):
    """Give every connected group of same-valued pixels its own ID; background = 0."""
    return measure.label(raw, connectivity=CC_CONNECTIVITY,
                         background=background).astype(np.int32)


# --- Measuring ---------------------------------------------------------------
def gamma2_angle(mask, centroid):
    """
    Orientation angle calculated from the complex shape parameter gamma_2.

    Builds the 2x2 second-moment tensor of the contour VERTICES about the centroid, 
    then extracts the major-axis direction. Writing the tensor as a complex number, 
    gamma_2 = (Mcc - Mrr + 2i*Mcr) / (Mcc + Mrr), 
    its argument is twice the major-axis angle, hence the factor 0.5. 
    The doubling is what makes this a headless director: the result is defined modulo 180 deg, 
    not 360 deg, which is correct for an axis with no head or tail.

    Returns radians in image coordinates, i.e. with the row axis pointing down.
    """
    contours = measure.find_contours(mask.astype(float), 0.5)
    if not contours:
        return 0.0
    # Longest contour = outer boundary

    c = max(contours, key=len)
    dc = c[:, 1] - centroid[1]          # x offsets (columns)
    dr = c[:, 0] - centroid[0]          # y offsets (rows)
    Mcc, Mrr, Mcr = np.sum(dc**2), np.sum(dr**2), np.sum(dc * dr)
    if Mcc + Mrr <= 1e-10:
        return 0.0
    real = (Mcc - Mrr) / (Mcc + Mrr)
    imag = 2 * Mcr / (Mcc + Mrr)
    return float(0.5 * np.arctan2(imag, real))


def measure_cells(labels, pixel_size_um):
    """Measure every cell and return a dict keyed by cell ID."""
    cells = {}
    for p in measure.regionprops(labels):
        mask = np.pad(p.image, 1)
        centroid = (p.centroid_local[0] + 1, p.centroid_local[1] + 1)

        area_px   = p.area
        # Crofton perimeter: on rasterised ellipses at this cell size. Circularity and shape index are derived from it.
        perim_px  = measure.perimeter_crofton(mask)
        angle_rad = gamma2_angle(mask, centroid)

        cells[p.label] = dict(
            area_um2        = area_px * pixel_size_um**2,
            perimeter_um    = perim_px * pixel_size_um,
            feret_max_um    = p.feret_diameter_max * pixel_size_um,
            circularity     = 4 * np.pi * area_px / perim_px**2,
            shape_index     = perim_px / np.sqrt(area_px),
            eccentricity    = p.eccentricity,
            # NaN, not 1.0: a degenerate sliver is maximally elongated, not round.
            aspect_ratio    = (p.axis_major_length / p.axis_minor_length
                               if p.axis_minor_length > 1e-3 else np.nan),
            solidity        = p.solidity,
            centroid_row    = float(p.centroid[0]),
            centroid_col    = float(p.centroid[1]),
            centroid_x_um   = float(p.centroid[1]) * pixel_size_um,
            centroid_y_um   = float(p.centroid[0]) * pixel_size_um,
            angle_rad       = angle_rad,
            # Reported angle: degrees counter-clockwise from +x with y UP.
            # This is the mirror of angle_rad, which is in image coordinates,
            # so a bar leaning right in plot 05 has a POSITIVE angle here.
            orientation_deg = float(np.degrees(-angle_rad)),
            n_neighbours    = 0,
            free_edge       = False,
        )
    return cells


def add_neighbours_and_free_edge(labels, cells):
    """
    For each cell, count touching cells and check whether any part of its
    border faces background instead of another cell.

    A cell with a free edge has that part of its outline drawn by the
    segmentation rather than by contact with neighbours.
    """
    struct = np.ones((3, 3), dtype=bool)          # 8-connected ring
    for p in measure.regionprops(labels):
        r0, c0, r1, c1 = p.bbox
        r0, c0 = max(r0 - 1, 0), max(c0 - 1, 0)
        r1 = min(r1 + 1, labels.shape[0])
        c1 = min(c1 + 1, labels.shape[1])
        crop = labels[r0:r1, c0:c1]

        mask = crop == p.label
        ring = binary_dilation(mask, struct) & ~mask
        touching = crop[ring]

        cells[p.label]["n_neighbours"] = len({int(v) for v in np.unique(touching) if v != 0})
        cells[p.label]["free_edge"] = bool(np.any(touching == 0))


# --- Plot helpers ------------------------------------------------------------

def draw_boundaries(ax, labels, lw=0.4, color="black", alpha=0.5):
    bnd = segmentation.find_boundaries(labels, mode="thick")
    ax.contour(bnd, levels=[0.5], colors=color, linewidths=lw, alpha=alpha)


def cell_value_image(labels, cells, key):
    """An image where each cell is painted with one metric value."""
    img = np.full(labels.shape, np.nan, dtype=float)
    for cid, c in cells.items():
        img[labels == cid] = c[key]
    return img


def frame_axes(ax):
    ax.set_xticks([]); ax.set_yticks([])


def save(fig, output_dir, name):
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, name), dpi=DPI, bbox_inches="tight",
                facecolor="white", transparent=False)
    plt.close(fig)
    print(f"  saved {name}")


def orientation_colormap():
    """
    orientation_deg is a headless director in [-90, 90]: 
    0 = along +x, 
    +/-90 = along +y. 
    A cyclic colormap added so -90 and +90 receive the same colour, so the wraparound is
    invisible and two cells 2 degrees apart are always 2 degrees apart in colour 
    """
    cmap = matplotlib.colormaps[ORIENT_CMAP]
    norm = mcolors.Normalize(vmin=-90, vmax=90)
    return cmap, norm


def draw_orientation_bars(ax, cells, lw=ORIENT_LINE_WIDTH, color="black"):
    """Draw the black long-axis bar for every cell onto ax (image coords)."""
    for c in cells.values():
        r, col = c["centroid_row"], c["centroid_col"]
        L = max(np.sqrt(c["area_um2"] / np.pi) * ORIENT_LINE_SCALE, ORIENT_LINE_FLOOR)
        dr = np.sin(c["angle_rad"]) * L
        dc = np.cos(c["angle_rad"]) * L
        ax.plot([col - dc, col + dc], [r - dr, r + dr], color=color, lw=lw)


# --- Plots -------------------------------------------------------------------

def plot_area_map(labels, cells, output_dir):
    norm = mcolors.BoundaryNorm(np.linspace(0, SCALE_AREA_MAX, 13), 12)
    CMAP_AREA.set_bad("white")
    fig, ax = plt.subplots(figsize=(7, 7))
    im = ax.imshow(cell_value_image(labels, cells, "area_um2"),
                   cmap=CMAP_AREA, norm=norm, interpolation="nearest")
    draw_boundaries(ax, labels, color="white")
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Area (um^2)")
    frame_axes(ax)
    save(fig, output_dir, "01_cell_area_map.png")


def plot_area_and_id_labels(labels, cells, output_dir):
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    for ax, mode, title in [(axes[0], "area", "Cell area (um^2)"),
                            (axes[1], "id",   "Cell ID")]:
        draw_boundaries(ax, labels, lw=0.55, color="#333333", alpha=0.75)
        for cid, c in cells.items():
            txt = str(int(round(c["area_um2"]))) if mode == "area" else str(cid)
            ax.text(c["centroid_col"], c["centroid_row"], txt,
                    ha="center", va="center", fontsize=5, fontweight="bold")
        ax.set_title(title, fontweight="bold")
        ax.invert_yaxis()
        ax.set_aspect("equal")
        frame_axes(ax)
    save(fig, output_dir, "02_cell_area_labels.png")


def plot_small_cells_orange(labels, cells, output_dir):
    rgba = np.zeros((*labels.shape, 4))
    for cid, c in cells.items():
        if c["area_um2"] < SMALL_AREA_UM2:
            rgba[labels == cid] = mcolors.to_rgba("#EE7733")

    background = np.ones((*labels.shape, 3))
    background[labels == 0] = (0.93, 0.93, 0.93)   # background is always label 0

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(background)
    ax.imshow(rgba)
    draw_boundaries(ax, labels, lw=0.45, color="#888888", alpha=0.6)
    frame_axes(ax)
    save(fig, output_dir, "03_small_cells_orange.png")


def plot_small_cells_magma(labels, cells, output_dir):
    cmap, norm = plt.cm.magma, mcolors.Normalize(0, SMALL_AREA_UM2)
    rgba = np.zeros((*labels.shape, 4))
    for cid, c in cells.items():
        if c["area_um2"] < SMALL_AREA_UM2:
            rgba[labels == cid] = cmap(norm(c["area_um2"]))

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(np.ones((*labels.shape, 3)))
    ax.imshow(rgba)
    draw_boundaries(ax, labels, lw=0.45, color="#888888", alpha=0.6)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Area (um^2)")
    frame_axes(ax)
    save(fig, output_dir, "04_small_cells_magma.png")


def plot_orientation_field(labels, cells, output_dir):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(np.ones((*labels.shape, 3)))
    draw_boundaries(ax, labels, lw=0.5, color="#BBBBBB", alpha=0.9)
    draw_orientation_bars(ax, cells)
    frame_axes(ax)
    save(fig, output_dir, "05_orientation_field.png")

def plot_orientation_rose(cells, output_dir):
    angles = [c["orientation_deg"] for c in cells.values()]
    edges = np.linspace(-90, 90, 19)          # bin width fixed at 10 degrees
    counts, _ = np.histogram(angles, bins=edges)
    centers = np.deg2rad((edges[:-1] + edges[1:]) / 2)

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(4.8, 4.8))
    ax.bar(centers, counts, width=np.deg2rad(10),
           color=ROSE_BAR_COLOR, edgecolor="none", alpha=0.85)
    ax.set_theta_zero_location("E")
    ax.set_thetamin(-90); ax.set_thetamax(90)

    # Angular gridlines every 30 degrees only. Bin width stays at 10 degrees as set above) 
    theta_ticks_deg = list(range(-90, 91, 30))
    ax.set_thetagrids(theta_ticks_deg)
    ax.xaxis.set_major_locator(mticker.FixedLocator(np.deg2rad(theta_ticks_deg)))

    # Radial (count) axis every 2 cells, with rlim pinned 
    max_count = int(counts.max()) if len(counts) else 0
    rmax = max(2, int(np.ceil((max_count + 1) / 2.0) * 2))
    ax.set_rlim(0, rmax)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(2))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, pos: f"{int(round(v))}"))
    ax.tick_params(axis="y", labelsize=7)

    ax.set_title(f"n = {len(angles)} cells", fontsize=9)
    save(fig, output_dir, "06_orientation_rose.png")


# --- Spreadsheet -------------------------------------------------------------

def export_excel(cells, output_dir):
    """
    One row per cell. 
    Perimeter_um is the Crofton estimate (skimage.measure.perimeter_crofton).
    Centroid_X_um / Centroid_Y_um are image coordinates: X across, Y DOWN from the top-left corner.
    Orientation_deg is counter-clockwise from +x with y UP, so it is the mirror of what plot 05 draws.
    """
    df = pd.DataFrame([{
        "Cell_ID":         cid,
        "Area_um2":        c["area_um2"],
        "Perimeter_um":    c["perimeter_um"],
        "Feret_Max_um":    c["feret_max_um"],
        "Circularity":     c["circularity"],
        "Shape_Index_q":   c["shape_index"],
        "Eccentricity":    c["eccentricity"],
        "Aspect_Ratio":    c["aspect_ratio"],
        "Solidity":        c["solidity"],
        "Orientation_deg": c["orientation_deg"],
        "Centroid_X_um":   c["centroid_x_um"],
        "Centroid_Y_um":   c["centroid_y_um"],
        "N_Neighbours":    c["n_neighbours"],
        "Free_Edge":       c["free_edge"],
    } for cid, c in sorted(cells.items())])

    out = os.path.join(output_dir, "cell_metrics.xlsx")
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="Cell Metrics", index=False)
        ws = writer.sheets["Cell Metrics"]
        from openpyxl.utils import get_column_letter
        from openpyxl.styles import PatternFill, Font
        for i, col in enumerate(df.columns, 1):
            ws.column_dimensions[get_column_letter(i)].width = len(col) + 4
            for cell in ws[get_column_letter(i)][1:]:
                if isinstance(cell.value, float):
                    cell.number_format = "0.0000"
        for cell in ws[1]:
            cell.fill = PatternFill("solid", fgColor="1F4E79")
            cell.font = Font(color="FFFFFF", bold=True)

    print(f"  saved cell_metrics.xlsx  ({len(df)} cells, "
          f"{int(df.Free_Edge.sum())} with a free edge)")
    return df


# --- Run one image -----------------------------------------------------------

def analyse_image(image_path, pixel_size_um):
    image_path = Path(image_path)
    # Results go next to the image, not wherever the IDE happens to be pointing.
    output_dir = image_path.parent / f"output_{image_path.stem}"
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n{'='*60}\n  {image_path.name}\n{'='*60}")

    raw = load_label_image(image_path)
    bg = detect_background(raw)
    print(f"[1] Relabelling (background value = {bg})...")
    labels = relabel_connected_components(raw, bg)
    n_cells = int(labels.max())

    print("[2] Measuring cells...")
    cells = measure_cells(labels, pixel_size_um)
    add_neighbours_and_free_edge(labels, cells)

    print("[3] Writing spreadsheet...")
    export_excel(cells, output_dir)

    print("[4] Drawing plots...")
    plot_area_map(labels, cells, output_dir)
    plot_area_and_id_labels(labels, cells, output_dir)
    plot_small_cells_orange(labels, cells, output_dir)
    plot_small_cells_magma(labels, cells, output_dir)
    plot_orientation_field(labels, cells, output_dir)
    plot_orientation_rose(cells, output_dir)
    print(f"  done: {n_cells} cells")


# --- Start here --------------------------------------------------------------

if __name__ == "__main__":
    tif_files, pixel_size = [], PIXEL_SIZE_UM

    for arg in sys.argv[1:]:
        if arg.lower().endswith((".tif", ".tiff")):
            tif_files.append(Path(arg))
        else:
            pixel_size = float(arg)

    if pixel_size is None:
        raise SystemExit(
            "Pixel size is required.\n"
            "Set PIXEL_SIZE_UM near the top of this script, or run:\n"
            "  python aniso.py image.tif 0.375")

    if not tif_files:
        # No file named: process every .tif in the folder this script lives in.
        folder = Path(__file__).resolve().parent
        tif_files = sorted(p for p in folder.iterdir()
                           if p.suffix.lower() in (".tif", ".tiff"))
        print(f"Looking in {folder}")

    if not tif_files:
        raise SystemExit("No .tif or .tiff files found.")

    print(f"Cell shape analysis  |  {len(tif_files)} file(s)  |  "
          f"pixel size = {pixel_size} um")

    for path in tif_files:
        analyse_image(path, pixel_size)

    print(f"\nDone. {len(tif_files)} of {len(tif_files)} processed.")
