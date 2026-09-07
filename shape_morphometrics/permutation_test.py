"""
Step 1  For each field of view (FOV), compute R, (R = 0, cells point every which way. R = 1, cells are all parallel.)
Step 2  For each FOV, compute a Rayleigh p-value.
Step 3  Call the FOV "aligned" (non-isotropic) if that p < 0.05.
Step 4  Compare the PROPORTION of aligned FOVs between conditions with a permutation test in which the ANIMAL is the unit that gets shuffled.

Restricted (whole-block) permutation
scipy.stats.permutation_test is handed arrays of ANIMAL NUMBERS, not FOV outcomes. 
permutation_type="independent" pools those numbers and re-deals them into groups of the original sizes. 

The statistic reaches a FOV only by looking it up through its animal, so an animal's FOVs can never be split apart.
n_resamples=np.inf SciPy tries every possible deal, so the pairwise p-values are exact rather than estimated.

NOTE: Orientation.xlsx contains one FOV per column, angles in degrees.

Requires: numpy, pandas, scipy, matplotlib, openpyxl
================================================================================
"""

import os
import re
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import permutation_test
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


INPUT_FILE = "Orientation.xlsx"
SHEET      = 0
OUTDIR     = "results"

ALPHA      = 0.05        # a FOV is "aligned" if its Rayleigh p is below this
SEED       = 1           # only used by the optional omnibus test; put it in Methods

CONDITIONS = ["naive", "vehicle", "cCpE"]          # left-to-right order in the plot
COLOURS    = {"naive": "#7f7f7f", "vehicle": "#4472c4", "cCpE": "#e0a138"}


# Read the file and work out which animal each FOV came from
def condition_of(column_name):
    """naive / vehicle / cCpE, from the column name."""
    name = column_name.lower()
    if "naive" in name:
        return "naive"
    if "ccpe" in name:
        return "cCpE"
    if "veh" in name:
        return "vehicle"
    return "UNKNOWN"


def animal_of(column_name):
    """Which animal a FOV belongs to. Two FOVs give the same answer when
    they are the left and right side of one animal, or a repeated field.
    """
    name = re.sub(r"^output_", "", column_name.lower())
    name = re.sub(r"_[lr]_", "_", name)          # side marker in the middle
    name = re.sub(r"(\d)[lr]$", r"\1", name)     # side marker at the end
    name = re.sub(r"(no\d+)_\d+$", r"\1", name)  # repeated field
    return name


def read_fov_columns(path, sheet):
    """Return {fov name: array of angles in degrees}, blanks dropped."""
    table = pd.read_excel(path, sheet_name=sheet)
    angles_by_fov = {}
    for column in table.columns:
        values = pd.to_numeric(table[column], errors="coerce").dropna().to_numpy()
        if len(values) > 0:
            angles_by_fov[str(column)] = values
    return angles_by_fov


# Circular statistics calculations per FOV
def order_parameter(angles_degrees):
    # R, Equivalent to MATLAB CircStat's circ_r(2*theta)
    radians = np.deg2rad(angles_degrees)
    return float(np.abs(np.mean(np.exp(2j * radians))))


def rayleigh_p_value(angles_degrees):
    # Rayleigh test of uniformity 
    n = len(angles_degrees)
    Z = n * order_parameter(angles_degrees) ** 2
    return float(np.exp(-Z))


# The permutation test ---------------------------------------------------------------------------------------------------------------------------
def proportion_aligned(animal_numbers, aligned_by_animal):
    #Proportion of aligned FOVs across these animals, pooling their FOVs.
    numbers = np.asarray(animal_numbers, dtype=int).ravel()
    return float(np.concatenate([aligned_by_animal[i] for i in numbers]).mean())


def pairwise_test(animals_in_group, aligned_by_animal, group_a, group_b):
    # Exact permutation test of group_a against group_b.
    a = np.asarray(animals_in_group[group_a], dtype=float)
    b = np.asarray(animals_in_group[group_b], dtype=float)

    def gap(animals_a, animals_b):
        #The size of the difference, ignoring its direction. Paired with alternative='greater'. This is a two-sided test.
        return abs(proportion_aligned(animals_a, aligned_by_animal)
                   - proportion_aligned(animals_b, aligned_by_animal))

    result = permutation_test(
        (a, b), gap,
        permutation_type="independent",   # re-deal the animals between groups
        vectorized=False,
        n_resamples=np.inf,               # try every deal -> exact p-value
        alternative="greater",
    )

    # The outcome is a proportion over 1-2 FOVs per animal, so many different deals give identical results. 
    null = np.asarray(result.null_distribution, dtype=float)
    floor = float((null >= null.max() - 1e-12).sum() / len(null))

    return dict(
        comparison=f"{group_a} vs {group_b}",
        n_animals_a=len(a),
        n_animals_b=len(b),
        difference=(proportion_aligned(a, aligned_by_animal)
                    - proportion_aligned(b, aligned_by_animal)),
        p_value=float(result.pvalue),
        smallest_possible_p=floor,
    )


def holm_correction(p_values):
    # Compare the corrected values against ALPHA 0.05
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    n = len(p_values)
    adjusted = np.empty(n)
    running_max = 0.0
    for rank, index in enumerate(order):
        running_max = max(running_max, (n - rank) * p_values[index])
        adjusted[index] = min(running_max, 1.0)
    return adjusted


# Figure
def make_figure_blocks(fov_table, animal_table, groups):
    #One block per FOV: filled = FOV aligned, hollow = if aligned.
    #bar spans 0 to 100, so one block is worth 100/n_fovs.

    fig, ax = plt.subplots(figsize=(4.4, 5.0))
    bar_width = 0.62

    for x, group in enumerate(groups):
        colour = COLOURS.get(group, "grey")
        fovs = fov_table[fov_table.condition == group]

        stack = sorted(fovs.aligned, reverse=True)   # aligned FOVs at the bottom

        slot = 100 / len(stack)          # what one FOV is worth, in percent
        left, right = x - bar_width / 2, x + bar_width / 2

        # the whole bar, then the aligned FOVs filled in on top of it
        ax.add_patch(plt.Rectangle((left, 0), bar_width, 100, facecolor="white",
                                   edgecolor=colour, linewidth=1.8, zorder=2))
        for i, aligned in enumerate(stack):
            if aligned:
                ax.add_patch(plt.Rectangle((left, i * slot), bar_width, slot,
                                           facecolor=colour, edgecolor="none",
                                           zorder=3))

        for i in range(len(stack) - 1):
            if stack[i] != stack[i + 1]:
                continue                     # the fill already draws that edge
            ax.plot([left, right], [(i + 1) * slot] * 2,
                    color="white" if stack[i] else colour, linewidth=1.5, zorder=4)

    labels = [f"{g}\n{(animal_table.condition == g).sum()} animals\n"
              f"{(fov_table.condition == g).sum()} FOVs" for g in groups]
    ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylim(-2, 102)
    ax.set_ylabel("Fields with significant alignment (%)", fontsize=12)
    ax.tick_params(labelsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["left"].set_bounds(0, 100)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "figure_blocks.pdf"), dpi=600)
    plt.close(fig)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    show = lambda v: f"{v:.4f}"

    # --- steps 1 to 3, one row per FOV -------------------------------------
    angles_by_fov = read_fov_columns(INPUT_FILE, SHEET)
    fov_table = pd.DataFrame([
        dict(fov=fov,
             animal=animal_of(fov),
             condition=condition_of(fov),
             n_cells=len(angles),
             R=order_parameter(angles),
             rayleigh_p=rayleigh_p_value(angles))
        for fov, angles in angles_by_fov.items()])
    fov_table["aligned"] = (fov_table.rayleigh_p < ALPHA).astype(int)

    print(fov_table[["fov", "animal", "condition", "n_cells", "R", "rayleigh_p",
                     "aligned"]].to_string(index=False, float_format=show))

    animal_table = fov_table.groupby(["condition", "animal"], as_index=False).agg(
        n_fovs=("fov", "size"), mean_R=("R", "mean"),
        proportion_aligned=("aligned", "mean"))

    groups = [g for g in CONDITIONS if g in set(fov_table.condition)]
    print("\nPer condition:")
    print(fov_table.groupby("condition").agg(
        n_animals=("animal", "nunique"), n_fovs=("fov", "size"),
        n_aligned=("aligned", "sum"), proportion=("aligned", "mean")
    ).round(3).to_string())

    # --- step 4, animals are the unit --------------------------------------
    # animals are numbered 0, 1, 2... and it is those numbers that get shuffled
    animal_names = animal_table.animal.tolist()
    aligned_by_animal = {i: fov_table.loc[fov_table.animal == name, "aligned"]
                                     .to_numpy(dtype=float)
                         for i, name in enumerate(animal_names)}
    animals_in_group = {g: [i for i, name in enumerate(animal_names)
                            if animal_table.condition.iloc[i] == g]
                        for g in groups}

    pairwise = pd.DataFrame([
        pairwise_test(animals_in_group, aligned_by_animal, a, b)
        for a, b in combinations(groups, 2)])
    pairwise["p_holm"] = holm_correction(pairwise.p_value.values)

    print("\nPairwise permutation tests:")
    print(pairwise.to_string(index=False, float_format=show))

    fov_table.to_csv(os.path.join(OUTDIR, "1_fov_results.csv"), index=False)
    animal_table.to_csv(os.path.join(OUTDIR, "2_animal_results.csv"), index=False)
    pairwise.to_csv(os.path.join(OUTDIR, "3_tests.csv"), index=False)
    make_figure_blocks(fov_table, animal_table, groups)
    print(f"\nSaved to {OUTDIR}/")


if __name__ == "__main__":
    main()