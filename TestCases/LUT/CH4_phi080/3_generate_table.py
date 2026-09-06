import os
os.environ.setdefault("MPLBACKEND", "Agg")

import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection

from Manifold_Generation.LUT.LUTGenerators import SU2TableGenerator_FGM
from Common.DataDrivenConfig import Config_FGM
from Data_Generation.DataGenerator_FGM import FlameletSolverDict

# Loading configuration.
Config = Config_FGM("TableGeneration.cfg")

Tgen = SU2TableGenerator_FGM(Config)

# Iterate to achieve a target number of nodes, or manually specify a coarse cell size.
# Tgen.setTargetNodeCount(4000)
Tgen.setMaximumCellSize(1.5e-2)

# Resolution of the uniform reference lattice spanning the table level. The lattice supplies the
# reference point cloud from which the gradient normalization factors below are derived.
Tgen.setPointCloudResolution(400)

# Bins along the progress variable used to extract the envelope. More bins follow the data more
# closely, but also start resolving the discrete spacing between individual flamelets.
Tgen.setBoundaryBinCount(100)

# Bins along the enthalpy used to resolve the boundary at maximum progress variable. That boundary
# is the equilibrium locus, which shifts with enthalpy and runs almost parallel to the enthalpy
# axis, so the sweep above confines it to its last few bins and cannot resolve it.
Tgen.setBoundaryTransverseBinCount(40)

# Apply refinement based on the values of thermochemical quantities:
Tgen.applyRefinementWithin("Temperature", lowerbound=270, upperbound=450.0, coef=0.5)

# Scale refinement based on the gradients of quantities:
Tgen.applyRefinementForGradientOf("ProdRateTot_PV", coef=0.3)

# Apply refinement in proximity of the reactants and products.
Tgen.refineEquilibrium(coef=0.5, margin=2e-2)

# Optionally: apply smoothing to table data to get rid of any waves or discontinuities.
# Tgen.setSmoothingParameter(0.1)

Tgen.generateTable()

# Export table in vtk format (view in ParaView).
Tgen.writeParaviewTable("LUT_vtk")

# ── Diagnostics: quick matplotlib views of the single 2D table level ────────
# SU2TableGenerator_FGM has no built-in plotting (writeParaviewTable/ParaView
# is the intended way to inspect it); build simple ones here from its exposed
# table data so we still get a quick sanity check without leaving Python.
nodes = Tgen._table_nodes[0]          # (Np, 2): [ProgressVariable, EnthalpyTot]
connectivity = Tgen._table_connectivity[0]
table_data = Tgen._data_in_table[0]
hull_idx = Tgen._table_hullnodes[0]
triang = mtri.Triangulation(nodes[:, 0], nodes[:, 1], connectivity)


def _plot_field(varname, save_path):
    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
    tc = ax.tripcolor(triang, table_data[varname].to_numpy(), shading='gouraud', cmap='inferno')
    cb = fig.colorbar(tc, ax=ax, pad=0.02)
    cb.set_label(varname, fontsize=12)
    ax.set_xlabel("ProgressVariable", fontsize=14)
    ax.set_ylabel("EnthalpyTot", fontsize=14)
    ax.set_title(varname, fontsize=14)
    fig.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved: %s" % save_path)


# Flamelet families behind the table. The concatenated manifold carries no record of which
# flamelet a data point came from, so the families are read back from their own folders. These
# are the Cantera solution nodes: the manifold resamples each flamelet onto a fixed number of
# points along its own arc length (see 2_collect_flamelet_data.py), so the points below lie on
# the same curves but not at the same spacing.
FLAMELET_FAMILIES = [("FREEFLAME",       "red",    "free flame"),
                     ("BURNERFLAME",     "blue",   "burner-stabilized"),
                     ("INT_BURNERFLAME", "purple", "interpolated burner flame"),
                     ("EQUILIBRIUM",     "green",  "equilibrium")]

# Rows taken from each equilibrium products file. Keep this equal to the value passed to
# SetNEquilibriumNodes in 2_collect_flamelet_data.py: the products enter the manifold as that
# many states rather than as a curve, and plotting the whole curve would show data the table
# never saw.
N_EQUILIBRIUM_NODES = 1


def _flamelet_data_by_family():
    """Controlling variables of the flamelet solutions behind the table, per flamelet family."""
    cv_names = Config.GetControllingVariables()
    data_per_family = {}
    for flamelet_type, colour, label in FLAMELET_FAMILIES:
        solver = FlameletSolverDict[flamelet_type](Config)
        flamelet_files = sorted(glob.glob(os.sep.join((Config.GetOutputDir(),
                                                       solver.getFlameletFolder(), "*", "*.csv"))))
        cv_data = []
        for flamelet_file in flamelet_files:
            solution = pd.read_csv(flamelet_file)
            if "Products" in os.path.basename(flamelet_file):
                solution = solution.iloc[:N_EQUILIBRIUM_NODES, :]
            # The progress variable is not stored with the solution; it is evaluated from the
            # mass fractions exactly as the data collection step does.
            progress_variable = Config.ComputeProgressVariable(list(solution.keys()), solution.values)
            cv_data.append(np.column_stack((progress_variable, solution[cv_names[1]].to_numpy())))
        if cv_data:
            data_per_family[flamelet_type] = (np.vstack(cv_data), colour, label)
    return data_per_family


def _plot_mesh_with_data(save_path):
    """Mesh and its boundary against the flamelet families the table is built from."""
    data_per_family = _flamelet_data_by_family()

    # An edge shared by a single cell lies on the boundary of the mesh. The perimeter node
    # indices are not ordered along the perimeter, so the boundary is collected from the
    # connectivity instead of by walking those nodes.
    edges = np.sort(np.vstack((connectivity[:, [0, 1]],
                               connectivity[:, [1, 2]],
                               connectivity[:, [2, 0]])), axis=1)
    unique_edges, edges_per_cell = np.unique(edges, axis=0, return_counts=True)
    boundary_segments = nodes[unique_edges[edges_per_cell == 1]][:, :, :2]

    # At table resolution the data points and the cells overlap into solid colour, so the full
    # view is paired with a detail of the burner-stabilized limit at mid-domain, where
    # individual cells, boundary edges and data points are separable.
    x_span = nodes[:, 0].max() - nodes[:, 0].min()
    y_span = nodes[:, 1].max() - nodes[:, 1].min()
    x_detail = 0.5*(nodes[:, 0].min() + nodes[:, 0].max())
    near_middle = np.absolute(boundary_segments[:, :, 0].mean(axis=1) - x_detail) < 0.05*x_span
    y_detail = boundary_segments[near_middle][:, :, 1].min()

    fig, axes = plt.subplots(1, 2, figsize=(17, 7), constrained_layout=True)
    for ax, detail in zip(axes, (False, True)):
        ax.triplot(nodes[:, 0], nodes[:, 1], connectivity,
                   lw=0.4 if detail else 0.2, color='0.6', zorder=2)
        for family_data, colour, label in data_per_family.values():
            # A family contributing only a handful of states would vanish at the marker size the
            # dense families need.
            marker_size = 4.0 if len(family_data) < 500 else 1.6
            ax.plot(family_data[:, 0], family_data[:, 1], '.', color=colour, zorder=1,
                    ms=2.5*marker_size if detail else marker_size,
                    label="%s (%i points)" % (label, len(family_data)))
        ax.add_collection(LineCollection(boundary_segments, colors='k',
                                         linewidths=2.0 if detail else 1.4, zorder=3))
        ax.plot([], [], '-', color='k', lw=1.6, label="mesh boundary")
        ax.plot([], [], '-', color='0.6', lw=1.0, label="mesh (%i nodes)" % len(nodes))
        if detail:
            ax.set_xlim(x_detail - 0.06*x_span, x_detail + 0.06*x_span)
            ax.set_ylim(y_detail - 0.02*y_span, y_detail + 0.10*y_span)
            ax.set_title("detail of the burner-stabilized limit", fontsize=12)
        else:
            ax.set_title("full table level", fontsize=12)
            ax.legend(markerscale=6, fontsize=9, loc='lower left')
        ax.set_xlabel("ProgressVariable", fontsize=13)
        ax.set_ylabel("EnthalpyTot", fontsize=13)

    fig.suptitle("Mesh boundary and the flamelet data the table is built from", fontsize=14)
    fig.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved: %s" % save_path)


_plot_mesh_with_data("mesh_boundary_data.png")

fig, ax = plt.subplots(figsize=(10, 10))
ax.triplot(nodes[:, 0], nodes[:, 1], connectivity, linewidth=0.5)
ax.plot(nodes[hull_idx, 0], nodes[hull_idx, 1], 'ko', ms=3, label="Hull nodes")
ax.set_xlabel("ProgressVariable", fontsize=14)
ax.set_ylabel("EnthalpyTot", fontsize=14)
ax.set_title("Mesh (phi = 0.80)", fontsize=14)
ax.legend(fontsize=12)
fig.savefig("mesh.png", dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved: mesh.png")

_plot_field("ProdRateTot_PV", "prodrate.png")
_plot_field("Temperature", "temperature.png")
_plot_field("Y_dot_net-CO", "co_net.png")
_plot_field("Y_dot_pos-CO", "co_pos.png")
_plot_field("Y_dot_neg-CO", "co_neg.png")
_plot_field("Y_dot_net-NOx", "nox_net.png")

# Write SU2 .drg table file (Dragon v1.0.1 format, 2D).
Tgen.ClampSourceTerms(species_list=["CO", "H2", "CO2", "H2O"], pv_frac=0.99, abs_tol=0.1)
Tgen.writeSU2Table()


