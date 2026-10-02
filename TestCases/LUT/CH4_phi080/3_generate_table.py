# Generate a 2D (ProgressVariable, EnthalpyTot) LUT for a single equivalence ratio phi = 0.80,
# from free flames, burner-stabilized flames and their cooled continuation to equilibrium.

# This script only saves figures to PNG, so select a non-interactive backend before the imports
# below pull matplotlib in.
import os
os.environ.setdefault("MPLBACKEND", "Agg")

from Common.DataDrivenConfig import Config_FGM
from Manifold_Generation.LUT.LUTGenerators import SU2TableGenerator_FGM

# Loading configuration. The mixture bounds are equal (see 0_generate_config.py), so the manifold
# is described by progress variable and total enthalpy and the table is a single 2D level.
Config = Config_FGM("TableGeneration.cfg")

Tgen = SU2TableGenerator_FGM(Config)

# Iterate to achieve a target number of nodes, or manually specify a coarse cell size.
# Tgen.setTargetNodeCount(4000)


# 100x100 = 0.01 = 1e-2
Tgen.setMaximumCellSize(2.0e-2)

# Bins along the enthalpy used to resolve the boundary at maximum progress variable. That boundary
# is the equilibrium locus, which shifts with enthalpy and runs almost parallel to the enthalpy
# axis, so the progress variable sweep confines it to its last few bins and cannot resolve it.
Tgen.setBoundaryTransverseBinCount(40)

# Apply refinement based on the values of thermochemical quantities:
# coef = 0.5 halves the edge lengths within the bounds
#Tgen.applyRefinementWithin("Temperature", lowerbound=270, upperbound=450.0, coef=0.5)

# Scale refinement based on the gradients of quantities:
Tgen.applyRefinementForGradientOf("ProdRateTot_PV", coef=0.3)

# Apply refinement in proximity of the reactants and products.
Tgen.refineEquilibrium(coef=0.5, margin=2e-2)

# Near the reactant and product boundaries, set the source terms of CO, H2, CO2 and H2O to zero
# where the absolute value of the source terms is below the tolerance.
Tgen.setSourceTermClamping(species_list=["CO", "H2", "CO2", "H2O"], pv_margin=1e-2, abs_tol=0.1)

# Optionally: apply smoothing to table data to get rid of any waves or discontinuities.
# Tgen.setSmoothingParameter(0.1)

Tgen.generateTable()

# Export table in vtk format (view in ParaView).
Tgen.writeParaviewTable("LUT_vtk")

# ── Diagnostics: the generator's own views of the single 2D table level ───────
# Mesh and its perimiter nodes.
Tgen.visualizeTableLevel(0, save_path="mesh.png")

# Mesh perimiter against the flamelet data it was extracted from: data outside the perimiter is
# dropped by the table, perimiter enclosing no data is extrapolated into. The counts of enclosed
# and discarded data points are reported to the terminal.
Tgen.visualizeTableLevelPerimiter(0, save_path="mesh_boundary_data.png")

# Table variables over the table plane.
for table_var, figure_name in (("ProdRateTot_PV", "prodrate.png"),
                               ("Temperature",    "temperature.png"),
                               ("Y_dot_net-CO",   "co_net.png"),
                               ("Y_dot_pos-CO",   "co_pos.png"),
                               ("Y_dot_neg-CO",   "co_neg.png"),
                               ("Y_dot_net-NOx",  "nox_net.png")):
    Tgen.visualizeTableLevel(0, var_to_plot=table_var, show_grid=False, save_path=figure_name)

# Slices at fixed enthalpy, sampled from the interpolated manifold rather than from the generated
# table, so these show what the table was built to represent. Query points outside the flamelet
# data support are left blank instead of extrapolated.
Tgen.plotTableSlices("ProgressVariable", "Temperature", "EnthalpyTot", n_slices=10,
                     save_path="T_vs_PV_slices.png")

# Write SU2 .drg table file (Dragon v1.0.1 format, 2D).
Tgen.writeSU2Table()
