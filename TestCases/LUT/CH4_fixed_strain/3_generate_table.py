# Generate a 2D (MixtureFraction, EnthalpyTot) LUT from counterflow diffusion
# flames at fixed strain rate.

# Limit inner thread pools BEFORE any library imports. gmsh's meshing kernel
# and the KD-tree/BLAS backends can spawn competing thread pools; under CPU
# contention this combination has been observed to segfault natively inside
# gmsh (not a catchable Python exception), so keep everything single-threaded.
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

# This script only saves figures to PNG; it never needs an interactive display.
# Force a non-interactive backend so it doesn't try (and fail) to initialize a
# Qt/X11 GUI in headless runs.
os.environ.setdefault("MPLBACKEND", "Agg")

from Common.DataDrivenConfig import Config_FGM
from Manifold_Generation.LUT.LUTGenerators import SU2TableGenerator_FGM

# Loading configuration.
Config = Config_FGM("TableGeneration.cfg")

Tgen = SU2TableGenerator_FGM(Config)

# The manifold is described by three controlling variables (see 0_generate_config.py),
# but this table is spanned by mixture fraction and total enthalpy alone: at fixed
# strain rate the counterflow state is a function of those two. Selecting two plane
# controlling variables and no level controlling variable yields a single 2D table;
# ProgressVariable is then carried along as an ordinary table variable.
Tgen.setTableAxes(["MixtureFraction", "EnthalpyTot"])

# Coarse cell size in the scaled table plane.
Tgen.setMaximumCellSize(5e-3)

# Bins along the mixture fraction used to extract the envelope bounding the table.
# The flamelets run diagonally across this plane, so within a bin the enthalpy of
# the boundary varies over the bin width; the envelope is placed at the extreme
# within each bin and is therefore biased outward by roughly half of that
# variation. More bins shrink the bias (100 bins: +3.9% enclosed area, 200: +1.9%,
# 400: +0.9%) at the cost of more boundary vertices.
Tgen.setBoundaryBinCount(200)

# Refine towards the reaction zone, which in this plane is a thin diagonal band.
Tgen.applyRefinementForGradientOf("Heat_Release", coef=0.3)

Tgen.generateTable()

# Export table in vtk format (view in ParaView).
Tgen.writeParaviewTable("LUT_vtk")

# ── Diagnostics: the generator's own views of the single 2D table level ───────
# Mesh and its perimiter nodes.
Tgen.visualizeTableLevel(0, save_path="mesh_ZH.png")

# Mesh perimiter against the flamelet data it was extracted from: data outside the
# perimiter is dropped by the table, perimiter enclosing no data is extrapolated into.
Tgen.visualizeTableLevelPerimiter(0, save_path="mesh_boundary_data.png")

# Table variables over the table plane.
Tgen.visualizeTableLevel(0, var_to_plot="Temperature", show_grid=False, save_path="T_ZH_2d.png")
Tgen.visualizeTableLevel(0, var_to_plot="Heat_Release", show_grid=False, save_path="HRR_ZH_2d.png")
Tgen.visualizeTableLevel(0, var_to_plot="Temperature", plot_3d=True, save_path="T_ZH_3d.png")

# Slices at fixed enthalpy, sampled from the interpolated manifold rather than from the
# generated table, so these show what the table was built to represent. Query points
# outside the flamelet data support are left blank instead of extrapolated.
Tgen.plotTableSlices("MixtureFraction", "Temperature", "EnthalpyTot", n_slices=10,
                     save_path="T_vs_Z_slices.png")
Tgen.plotTableSlices("MixtureFraction", "Heat_Release", "EnthalpyTot", n_slices=10,
                     save_path="HRR_vs_Z_slices.png")

# Write SU2 .drg table file (Dragon v1.0.1 format, 2D, MixtureFraction + EnthalpyTot).
Tgen.writeSU2Table()
