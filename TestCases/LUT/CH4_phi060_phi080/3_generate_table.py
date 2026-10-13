# Generate a 3D (ProgressVariable, EnthalpyTot, MixtureFraction) LUT for equivalence ratios
# phi = 0.60 to 0.80, with one 2D table level per mixture fraction.

# This script only saves figures to PNG, so select a non-interactive backend before the imports
# below pull matplotlib in.
import os
os.environ.setdefault("MPLBACKEND", "Agg")

from Common.DataDrivenConfig import Config_FGM
from Manifold_Generation.LUT.LUTGenerators import SU2TableGenerator_FGM

Config = Config_FGM("TableGeneration.cfg")

Tgen = SU2TableGenerator_FGM(Config)
Tgen.setTableAxes(plane_cv_names=["ProgressVariable", "EnthalpyTot"], level_cv_name="MixtureFraction")
Tgen.setNTableLevels(5)

Tgen.setMaximumCellSize(2.0e-2)
Tgen.setBoundaryTransverseBinCount(40)

Tgen.applyRefinementForGradientOf("ProdRateTot_PV", coef=0.3)
Tgen.refineEquilibrium(coef=0.5, margin=2e-2)

# Near the reactant and product boundaries, set the source terms of CO, H2, CO2 and H2O to zero
# where the absolute value of the source terms is below the tolerance.
Tgen.setSourceTermClamping(species_list=["CO", "H2", "CO2", "H2O"], pv_margin=1e-2, abs_tol=0.1)

Tgen.generateTable()

Tgen.writeParaviewTable("LUT_vtk")

# ── Diagnostics: the generator's own views of the first and last table level ──
for level_index in (0, Tgen.getTableLevels().size - 1):
    Tgen.visualizeTableLevel(level_index, save_path="mesh_level%i.png" % level_index)
    Tgen.visualizeTableLevelPerimiter(level_index, save_path="mesh_boundary_data_level%i.png" % level_index)
    for table_var, figure_name in (("ProdRateTot_PV", "prodrate"),
                                   ("Temperature",    "temperature")):
        Tgen.visualizeTableLevel(level_index, var_to_plot=table_var, show_grid=False, \
                                 save_path="%s_level%i.png" % (figure_name, level_index))

# Write SU2 .drg table file (Dragon v1.0.1 format, 3D).
Tgen.writeSU2Table()
