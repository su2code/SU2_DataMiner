# Generate a 3D (ProgressVariable, EnthalpyTot, MixtureFraction) LUT from premixed
# hydrogen-air flamelets.

# The table generator draws its figures with matplotlib; this script only saves them
# to PNG, so select a non-interactive backend before the import pulls matplotlib in.
import os
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

from Common.DataDrivenConfig import Config_FGM
from Manifold_Generation.LUT.LUTGenerators import SU2TableGenerator_FGM

# Loading configuration.
Config = Config_FGM("TableGeneration.cfg")

# Initializing table module. The table is spanned by the controlling variables of the
# configuration: progress variable and total enthalpy in the plane of each table level,
# mixture fraction across the levels.
Tgen = SU2TableGenerator_FGM(Config)

# Interpolator settings for evaluating the flamelet data on the table nodes.
Tgen.setNNearestNeighbors(14)
Tgen.setInverseDistanceExponent(3)

# Manually set mixture fraction limits.
Tgen.setTableLimits(0.00939225575395504, 0.0144703624619207)
Tgen.setNTableLevels(20)

# Generate table connectivity and interpolate flamelet data.
Tgen.generateTable()

# Visualize the table level closest to equivalence ratio 0.5.
cv_target = Config.GetUnburntScalars(equivalence_ratio=0.5, temperature=300.0)
pv_target = cv_target[0]
z_target = cv_target[2]
print("Target unburnt progress variable: ", pv_target)

level_index = int(np.argmin(np.absolute(Tgen.getTableLevels() - z_target)))

# Interpolated reaction rate on that level.
Tgen.visualizeTableLevel(level_index, var_to_plot="ProdRateTot_PV",
                         save_path="prodrate_level_%02i.png" % level_index)

# Table level connectivity on that level.
Tgen.visualizeTableLevel(level_index, save_path="mesh_level_%02i.png" % level_index)

# Write SU2 .drg table file.
Tgen.writeSU2Table()
