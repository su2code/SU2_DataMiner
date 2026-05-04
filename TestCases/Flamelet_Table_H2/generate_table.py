from su2dataminer.manifold import SU2TableGenerator_FGM
from su2dataminer.config import Config_FGM

# Loading configuration.
Config = Config_FGM("TableGeneration.cfg")

# Initializing table module and pre-process interpolator.
Tgen = SU2TableGenerator_FGM(Config, n_near=10, p_fac=2)
Tgen.SetNCores(4)
# Manually set mixture fraction limits.
Tgen.SetMixtureFractionLimits(mix_frac_max=0.02, mix_frac_min=0.00939225575395504)
Tgen.SetNTableLevels(10)

Tgen.SetNnodes_Target(6000)
# Visualize the interpolated reaction rate at equivalence ratio 0.5.
cv_target = Config.GetUnburntScalars(equivalence_ratio=0.5, temperature=300.0)
pv_target = cv_target[0]
z_target = cv_target[2]
# # Visualize table level connectivity at equivalence ratio 0.5.
Tgen.InsertMixtureFractionLevel(z_target)
# Visualize table level connectivity at equivalence ratio 0.5.
Tgen.VisualizeTableLevel(z_target)

# # Generate table connectivity and interpolate flamelet data.
Tgen.GenerateTableNodes()

# Write SU2 .drg table file.
Tgen.WriteOutParaview("table_paraview")
Tgen.WriteTableFile()