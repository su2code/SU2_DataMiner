from Manifold_Generation.LUT.FlameletTableGeneration import SU2TableGenerator
from Common.DataDrivenConfig import FlameletAIConfig

# Loading configuration.
Config = Config_FGM("TableGeneration.cfg")

# Initializing table module and pre-process interpolator.
Tgen = SU2TableGenerator_FGM(Config)

# Distribute tabulation process over 4 cores.
Tgen.SetNCores(4)

# Manually set mixture fraction limits.
Tgen.SetMixtureFractionLimits(mix_frac_max=0.02, mix_frac_min=0.009392)

# Use 50 table levels and approximately 2k nodes per table level.
Tgen.SetNTableLevels(10)
Tgen.SetNnodes_Target(5000)
Tgen.ConditionalRefinement("Temperature", lowerbound=600,upperbound=700, coef=0.5)
Tgen.RefineReactants(0.1)
# Tgen.ConditionalRefinement("Temperature", lowerbound=1000, coef=0.2)

# Reaction rate and preferential diffusion closure sources are zeroed within thin dead
# bands in normalized reaction progress at the unburnt and burnt ends of each table level
# (defaults shown; applied automatically during GenerateTableNodes).
Tgen.SetSourceDeadBands(c_unburnt=1e-3, c_burnt=1e-2)


# Insert mixture fraction level at equivalence ratio of 0.5.
cv_target = Config.GetUnburntScalars(equivalence_ratio=0.5, temperature=300.0)
z_target = cv_target[2]
Tgen.InsertMixtureFractionLevel(z_target)

# Visualize table level connectivity at equivalence ratio 0.5.
Tgen.VisualizeTableLevel(z_target)

# Generate table connectivity and interpolate flamelet data.
Tgen.generateTableNodes()

# Write SU2 .drg table file.
Tgen.writeSU2Table()