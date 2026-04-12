from Manifold_Generation.LUT.FlameletTableGeneration import SU2TableGenerator
from Common.DataDrivenConfig import Config_FGM 

# Loading configuration.
Config = Config_FGM("TableGeneration.cfg")

# Initializing table module and pre-process interpolator.
Tgen = SU2TableGenerator(Config, n_near=14, p_fac=3)

# Manually set mixture fraction limits.
# here: phi = [0.5, 1.0] so use between [0.55, 0.95]

# 0.8 - 0.81
Tgen.SetMixtureFractionLimits(mix_frac_max=0.0451748, mix_frac_min=0.044642)
#Tgen.SetMixtureFractionLimits(mix_frac_max=0.0525, mix_frac_min=0.0312)
Tgen.SetNTableLevels(2)

# Visualize the interpolated reaction rate at equivalence ratio 0.8.
cv_target = Config.GetUnburntScalars(equivalence_ratio=0.805, temperature=270.0)
pv_target = cv_target[0]
z_target = cv_target[2]
print("Target unburnt progress variable: ", pv_target)
Tgen.VisualizeTableLevel(z_target, "ProdRateTot_PV")

# Visualize the interpolated reaction rate at equivalence ratio 0.8.
cv_target = Config.GetUnburntScalars(equivalence_ratio=0.805, temperature=270.0)
pv_target = cv_target[0]
z_target = cv_target[2]
print("Target unburnt progress variable: ", pv_target)
Tgen.VisualizeTableLevel(z_target, "Heat_Release")

cv_target = Config.GetUnburntScalars(equivalence_ratio=0.805, temperature=270.0)
pv_target = cv_target[0]
z_target = cv_target[2]
print("Target unburnt progress variable: ", pv_target)
Tgen.VisualizeTableLevel(z_target, "Temperature")

#cv_target = Config.GetUnburntScalars(equivalence_ratio=0.70, temperature=270.0)
#pv_target = cv_target[0]
#z_target = cv_target[2]
#print("Target unburnt progress variable: ", pv_target)
#Tgen.VisualizeTableLevel(z_target, "Temperature")

#cv_target = Config.GetUnburntScalars(equivalence_ratio=0.60, temperature=270.0)
#pv_target = cv_target[0]
#z_target = cv_target[2]
#print("Target unburnt progress variable: ", pv_target)
#Tgen.VisualizeTableLevel(z_target, "Temperature")

# Visualize table level connectivity at equivalence ratio 0.5.
Tgen.VisualizeTableLevel(z_target)

# Generate table connectivity and interpolate flamelet data.
Tgen.GenerateTableNodes()

# Write SU2 .drg table file.
Tgen.WriteTableFile()
