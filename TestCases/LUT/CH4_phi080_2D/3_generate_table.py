from Manifold_Generation.LUT.FlameletTableGeneration import SU2TableGenerator
from Common.DataDrivenConfig import Config_FGM
import matplotlib.pyplot as plt

# Loading configuration.
Config = Config_FGM("TableGeneration.cfg")

# Initializing table module and pre-process interpolator.
Tgen = SU2TableGenerator(Config, n_near=14, p_fac=1)

# Generate a 2D (ProgressVariable, EnthalpyTot) LUT for a single equivalence
# ratio phi = 0.80.  Because phi_min == phi_max, the table generator writes a
# Dragon v1.0.1 file without mixture-fraction levels.
Tgen.SetEquivalenceRatioLimits(phi_min=0.80, phi_max=0.80)
Tgen.SetNTableLevels(1)
Tgen.SetRefinementFields(["ProdRateTot_PV", "Y_dot_net-CO", "Y_dot_pos-CO","Y_dot_neg-CO"])

# small
#Tgen.SetBaseCellSize(1e-2)
#Tgen.SetRefinedCellSize(1e-2)
#Tgen.SetRefinementRadius(1e-2)
#Tgen.SetRefinementMethod("gradient")
#Tgen.SetMaxRefinementSeeds(500)
#Tgen.SetHullCellSize(1.0e-2)

# medium 
Tgen.SetBaseCellSize(5e-3)
Tgen.SetRefinedCellSize(5e-3)
Tgen.SetRefinementRadius(5e-3)
Tgen.SetRefinementMethod("gradient")
Tgen.SetMaxRefinementSeeds(500)
Tgen.SetHullCellSize(5.0e-3)

# fine
#Tgen.SetBaseCellSize(2.5e-3)
#Tgen.SetRefinedCellSize(1.0e-3)
#Tgen.SetRefinementRadius(2.5e-3)
#Tgen.SetRefinementMethod("gradient")
#Tgen.SetMaxRefinementSeeds(500)
#Tgen.SetHullCellSize(2.5e-3)




Tgen.SetTableAxes(level_cv_name="MixtureFraction",
                       plane_cv_names=["ProgressVariable", "EnthalpyTot"])

# Visualize the table mesh and reaction rate at phi = 0.80.
cv_target = Config.GetUnburntScalars(equivalence_ratio=0.80, temperature=300.0)
pv_target = cv_target[0]
z_target  = cv_target[2]
print("Target unburnt progress variable:", pv_target)
Tgen.VisualizeTableLevel(z_target, "ProdRateTot_PV", save_path="prodrate.png", show=True, show_grid=False)
Tgen.VisualizeTableLevel(z_target, "Temperature", save_path="temperature.png",show=True, show_grid=False)
Tgen.VisualizeTableLevel(z_target, "Y_dot_net-CO", save_path="co_net.png",show=True, show_grid=False)
Tgen.VisualizeTableLevel(z_target, "Y_dot_pos-CO", save_path="co_pos.png",show=True, show_grid=False)
Tgen.VisualizeTableLevel(z_target, "Y_dot_neg-CO", save_path="co_neg.png",show=True, show_grid=False)
Tgen.VisualizeTableLevel(z_target, "Y_dot_net-NOx", save_path="nox_net.png",show=True, show_grid=False)
Tgen.VisualizeTableLevel(z_target, show=True)

# Show all plots
plt.show()

# Generate table connectivity and interpolate flamelet data.
Tgen.GenerateTableNodes()
# from 99% of the max progress variable, set the source terms of CO and H2 to zero
# if the absolute value of the source terms is |S| < 0.1
Tgen.ClampSourceTerms(species_list=["CO", "H2", "CO2", "H2O"], pv_frac=0.99, abs_tol=0.1)

# Write SU2 .drg table file (Dragon v1.0.1 format, 2D).
Tgen.WriteTableFile()
