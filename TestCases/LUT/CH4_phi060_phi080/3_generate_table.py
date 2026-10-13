from Manifold_Generation.LUT.FlameletTableGeneration import SU2TableGenerator
from Common.DataDrivenConfig import Config_FGM
import matplotlib.pyplot as plt


# ── Output options ──────────────────────────────────────────────────────────
SHOW_FIGURES = True    # keep all figures open (non-blocking) while script runs
SAVE_FIGURES = True    # save all figures as PNG files in the working directory

def _save(name):
    return name if SAVE_FIGURES else None

def _show(fig_fn, *args, **kwargs):
    """Call a visualisation method once with save_path and show set from the global flags."""
    if not SHOW_FIGURES and not SAVE_FIGURES:
        return
    kwargs['save_path'] = kwargs.get('save_path') if SAVE_FIGURES else None
    kwargs['show'] = SHOW_FIGURES
    fig_fn(*args, **kwargs)



# Loading configuration.
Config = Config_FGM("TableGeneration.cfg")

# Initializing table module and pre-process interpolator.
Tgen = SU2TableGenerator(Config, n_near=14, p_fac=1)

# Performance optimizations
Tgen.SetNCores(4)  # Enable parallel processing with 4 cores
Tgen.SetCurvatureGridResolution(300)  # Reduce grid resolution for faster curvature computation (default: 800)

# Generate a 2D (ProgressVariable, EnthalpyTot) LUT for a single equivalence
# ratio phi = 0.80.  Because phi_min == phi_max, the table generator writes a
# Dragon v1.0.1 file without mixture-fraction levels.
Tgen.SetEquivalenceRatioLimits(phi_min=0.7, phi_max=0.7)
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

# Generate table connectivity and interpolate flamelet data.
Tgen.GenerateTableNodes()

# Visualize the table mesh and reaction rate at phi = 0.80.
cv_target = Config.GetUnburntScalars(equivalence_ratio=0.70, temperature=300.0)
pv_target = cv_target[0]
z_target  = cv_target[2]
print("Target unburnt progress variable:", pv_target)

# Mesh only:
_show(Tgen.VisualizeTableLevel, z_target,
      save_path=_save("mesh_phi070.png"))
# 2D colour map of Temperature (no wireframe):
_show(Tgen.VisualizeTableLevel, z_target, "Temperature",
      plot_3d=False, show_grid=False, save_path=_save("T_phi070_2d.png"))
# 3D surface of Temperature (with wireframe):
_show(Tgen.VisualizeTableLevel, z_target, "Temperature",
      plot_3d=True, show_grid=True, save_path=_save("T_phi070_3d.png"))
# Production rate of progress variable:
_show(Tgen.VisualizeTableLevel, z_target, "ProdRateTot_PV",
      plot_3d=False, show_grid=False, save_path=_save("ProdRate_phi070_2d.png"))
_show(Tgen.VisualizeTableLevel, z_target, "ProdRateTot_PV",
      plot_3d=True, show_grid=False, save_path=_save("ProdRate_phi070_3d.png"))
# Optional: CO source term visualizations (uncomment if needed)
#_show(Tgen.VisualizeTableLevel, z_target, "Y_dot_net-CO",
#      plot_3d=False, show_grid=False, save_path=_save("Y_dot_net_CO_2d.png"))
#_show(Tgen.VisualizeTableLevel, z_target, "Y_dot_pos-CO",
#      plot_3d=False, show_grid=False, save_path=_save("Y_dot_pos_CO_2d.png"))
#_show(Tgen.VisualizeTableLevel, z_target, "Y_dot_neg-CO",
#      plot_3d=False, show_grid=False, save_path=_save("Y_dot_neg_CO_2d.png"))

# from 99% of the max progress variable, set the source terms of CO and H2 to zero
# if the absolute value of the source terms is |S| < 0.1
Tgen.ClampSourceTerms(species_list=["CO", "H2", "CO2", "H2O"], pv_frac=0.99, abs_tol=0.1)

# Write SU2 .drg table file (Dragon v1.0.1 format, 2D).
Tgen.WriteTableFile()

# Keep all figure windows open until the user closes them.
if SHOW_FIGURES:
    plt.show()
