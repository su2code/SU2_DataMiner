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

#
Tgen.SetEquivalenceRatioLimits(phi_min=0.30, phi_max=0.90)
Tgen.SetNTableLevels(3)
Tgen.SetRefinementFields(["ProdRateTot_PV","Heat_Release"])
Tgen.SetRefinementMethod("gradient")

# coarse
#Tgen.SetBaseCellSize(1e-2)
#Tgen.SetRefinedCellSize(0.5e-2)
#Tgen.SetRefinementRadius(2.0e-2)
#Tgen.SetMaxRefinementSeeds(500)
#Tgen.SetHullCellSize(1.0e-2)

# Generate table connectivity and interpolate flamelet data.
Tgen.GenerateTableNodes()

Tgen.SetTableAxes(level_cv_name="MixtureFraction",
                       plane_cv_names=["ProgressVariable", "EnthalpyTot"])

# ── Diagnostics: inspect table quality before writing ──────────────────────
# T(Z) for 10 enthalpy levels spanning the full data range.
_show(Tgen.PlotTableSlices,
    x_cv        = "MixtureFraction",
    y_var       = "Temperature",
    slice_cv    = "EnthalpyTot",
    slice_range = None,
    n_slices    = 10,
    n_x_points  = 300,
    save_path   = _save("T_vs_Z_slices.png"))

# Heat_Release(Z) for the same 10 enthalpy levels.
_show(Tgen.PlotTableSlices,
    x_cv        = "MixtureFraction",
    y_var       = "Heat_Release",
    slice_cv    = "EnthalpyTot",
    slice_range = None,
    n_slices    = 10,
    n_x_points  = 300,
    save_path   = _save("HRR_vs_Z_slices.png"))

# Visualize the table mesh and reaction rate at phi = 0.55.
cv_target = Config.GetUnburntScalars(equivalence_ratio=0.60, temperature=300.0)
pv_target = cv_target[0]
z_target  = cv_target[2]
print("Target unburnt progress variable:", pv_target)

# Mesh only:
_show(Tgen.VisualizeTableLevel, z_target,
      save_path=_save("mesh_cH.png"))

# 2D colour map of Temperature (no wireframe):
_show(Tgen.VisualizeTableLevel, z_target, "Temperature",
      plot_3d=False, show_grid=False, save_path=_save("T_cH_2d.png"))

# 3D surface of Temperature (with wireframe):
_show(Tgen.VisualizeTableLevel, z_target, "Temperature",
      plot_3d=True, show_grid=True, save_path=_save("T_phi060_3d.png"))
# Production rate of progress variable:
_show(Tgen.VisualizeTableLevel, z_target, "ProdRateTot_PV",
      plot_3d=False, show_grid=False, save_path=_save("ProdRate_phi060_2d.png"))
_show(Tgen.VisualizeTableLevel, z_target, "ProdRateTot_PV",
      plot_3d=True, show_grid=False, save_path=_save("ProdRate_phi060_3d.png"))

# from 99% of the max progress variable, set the source terms of H2 to zero
# if the absolute value of the source terms is |S| < 0.1
Tgen.ClampSourceTerms(species_list=["H2", "H2O"], pv_frac=0.99, abs_tol=0.1)

# Write SU2 .drg table file (Dragon v1.0.1 format, 2D).
Tgen.WriteTableFile()

# Keep all figure windows open until the user closes them.
if SHOW_FIGURES:
    plt.show()
