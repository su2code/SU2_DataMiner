from Manifold_Generation.LUT.FlameletTableGeneration import SU2TableGenerator
from Common.DataDrivenConfig import Config_FGM
import matplotlib
import matplotlib.pyplot as plt

# ── Output options ──────────────────────────────────────────────────────────
SHOW_FIGURES = True    # keep all figures open (non-blocking) while script runs
SAVE_FIGURES = True    # save all figures as PNG files in the working directory

if SHOW_FIGURES:
    plt.ion()   # interactive mode: plt.pause() displays without blocking

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
Tgen = SU2TableGenerator(Config, n_near=14, p_fac=3)

# Generate a 2D (MixtureFraction, EnthalpyTot) LUT from counterflow diffusion
# flames at fixed strain rate.  The plane axes are (Z, h); ProgressVariable is
# the nominal "level" variable but is unused (single-level 2D table).
Tgen.SetTableAxes(level_cv_name="ProgressVariable",
                  plane_cv_names=["MixtureFraction", "EnthalpyTot"])

# min == max triggers 2D mode: one mesh spanning the full (Z, h) data cloud.
Tgen.SetMixtureFractionLimits(mix_frac_min=0.0, mix_frac_max=0.0)
Tgen.SetNTableLevels(1)

# Refinement: use Temperature and Heat_Release as indicators.
# Heat_Release concentrates points near the reaction zone in (Z, h) space.
Tgen.SetRefinementFields(["Heat_Release", "Temperature"])
Tgen.SetRefinementMethod("gradient")

# medium resolution
Tgen.SetBaseCellSize(5e-3)
Tgen.SetRefinedCellSize(2e-3)
Tgen.SetRefinementRadius(5e-3)
Tgen.SetMaxRefinementSeeds(500)
Tgen.SetHullCellSize(5e-3)

# Generate table connectivity and interpolate flamelet data onto mesh nodes.
Tgen.GenerateTableNodes()

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

# Mesh only:
_show(Tgen.VisualizeTableLevel, 0.0,
      save_path=_save("mesh_ZH.png"))
# 2D colour map of Temperature (no wireframe):
_show(Tgen.VisualizeTableLevel, 0.0, "Temperature",
      plot_3d=False, show_grid=False, save_path=_save("T_ZH_2d.png"))
# 3D surface of Temperature (with wireframe):
_show(Tgen.VisualizeTableLevel, 0.0, "Temperature",
      plot_3d=True, show_grid=True, save_path=_save("T_ZH_3d.png"))

# Write SU2 .drg table file (Dragon v1.0.1 format, 2D, MixtureFraction + EnthalpyTot).
Tgen.WriteTableFile()

# Keep all figure windows open until the user closes them.
if SHOW_FIGURES:
    plt.ioff()
    plt.show(block=True)
