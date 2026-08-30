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

# This script only saves figures to PNG (PlotTableSlices with show=False); it
# never needs an interactive display. Force a non-interactive backend so it
# doesn't try (and fail/abort) to initialize a Qt/X11 GUI in headless runs.
os.environ.setdefault("MPLBACKEND", "Agg")

from Manifold_Generation.LUT.FlameletTableGeneration import SU2TableGenerator
from Common.DataDrivenConfig import Config_FGM
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
# flames at fixed strain rate. ProgressVariable is the nominal "level" variable
# but is unused here (single-level 2D table).
Tgen.SetTableAxes(level_cv_name="ProgressVariable",
                  plane_cv_names=["MixtureFraction", "EnthalpyTot"])

# min == max triggers 2D mode: one mesh spanning the full (Z, h) data cloud.
Tgen.SetMixtureFractionLimits(mix_frac_min=0.0, mix_frac_max=0.0)
Tgen.SetNTableLevels(1)

# Refinement: use Temperature and Heat_Release as indicators.
# Heat_Release concentrates points near the reaction zone in (Z, h) space.
Tgen.SetRefinementFields(["Heat_Release", "Temperature"])
Tgen.SetRefinementMethod("gradient")

# Medium resolution.
Tgen.SetBaseCellSize(5e-3)
Tgen.SetRefinedCellSize(2e-3)
Tgen.SetRefinementRadius(5e-3)
Tgen.SetMaxRefinementSeeds(500)
Tgen.SetHullCellSize(5e-3)

# Generate table connectivity and interpolate flamelet data onto mesh nodes.
Tgen.generateTableNodes()

# ── Diagnostics: inspect table quality before writing ──────────────────────
# T(Z) for 10 enthalpy levels spanning the full data range.
Tgen.PlotTableSlices(
    x_cv="MixtureFraction", y_var="Temperature", slice_cv="EnthalpyTot",
    n_slices=10, n_x_points=300, save_path="T_vs_Z_slices.png")

# Heat_Release(Z) for the same 10 enthalpy levels.
Tgen.PlotTableSlices(
    x_cv="MixtureFraction", y_var="Heat_Release", slice_cv="EnthalpyTot",
    n_slices=10, n_x_points=300, save_path="HRR_vs_Z_slices.png")



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
Tgen.writeSU2Table()
