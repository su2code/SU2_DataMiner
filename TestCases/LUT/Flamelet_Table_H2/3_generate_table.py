import numpy as np
import matplotlib.pyplot as plt

from Common.DataDrivenConfig import Config_FGM
from Manifold_Generation.LUT.LUTGenerators import SU2TableGenerator_FGM


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
# n_near = number of nearest neighbors
Tgen = SU2TableGenerator_FGM(Config)
Tgen.setNNearestNeighbors(10)
Tgen.setInverseDistanceExponent(3)

# Performance optimizations
#Tgen.setNProcessors(4)  # Enable parallel processing with 4 cores

# The table plane is spanned by progress variable and total enthalpy, the levels by mixture fraction.
Tgen.setTableAxes(plane_cv_names=["ProgressVariable", "EnthalpyTot"],
                  level_cv_name="MixtureFraction")

# Table levels are set in mixture fraction, so the equivalence ratio limits are converted first.
# Deliberately narrower than the generated range: the manifold thins out towards its extremes.
def _mixture_fraction_at(equivalence_ratio):
    return Config.GetUnburntScalars(equivalence_ratio=equivalence_ratio, temperature=300.0)[2]

Tgen.setTableLimits(_mixture_fraction_at(0.25), _mixture_fraction_at(1.25))
#Tgen.setNTableLevels(50)
Tgen.setNTableLevels(10)

# Refinement is scaled by the gradients of these quantities rather than seeded around them, so the
# refined cell size enters as a factor on the coarse cell size (2.5e-3 / 5e-3 = 0.5) and the
# refinement radius and seed count have no counterpart.
Tgen.applyRefinementForGradientOf("ProdRateTot_PV", coef=0.5)
Tgen.applyRefinementForGradientOf("Heat_Release", coef=0.5)

# coarse
# note: mesh size is [1x1] so base cell size of 1e-2 means 100 points
#Tgen.setMaximumCellSize(1e-2)

Tgen.setMaximumCellSize(5e-3)

# Generate table connectivity and interpolate flamelet data. The diagnostics below draw on the
# generated table, so unlike in the previous table module they come after this rather than before.
Tgen.generateTable()

def _level_at(equivalence_ratio):
    """Index of the table level closest to a given equivalence ratio."""
    return int(np.argmin(np.absolute(Tgen.getTableLevels() - _mixture_fraction_at(equivalence_ratio))))


# ── Diagnostics: inspect table quality before writing ──────────────────────
# T(Z) for 10 enthalpy levels spanning the full data range.
_show(Tgen.plotTableSlices,
    x_cv        = "MixtureFraction",
    y_var       = "Temperature",
    slice_cv    = "EnthalpyTot",
    slice_range = None,
    n_slices    = 10,
    n_x_points  = 600,
    save_path   = _save("T_vs_Z_slices.png"))

# Heat_Release(Z) for the same 10 enthalpy levels.
_show(Tgen.plotTableSlices,
    x_cv        = "MixtureFraction",
    y_var       = "Heat_Release",
    slice_cv    = "EnthalpyTot",
    slice_range = None,
    n_slices    = 10,
    n_x_points  = 600,
    save_path   = _save("HRR_vs_Z_slices.png"))

# Visualize the table mesh and reaction rate at phi = 0.3.
level_index = _level_at(0.3)
print("Table level at phi = 0.30:", level_index)

# Mesh only:
_show(Tgen.visualizeTableLevel, level_index,
      save_path=_save("mesh_cH.png"))

# Mesh perimiter against the flamelet data it was built from:
_show(Tgen.visualizeTableLevelPerimiter, level_index,
      save_path=_save("mesh_boundary_data.png"))

# 2D colour map of Temperature (no wireframe):
_show(Tgen.visualizeTableLevel, level_index, "Temperature",
      plot_3d=False, show_grid=False, save_path=_save("T_cH_2d.png"))

# 3D surface of Temperature (with wireframe):
_show(Tgen.visualizeTableLevel, level_index, "Temperature",
      plot_3d=True, show_grid=True, save_path=_save("T_phi030_3d.png"))
# Production rate of progress variable:
_show(Tgen.visualizeTableLevel, level_index, "ProdRateTot_PV",
      plot_3d=False, show_grid=False, save_path=_save("ProdRate_phi030_2d.png"))
_show(Tgen.visualizeTableLevel, level_index, "ProdRateTot_PV",
      plot_3d=True, show_grid=False, save_path=_save("ProdRate_phi030_3d.png"))


# Visualize the table mesh and reaction rate at phi = 0.6.
level_index = _level_at(0.6)
print("Table level at phi = 0.60:", level_index)
# 3D surface of Temperature (with wireframe):
_show(Tgen.visualizeTableLevel, level_index, "Temperature",
      plot_3d=True, show_grid=True, save_path=_save("T_phi060_3d.png"))
# Production rate of progress variable:
_show(Tgen.visualizeTableLevel, level_index, "ProdRateTot_PV",
      plot_3d=False, show_grid=False, save_path=_save("ProdRate_phi060_2d.png"))
_show(Tgen.visualizeTableLevel, level_index, "ProdRateTot_PV",
      plot_3d=True, show_grid=False, save_path=_save("ProdRate_phi060_3d.png"))


# Visualize the table mesh and reaction rate at phi = 0.9.
level_index = _level_at(0.9)
print("Table level at phi = 0.90:", level_index)
# 3D surface of Temperature (with wireframe):
_show(Tgen.visualizeTableLevel, level_index, "Temperature",
      plot_3d=True, show_grid=True, save_path=_save("T_phi090_3d.png"))


# from 99% of the max progress variable, set the source terms of H2 to zero
# if the absolute value of the source terms is |S| < 0.1
Tgen.ClampSourceTerms(species_list=["H2", "H2O"], pv_frac=0.99, abs_tol=0.1)

# Export table in vtk format (view in ParaView).
Tgen.writeParaviewTable("LUT_paraview")

# Write SU2 .drg table file.
Tgen.writeSU2Table()

# Keep all figure windows open until the user closes them.
if SHOW_FIGURES:
    plt.show()
