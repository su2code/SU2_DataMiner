# Plot the flamelet data together with the table level boundary that bounds it, for inspection.
#
# The boundary is the perimeter the mesher is given, so this shows exactly which part of the state
# space ends up in the table: flamelet data outside the perimeter is lost, and perimeter enclosing
# no data is extrapolated into.
import os
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path

from Common.DataDrivenConfig import Config_FGM
from Manifold_Generation.LUT.LUTGenerators import SU2TableGenerator_FGM

Config = Config_FGM("TableGeneration.cfg")

# Match the settings of 3_generate_table.py: the boundary depends on the coarse cell size.
Tgen = SU2TableGenerator_FGM(Config)
Tgen.setMaximumCellSize(1.5e-2)
Tgen.setPointCloudResolution(400)
Tgen.setBoundaryBinCount(100)
Tgen.setBoundaryTransverseBinCount(40)

# Prepare the table levels and the controlling variable scaling without generating the table itself.
Tgen._processTableLevels()
Tgen._defineFluidDataInterpolator()

controlling_variables = Config.GetControllingVariables()
flamelet_data = Tgen._getFluidDataForInterpolator()
cv_data = np.column_stack(tuple(flamelet_data[cv].to_numpy() for cv in controlling_variables))
cv_scaled = Tgen._scaler_controlling_variables.transform(cv_data)[:, :2]

for iLevel, level in enumerate(Tgen._table_levels):
    boundary = Tgen._createBoundaryPolylineForTableLevel(level)
    closed_boundary = np.vstack((boundary, boundary[:1]))

    perimeter = Path(boundary)
    enclosed = perimeter.contains_points(cv_scaled, radius=1e-9) | \
               perimeter.contains_points(cv_scaled, radius=-1e-9)
    area = 0.5*abs(np.sum(boundary[:, 0]*np.roll(boundary[:, 1], -1) -
                          boundary[:, 1]*np.roll(boundary[:, 0], -1)))

    print("Table level %i (%s = %.5f)" % (iLevel, controlling_variables[-1], level))
    print("  boundary vertices     : %i" % len(boundary))
    print("  flamelet data points  : %i" % len(cv_scaled))
    print("  points outside        : %i (%.3f%%)" % (np.count_nonzero(~enclosed),
                                                     100*np.mean(~enclosed)))
    print("  enclosed area         : %.5f of the scaled state space" % area)

    views = [("full", None, None),
             ("unburnt side (minimum progress variable)", (-0.02, 0.45), (0.50, 0.90)),
             ("burnt side (maximum progress variable)", (0.955, 1.01), (0.45, 1.02))]

    fig, axes = plt.subplots(1, len(views), figsize=(6.2*len(views), 6.5))
    for ax, (title, xlimits, ylimits) in zip(axes, views):
        ax.plot(cv_scaled[enclosed, 0], cv_scaled[enclosed, 1], '.', ms=1.4,
                color='0.72', label='flamelet data (enclosed)')
        if np.any(~enclosed):
            ax.plot(cv_scaled[~enclosed, 0], cv_scaled[~enclosed, 1], '.', ms=5,
                    color='red', label='flamelet data (outside)')
        ax.plot(closed_boundary[:, 0], closed_boundary[:, 1], '-', color='tab:blue', lw=1.6,
                label='table level boundary')
        ax.plot(boundary[:, 0], boundary[:, 1], 'o', ms=3.2, color='tab:blue',
                label='boundary vertices')
        if xlimits:
            ax.set_xlim(xlimits)
            ax.set_ylim(ylimits)
        ax.set_xlabel("%s (scaled)" % controlling_variables[0], fontsize=11)
        ax.set_ylabel("%s (scaled)" % controlling_variables[1], fontsize=11)
        ax.set_title(title, fontsize=11)
    axes[0].legend(markerscale=6, fontsize=9, loc='lower left')

    figure_name = "boundary_check.png" if len(Tgen._table_levels) == 1 \
                  else "boundary_check_level_%02i.png" % iLevel
    fig.tight_layout()
    fig.savefig(figure_name, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print("  Saved: %s" % figure_name)
