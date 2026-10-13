#!/usr/bin/env python3
# Regression test for the non-premixed FGM table: a 2D (MixtureFraction, EnthalpyTot) table built
# from counter-flow diffusion flamelets at a fixed global strain rate.
#
# The table plane is not spanned by the progress variable, so this exercises the branch that bounds
# a table level by the flamelet data alone rather than by evaluated reactant and product states.
# The manifold is read from LUT_data_full.csv, a reduced version of the CH4_fixed_strain test case.
import os
import sys

from su2dataminer.config import Config_FGM
from su2dataminer.manifold import SU2TableGenerator_FGM

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from drg_summary import writeTableSummary

config = Config_FGM()
config.SetFuelDefinition(fuel_species=["CH4", "N2"], fuel_weights=[0.23, 0.77])
config.SetOxidizerDefinition(oxidizer_species=["O2", "N2"], oxidizer_weights=[0.23, 0.77])
config.SetReactionMechanism('gri30.yaml')

# Unequal mixture bounds keep MixtureFraction among the controlling variables; they describe no
# premixed sweep here.
config.SetMixtureBounds(0.79, 0.81)
config.SetNpMix(1)
config.SetUnbTempBounds(300, 700)
config.SetProgressVariableDefinition(pv_species=['CO2', 'CO', 'H2', 'H2O'], pv_weights=[1, 1, 1, 1])
config.SetConcatenationFileHeader("LUT_data")
config.SetOutputDir(os.path.dirname(os.path.abspath(__file__)))

tablegen = SU2TableGenerator_FGM(config)

# Fixing the interpolator settings keeps the result reproducible and skips the tuning search.
tablegen.setNNearestNeighbors(8)
tablegen.setInverseDistanceExponent(2.0)

# The counter-flow state at fixed strain is a function of mixture fraction and enthalpy alone, so
# those two span the table and ProgressVariable is carried along as an ordinary table variable.
tablegen.setTableAxes(["MixtureFraction", "EnthalpyTot"])

# Deliberately coarse: the test exercises the code path, not the resolution, and a small table
# keeps the reference file small.
tablegen.setMaximumCellSize(8e-2)
tablegen.setPointCloudResolution(50)
tablegen.setBoundaryBinCount(16)
tablegen.setBoundaryTransverseBinCount(6)
tablegen.applyRefinementForGradientOf("Heat_Release", coef=0.5)

tablegen.generateTable()
tablegen.writeSU2Table("LUT_test")

# The harness cannot compare the tab separated data lines of the table file itself, so the
# extent of every table variable is written out alongside it.
writeTableSummary("LUT_test.drg", "LUT_summary.csv")
