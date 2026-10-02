#!/usr/bin/env python3
# Regression test for the premixed FGM table: a 2D (ProgressVariable, EnthalpyTot) table for a
# single equivalence ratio, built from free flames, burner-stabilized flames, their interpolated
# continuation and the equilibrium products.
#
# The table plane runs from the reactant to the product state, so this exercises the branch that
# bounds a table level by states evaluated directly, the equilibrium boundary at maximum progress
# variable, and the source term clamping. The manifold is read from LUT_data_full.csv, a reduced
# version of the CH4_phi080 test case.
import os
import sys

from su2dataminer.config import Config_FGM
from su2dataminer.manifold import SU2TableGenerator_FGM

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from drg_summary import writeTableSummary

config = Config_FGM()
config.SetFuelDefinition(fuel_species=["CH4"], fuel_weights=[1.0])
config.SetReactionMechanism('gri30.yaml')

# Equal mixture bounds mean a single mixture, described by progress variable and total enthalpy
# alone.
config.SetMixtureBounds(0.80, 0.80)
config.SetNpMix(1)
config.SetUnbTempBounds(300, 600)
config.SetProgressVariableDefinition(pv_species=['CO2', 'CO', 'H2', 'H2O'], pv_weights=[1, 1, 1, 1])
config.SetConcatenationFileHeader("LUT_data")
config.SetOutputDir(os.path.dirname(os.path.abspath(__file__)))

tablegen = SU2TableGenerator_FGM(config)

# Fixing the interpolator settings keeps the result reproducible and skips the tuning search.
tablegen.setNNearestNeighbors(8)
tablegen.setInverseDistanceExponent(2.0)

# Deliberately coarse: the test exercises the code path, not the resolution, and a small table
# keeps the reference file small.
tablegen.setMaximumCellSize(8e-2)
tablegen.setPointCloudResolution(50)
tablegen.setBoundaryBinCount(16)
tablegen.setBoundaryTransverseBinCount(6)
tablegen.applyRefinementForGradientOf("ProdRateTot_PV", coef=0.5)
tablegen.refineEquilibrium(coef=0.5, margin=2e-2)
tablegen.setSourceTermClamping(species_list=["CO", "CO2"], pv_margin=1e-2, abs_tol=0.1)

tablegen.generateTable()
tablegen.writeSU2Table("LUT_test")

# The harness cannot compare the tab separated data lines of the table file itself, so the
# extent of every table variable is written out alongside it.
writeTableSummary("LUT_test.drg", "LUT_summary.csv")
