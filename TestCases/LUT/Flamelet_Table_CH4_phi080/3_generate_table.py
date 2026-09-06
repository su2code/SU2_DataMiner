from Manifold_Generation.LUT.LUTGenerators import TableGenerator_FGM
from Common.DataDrivenConfig import Config_FGM

# Loading configuration.
Config = Config_FGM("TableGeneration.cfg")
Tgen = TableGenerator_FGM(Config)

# Iterate to achieve a target number of nodes
#Tgen.setTargetNodeCount(4000)
# Or manually specify the coarse cell size
Tgen.setMaximumCellSize(1.9e-2)

# Apply refinement based on the values of thermochemical quantities:
Tgen.applyRefinementWithin("Temperature", lowerbound=270, upperbound=450.0, coef=0.5)

# You can do this for any number of variables:
# Tgen.applyRefinementWithin("Cp", lowerbound=1000, upperbound=1200, coef=0.5)

# Scale refinement based on the gradients of quantities:
Tgen.applyRefinementForGradientOf("ProdRateTot_PV", coef=0.3)

# Apply refinement in proximity of the reactants and products
Tgen.refineEquilibrium(coef=0.5,margin=2e-2)

# Optionally: apply smoothing to table data to get rid of any waves or discontinuities
# Higher coefficient = more smoothing
#Tgen.setSmoothingParameter(0.1)

Tgen.generateTable()

# Export table in vtk format
Tgen.writeParaviewTable("LUT_vtk")


