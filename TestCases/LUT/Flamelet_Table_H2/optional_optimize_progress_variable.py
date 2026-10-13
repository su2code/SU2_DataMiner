# Optional: replace the progress variable definition of 0_generate_config.py by one optimized for
# this manifold. Run it between 1_generate_flamelet_data.py and 2_collect_flamelet_data.py, which is
# where the progress variable is evaluated from the flamelet mass fractions. Not part of the four
# step sequence: 0_generate_config.py sets the Efimov definition instead.
from su2dataminer.config import Config_FGM
from su2dataminer.process_data import PVOptimizer

config = Config_FGM("TableGeneration.cfg")

PVO = PVOptimizer(config)
PVO.SetAdditionalProgressVariables(["Temperature","Heat_Release"])
PVO.SetNGenerations(30)

PVO.OptimizePV()

config.SetProgressVariableDefinition(PVO.GetOptimizedSpecies(), PVO.GetOptimizedWeights())
config.PrintBanner()
config.SaveConfig()
