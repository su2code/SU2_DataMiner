# Generate flamelet data for pre-mixed hydrogen-air problems

from Common.DataDrivenConfig import FlameletAIConfig
from Data_Generation.DataGenerator_FGM import ComputeFlameletData

# Load FlameletAI configuration
Config = Config_FGM("TableGeneration.cfg")

# Distribute flamelet data generation process.
ComputeFlameletData(Config, run_parallel=True, N_processors=6)
