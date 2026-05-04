# Collect flamelet data into data sets for table generation
from su2dataminer.config import Config_FGM
from su2dataminer.process_data import FlameletConcatenator

Config = Config_FGM("TableGeneration.cfg")

Concat = FlameletConcatenator(Config)
Concat.SetNFlameletNodes(200)

# Read and concatenate flamelet data
Concat.ConcatenateFlameletData()
