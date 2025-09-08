
from su2dataminer.config import Config_FGM
from su2dataminer.generate_data import ComputeFlameletData,ComputeBoundaryData

c = Config_FGM("config_Hradical.cfg")

f = ComputeFlameletData(c, run_parallel=True,N_processors=4)
b = ComputeBoundaryData(c, run_parallel=True,N_processors=4)