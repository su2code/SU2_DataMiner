import numpy as np 
from su2dataminer.config import Config_FGM 
from su2dataminer.process_data import FlameletConcatenator
c = Config_FGM("config_methane.cfg")

Fc = FlameletConcatenator(c)
Fc.SetLookUpVars(["Density","Heat_Release"])
Fc.ConcatenateFlameletData()
