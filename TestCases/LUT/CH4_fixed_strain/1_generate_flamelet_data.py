# Generate flamelet data for a counterflow diffusion flame sweep

# Limit inner thread pools BEFORE any library imports to prevent oversubscription
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

from Common.DataDrivenConfig import Config_FGM
from Data_Generation.DataGenerator_FGM import ComputeFlameletData

# Load FGM configuration
Config = Config_FGM("TableGeneration.cfg")

ComputeFlameletData(Config, run_parallel=False, N_processors=1, loglevel=1,
                    counter_flame_refine={"ratio": 3, "slope": 0.05, "curve": 0.05, "prune": 0.01})

