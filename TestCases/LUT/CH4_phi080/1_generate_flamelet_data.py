# Generate flamelet data for a single phi = 0.80 methane-air flame

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

# refinement values:
# free_flame_refine={"ratio": 3.0, "slope": 0.1, "curve": 0.1, "prune":0.01},
# this leads to Np = 180

if __name__ == '__main__':
    ComputeFlameletData(Config, run_parallel=True, N_processors=4, loglevel=1,
                        free_flame_refine={"ratio": 3.0, "slope": 0.1, "curve": 0.1, "prune": 0.01},
                        burner_flame_refine={"ratio": 3.0, "slope": 0.15, "curve": 0.15, "prune": 0.01})
