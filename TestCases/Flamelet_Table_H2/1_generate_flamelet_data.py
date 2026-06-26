# Generate flamelet data for a premixed hydrogen flame sweep

# Limit inner thread pools BEFORE any library imports to prevent oversubscription
import os

# 8 cores, only one thread.
MPI_N = 8
MPI_NT = 1

os.environ.setdefault("OMP_NUM_THREADS", str(MPI_NT))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(MPI_NT))
os.environ.setdefault("MKL_NUM_THREADS", str(MPI_NT))
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", str(MPI_NT))
os.environ.setdefault("NUMEXPR_NUM_THREADS", str(MPI_NT))

from Common.DataDrivenConfig import Config_FGM
from Data_Generation.DataGenerator_FGM import ComputeFlameletData

# Load FGM configuration
Config = Config_FGM("TableGeneration.cfg")

# refinement values:
# free_flame_refine={"ratio": 3.0, "slope": 0.1, "curve": 0.1, "prune":0.01},
# this leads to Np = 180 (small mesh)



ComputeFlameletData(Config, run_parallel=True, N_processors=MPI_N, loglevel=0,
                    # medium
                    free_flame_refine={"ratio": 2.0, "slope": 0.02, "curve": 0.02, "prune":0.005},
                    burner_flame_refine={"ratio": 3.0, "slope": 0.02, "curve": 0.02, "prune": 0.01})

