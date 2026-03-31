# Generate flamelet data for pre-mixed methane-air problems

# Limit inner thread pools BEFORE any library imports to prevent oversubscription
# when joblib spawns parallel workers (OpenBLAS pthreads + MKL + llvm-openmp
# each default to all cores, multiplying with N_processors -> OOM kill)
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

# Distribute flamelet data generation process over 20 cores.
ComputeFlameletData(Config, run_parallel=True, N_processors=2)
