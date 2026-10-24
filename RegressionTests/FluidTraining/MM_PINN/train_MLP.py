#!/usr/bin/env python3

import os
import sys
import glob
import numpy as np
from su2dataminer.manifold import TrainMLP_NICFD
from Common.DataDrivenConfig import Config_NICFD

C = Config_NICFD()
C.SetFluid("MM")
C.SetEquationOfState("HEOS")
C.UseAutoRange(True)
C.UsePTGrid(False)
C.SetNpDensity(50)
C.SetNpEnergy(50)
C.SetAlphaExpo(-1.8269e+00)
C.SetLRDecay(+9.8959e-01)
C.SetOutputDir(os.getcwd()+"/")
C.SetStateVars(["s","T","p","c2"])
C.SetHiddenLayerArchitecture([10])
C.SetActivationFunction("exponential")
C.SetBatchExpo(4)
C.SaveConfig()

M = TrainMLP_NICFD(C)
M.SetTrainHardware("CPU",0)
M.SetNEpochs(20)
M.SetScaler("minmax")
M.CommenceTraining()
C.UpdateMLPHyperParams(M)
C.WriteSU2MLP("SU2_MLP")


# Weights depend on the CPU, so only the final validation losses are compared.
loss_ref = np.array([1.5741e-01, 6.2190e-01, 1.0006e+00, 3.3100e-01])
history_files = glob.glob(M.worker_dir + "/Model_*/TrainingHistory.csv")
loss = np.loadtxt(max(history_files, key=os.path.getmtime), delimiter=",", skiprows=1)[-1, 2:]
if not np.all(np.isfinite(loss)) or np.any(loss > 1.5*loss_ref):
    print("ERROR: final validation losses", loss, "exceed the reference", loss_ref)
    sys.exit(1)
