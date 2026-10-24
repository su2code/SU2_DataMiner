#!/usr/bin/env python3
import os
import sys
import glob
import numpy as np
from su2dataminer.config import Config_FGM
from su2dataminer.generate_data import ComputeFlameletData,ComputeBoundaryData
from su2dataminer.process_data import FlameletConcatenator
from su2dataminer.manifold import TrainMLP_FGM

alpha_expo = -2.6
lr_decay = 0.999
activation_function="swish"
N_H = [10,10,10,10]
batch_expo = 2

config = Config_FGM()
config.SetFuelDefinition(["H2"], [1.0])
config.SetReactionMechanism("h2o2.yaml")
config.SetTransportModel('multicomponent')
config.SetNpConcatenation(244)
config.DefineMixtureStatus(False)
config.SetMixtureBounds(0.5, 0.6)
config.SetNpMix(2)
config.SetUnbTempBounds(300, 350)
config.SetNpTemp(6)
config.setFlameletTypes(["FREEFLAME","EQUILIBRIUM"])
config.AddOutputGroup(['Temperature'])
config.SetAlphaExpo(alpha_expo,0)
config.SetLRDecay(lr_decay,0)
config.SetActivationFunction(activation_function,0)
config.SetHiddenLayerArchitecture(N_H, 0)
config.SetBatchExpo(batch_expo,0)
config.SaveConfig()

ComputeFlameletData(config)
ComputeBoundaryData(config)


F = FlameletConcatenator(config, verbose_level=0)
F.ConcatenateFlameletData()
F.CollectBoundaryData()

T = TrainMLP_FGM(config, 0)
T.EnableBCLoss(False)
T.SetTrainHardware("CPU")
consistent_hp_params = (config.GetAlphaExpo(0)==T.alpha_expo==alpha_expo) and \
                 (config.GetLRDecay(0)==T.lr_decay==lr_decay) and \
                 (config.GetActivationFunction(0)==T.activation_function==activation_function) and \
                 (config.GetBatchExpo(0) == T.batch_expo == batch_expo)

T.SetNEpochs(100)
T.CommenceTraining()

config.UpdateMLPHyperParams(T)
config.PrintBanner()
    
if not consistent_hp_params:
    print("ERROR: network hyperparameters are inconsistent between config and trainer!")
    sys.exit(1)

# Weights and loss histories depend on the CPU, so only the final validation loss is compared.
history_files = glob.glob("architectures_Group1/Worker_0/Model_*/TrainingHistory.csv")
val_loss = np.loadtxt(max(history_files, key=os.path.getmtime), delimiter=",", skiprows=1)[-1, 2]
val_loss_ref = np.loadtxt("TrainingHistory_ref.csv", delimiter=",", skiprows=1)[-1, 2]
if not (0.5 < val_loss/val_loss_ref < 2.0):
    print("ERROR: final validation loss %.4e is not within a factor 2 of the reference %.4e" % (val_loss, val_loss_ref))
    sys.exit(1)
sys.exit(0)