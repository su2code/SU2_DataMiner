from su2dataminer.config import Config_FGM
from su2dataminer.manifold import TrainMLP_FGM

c = Config_FGM("config_Hradical.cfg")

N_H = [[18,8,8],[12,26],[16,12,8,5],[22, 39, 34, 19],[12,13,9],[18, 34, 29, 14]]
alpha_expo = [-1.785,-1.778,-1.656,-2.779,-1.845,-2.279]
lr_decay = [0.9906,0.9911,0.9916,0.9918,0.9906,0.9918]
phi = "swish"
batch_expo = 6 

for iGroup in range(c.GetNMLPOutputGroups()):
    trainer = TrainMLP_FGM(c, iGroup)
    trainer.SetActivationFunction(phi)
    trainer.SetBatchExpo(batch_expo)
    trainer.SetHiddenLayers(N_H[iGroup])
    trainer.SetAlphaExpo(alpha_expo[iGroup])
    trainer.SetLRDecay(lr_decay[iGroup])
    trainer.SetVerbose(1)
    trainer.CommenceTraining()
    c.UpdateMLPHyperParams(trainer)
    c.SaveConfig()