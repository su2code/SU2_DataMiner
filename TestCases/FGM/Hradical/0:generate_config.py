import os 
from su2dataminer.config import Config_FGM

if not os.path.isdir("flamelet_data"):
    os.mkdir("flamelet_data")

pv_species = ["H2","H","O","O2","OH","H2O","HO2","H2O2"]
pv_weights = [-2.71,-8.54,-0.322, -0.783,-1.88,2.92,2.31,-6.36e-3]
c = Config_FGM()
c.SetFuelDefinition(["H2"],[1.0])
c.SetReactionMechanism("h2o2.yaml")
c.SetUnbTempBounds(300, 800)
c.SetControllingVariables(["ProgressVariable","EnthalpyTot","MixtureFraction","Y_H"])
c.RunFreeFlames(True)
c.RunBurnerFlames(True)
c.DefineMixtureStatus(False)
c.RunEquilibrium(True)
c.SetMixtureBounds(0.25, 20.0)
c.SetNpTemp(20)
c.SetNpMix(40)
c.SetTransportModel("multicomponent")
c.SetOutputDir("%s/flamelet_data/" % os.getcwd())
c.SetAverageLewisNumbers(0.5, 300.0)
c.SetProgressVariableDefinition(pv_species, pv_weights)
c.SetConfigName("config_Hradical")
c.PrintBanner()
c.SaveConfig()
