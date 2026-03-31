from Common.DataDrivenConfig import Config_FGM 
import os 

Config = Config_FGM()
Config.SetConfigName("TableGeneration")

# Methane-air flamelets with equivalence ratio 0.75
Config.SetFuelDefinition(fuel_species=["CH4"],fuel_weights=[1.0])
Config.SetReactionMechanism('gri30.yaml')
Config.SetMixtureBounds(0.74, 0.76)
Config.SetNpMix(5)
Config.SetUnbTempBounds(275, 350)
Config.SetNpTemp(4)

# Enable preferential diffusion through selecting the "multicomponent" transport model.
#Config.SetTransportModel('multicomponent')
Config.SetTransportModel('mixture-averaged')

Config.SetConcatenationFileHeader("LUT_data")

# Setting the Efimov progress variable definition.
Config.SetProgressVariableDefinition(pv_species=['H2O', 'CO2', 'CO', 'H2'],\
                                     pv_weights=[1.0, 1.0, 1.0, 1.0])

# Preparing flamelet output directory.
flamelet_data_dir = os.getcwd() + "/flamelet_data/"
if not os.path.isdir(flamelet_data_dir):
    os.mkdir(flamelet_data_dir)
Config.SetOutputDir(flamelet_data_dir) 
 
Config.PrintBanner()
Config.SaveConfig()
