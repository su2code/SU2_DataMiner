#!/usr/bin/env python3

###############################################################################################
#       #      _____ __  _____      ____        __        __  ____                   #        #
#       #     / ___// / / /__ \    / __ \____ _/ /_____ _/  |/  (_)___  ___  _____   #        #
#       #     \__ \/ / / /__/ /   / / / / __ `/ __/ __ `/ /|_/ / / __ \/ _ \/ ___/   #        #
#       #    ___/ / /_/ // __/   / /_/ / /_/ / /_/ /_/ / /  / / / / / /  __/ /       #        #
#       #   /____/\____//____/  /_____/\__,_/\__/\__,_/_/  /_/_/_/ /_/\___/_/        #        #
#       #                                                                            #        #
###############################################################################################

############################ FILE NAME: generate_config.py ####################################
#=============================================================================================#
# author: Evert Bunschoten                                                                    |
#    :PhD Candidate ,                                                                         |
#    :Flight Power and Propulsion                                                             |
#    :TU Delft,                                                                               |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#   Generate SU2 DataMiner configuration for hydrogen tabulation test case.                   |
# Version: 3.1.0                                                                              |
#                                                                                             |
#=============================================================================================#

from su2dataminer.config import Config_FGM
import os 

Config = Config_FGM()
Config.SetConfigName("TableGeneration")

# Hydrogen-air flamelets with equivalence ratio between 0.3 and 1.0
Config.SetFuelDefinition(fuel_species=["H2"],fuel_weights=[1.0])
Config.SetReactionMechanism('h2o2.yaml')
Config.SetMixtureBounds(0.3, 1.0)
Config.SetNpMix(30)
Config.SetUnbTempBounds(300, 800)
Config.SetNpTemp(30)

# Enable preferential diffusion through selecting the "multicomponent" transport model.
Config.SetTransportModel('multicomponent')

Config.SetConcatenationFileHeader("LUT_data")

# Setting the progress variable definition.
Config.SetProgressVariableDefinition(pv_species=['H2', 'H', 'O2', 'O', 'H2O', 'OH', 'H2O2', 'HO2'],\
                                     pv_weights=[-2.59, 8.51e-02, -1.10e+00, -3.21e-01, +2.65e+00, -1.91e+00, +8.86e-02, +1.40e+00])

# Preparing flamelet output directory.
flamelet_data_dir = os.getcwd() + "/flamelet_data/"
if not os.path.isdir(flamelet_data_dir):
    os.mkdir(flamelet_data_dir)
Config.SetOutputDir(flamelet_data_dir) 

# Display settings and save configuration
Config.PrintBanner()
Config.SaveConfig()
