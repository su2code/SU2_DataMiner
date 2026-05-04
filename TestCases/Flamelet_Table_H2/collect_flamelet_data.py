#!/usr/bin/env python3

###############################################################################################
#       #      _____ __  _____      ____        __        __  ____                   #        #
#       #     / ___// / / /__ \    / __ \____ _/ /_____ _/  |/  (_)___  ___  _____   #        #
#       #     \__ \/ / / /__/ /   / / / / __ `/ __/ __ `/ /|_/ / / __ \/ _ \/ ___/   #        #
#       #    ___/ / /_/ // __/   / /_/ / /_/ / /_/ /_/ / /  / / / / / /  __/ /       #        #
#       #   /____/\____//____/  /_____/\__,_/\__/\__,_/_/  /_/_/_/ /_/\___/_/        #        #
#       #                                                                            #        #
###############################################################################################

########################## FILE NAME: collect_flamelet_data.py ################################
#=============================================================================================#
# author: Evert Bunschoten                                                                    |
#    :PhD Candidate ,                                                                         |
#    :Flight Power and Propulsion                                                             |
#    :TU Delft,                                                                               |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#   Accumulate flamelet data for hydrogen FGM tabulation test case.                           |
# Version: 3.1.0                                                                              |
#                                                                                             |
#=============================================================================================#
from su2dataminer.config import Config_FGM
from su2dataminer.process_data import FlameletConcatenator

Config = Config_FGM("TableGeneration.cfg")

# Interpolate 200 points per flamelet
Concat = FlameletConcatenator(Config)
Concat.SetNFlameletNodes(200)

# Read and concatenate flamelet data
Concat.ConcatenateFlameletData()
