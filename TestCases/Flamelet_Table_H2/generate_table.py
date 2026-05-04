#!/usr/bin/env python3

###############################################################################################
#       #      _____ __  _____      ____        __        __  ____                   #        #
#       #     / ___// / / /__ \    / __ \____ _/ /_____ _/  |/  (_)___  ___  _____   #        #
#       #     \__ \/ / / /__/ /   / / / / __ `/ __/ __ `/ /|_/ / / __ \/ _ \/ ___/   #        #
#       #    ___/ / /_/ // __/   / /_/ / /_/ / /_/ /_/ / /  / / / / / /  __/ /       #        #
#       #   /____/\____//____/  /_____/\__,_/\__/\__,_/_/  /_/_/_/ /_/\___/_/        #        #
#       #                                                                            #        #
###############################################################################################

############################# FILE NAME: generate_table.py ####################################
#=============================================================================================#
# author: Evert Bunschoten                                                                    |
#    :PhD Candidate ,                                                                         |
#    :Flight Power and Propulsion                                                             |
#    :TU Delft,                                                                               |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#   Generate 3D table for hydrogen FGM tabulation test case.                                  |
# Version: 3.1.0                                                                              |
#                                                                                             |
#=============================================================================================#
from su2dataminer.manifold import SU2TableGenerator_FGM
from su2dataminer.config import Config_FGM

# Loading configuration.
Config = Config_FGM("TableGeneration.cfg")

# Initializing table module and pre-process interpolator.
Tgen = SU2TableGenerator_FGM(Config)

# Distribute tabulation process over 4 cores.
Tgen.SetNCores(4)

# Manually set mixture fraction limits.
Tgen.SetMixtureFractionLimits(mix_frac_max=0.02, mix_frac_min=0.009392)

# Use 50 table levels and approximately 2k nodes per table level.
Tgen.SetNTableLevels(50)
Tgen.SetNnodes_Target(2000)

# Insert mixture fraction level at equivalence ratio of 0.5.
cv_target = Config.GetUnburntScalars(equivalence_ratio=0.5, temperature=300.0)
z_target = cv_target[2]
Tgen.InsertMixtureFractionLevel(z_target)

# Visualize table level connectivity at equivalence ratio 0.5.
Tgen.VisualizeTableLevel(z_target)

# Visualize table data for temperature at equivalence ratio 0.5
Tgen.VisualizeTableLevel(z_target, var_to_plot="Temperature")

# Generate table connectivity and interpolate flamelet data for each table level.
Tgen.GenerateTableNodes()

# Save table levels in vtk format so table data can be inspected with ParaView
Tgen.WriteOutParaview("hydrogen_table")

# Write SU2 .drg table file.
Tgen.WriteTableFile("hydrogen_table")