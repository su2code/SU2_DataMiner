from Common.DataDrivenConfig import Config_FGM
import os

Config = Config_FGM()
Config.SetConfigName("TableGeneration")

# Fuel/oxidizer for a diluted CH4/air counterflow diffusion flame.
Config.SetFuelDefinition(fuel_species=["CH4", "N2"], fuel_weights=[0.23, 0.77])
Config.SetOxidizerDefinition(oxidizer_species=["O2", "N2"],
                             oxidizer_weights=[0.23, 0.77])
Config.SetReactionMechanism('gri30.yaml')

# Mixture bounds only affect premixed flamelet types, all disabled below.
# NOTE: bounds must NOT be equal here. Config_FGM derives the controlling
# variable list from these bounds, and equal bounds would drop MixtureFraction
# from that list -- which the non-premixed (Z, h) table needs as a plane axis.
Config.SetMixtureBounds(0.79, 0.81)
Config.SetNpMix(1)

Config.SetUnbTempBounds(250, 800)
Config.SetNpTemp(56)

Config.RunFreeFlames(False)
Config.RunBurnerFlames(False)
Config.RunEquilibrium(False)

Config.RunCounterFlames(True)
Config.SetCounterFlowFixedStrain(True)
Config.SetCounterFlowStrainRate(56.0)
Config.SetInitialGridLength(0.02)   # flame width in metres

Config.SetTransportModel('unity-Lewis-number')
Config.SetConcatenationFileHeader("LUT_data")
Config.SetSaveMoleFractions(True)   # needed for the X-species lookup variables in step 2

# Progress variable definition: required so flamelet data collection can compute
# a ProgressVariable column, even though it plays no physical role here -- the
# (Z, h) table uses ProgressVariable only as an unused, fixed "level" CV.
Config.SetProgressVariableDefinition(
    pv_species=['CO2', 'CO', 'H2', 'H2O'],
    pv_weights=[1, 1, 1, 1])

# Preparing flamelet output directory.
flamelet_data_dir = os.getcwd() + os.sep + "flamelet_data"
if not os.path.isdir(flamelet_data_dir):
    os.mkdir(flamelet_data_dir)
Config.SetOutputDir(flamelet_data_dir)

Config.PrintBanner()
Config.SaveConfig()
