from Common.DataDrivenConfig import Config_FGM
import os

Config = Config_FGM()
Config.SetConfigName("TableGeneration")

# Fuel/oxidizer match equil.py: diluted CH4 (23% CH4 + 77% N2) vs diluted air
Config.SetFuelDefinition(fuel_species=["CH4", "N2"], fuel_weights=[0.23, 0.77])
Config.SetOxidizerDefinition(oxidizer_species=["O2", "N2"],
                             oxidizer_weights=[0.23, 0.77])
Config.SetReactionMechanism('gri30.yaml')

# For counterflow diffusion flames, use mixture fraction bounds covering full range
Config.DefineMixtureStatus(True)  # Use mixture fraction, not equivalence ratio
Config.SetMixtureBounds(0.0, 1.0)  # Full mixture fraction range for non-premixed flames
Config.SetNpMix(2)  # Minimum value to avoid empty range (not used for counterflow flames)

Config.RunFreeFlames(False)
Config.SetUnbTempBounds(250, 800)
Config.SetNpTemp(56)

Config.RunBurnerFlames(False)
Config.RunEquilibrium(False)

Config.RunCounterFlames(True)
Config.SetCounterFlowFixedStrain(True)
Config.SetCounterFlowStrainRate(56.0)
Config.SetInitialGridLength(0.02)   # flame width in metres
Config.SetSaveMoleFractions(True)   # Save mole fractions (X-) alongside mass fractions (Y-)

Config.SetTransportModel('unity-Lewis-number')
Config.SetConcatenationFileHeader("LUT_data")


# progress variable definition (still needed, not used)
Config.SetProgressVariableDefinition(
    pv_species=['CO2', 'CO','H2','H2O'],
    pv_weights=[1, 1, 1, 1])

# Preparing flamelet output directory.
flamelet_data_dir = os.getcwd() + "/flamelet_data/"
if not os.path.isdir(flamelet_data_dir):
    os.mkdir(flamelet_data_dir)
Config.SetOutputDir(flamelet_data_dir)

Config.PrintBanner()
Config.SaveConfig()
