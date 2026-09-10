from Common.DataDrivenConfig import Config_FGM
import os

Config = Config_FGM()
Config.SetConfigName("TableGeneration")

# Methane-air flamelets at a single equivalence ratio phi = 0.80
Config.SetFuelDefinition(fuel_species=["CH4"], fuel_weights=[1.0])
Config.SetReactionMechanism('gri30.yaml')

# Single phi = 0.80 flame -> single mixture fraction point (min == max). Equal mixture bounds
# reduce the controlling variables to progress variable and total enthalpy, which is exactly the
# plane this table is spanned by, so no table axes have to be selected in 3_generate_table.py.
Config.SetMixtureBounds(0.80, 0.80)
Config.SetNpMix(1)

Config.SetUnbTempBounds(270, 750)
Config.SetNpTemp(49)
Config.SetNpMdot(40)          # burner flames across the mdot range
Config.SetMdotDHTarget(20000.0)   # J/kg target dH between flames
Config.SetNpMdotExtra(40)    # synthetic flames linearly interpolated from lowest-mdot burner flame to equilibrium
Config.SetInitialGridLength(0.2)  # Initial flamelet domain length in metres

# Explicitly select which flamelet types to generate.
Config.RunFreeFlames(True)
Config.RunBurnerFlames(True)
Config.RunExtraInterpolatedBurnerFlames(True)
Config.SetSrcInterpExponent(2.0)   # Decay of interpolated flamelets

Config.RunEquilibrium(True)

Config.SetTransportModel('unity-Lewis-number')

Config.SetConcatenationFileHeader("LUT_data")


# progress variable definition
Config.SetProgressVariableDefinition(
    pv_species=['CO2', 'CO', 'H2', 'H2O'],
    pv_weights=[1, 1, 1, 1])

# Preparing flamelet output directory.
flamelet_data_dir = os.getcwd() + "/flamelet_data/"
if not os.path.isdir(flamelet_data_dir):
    os.mkdir(flamelet_data_dir)
Config.SetOutputDir(flamelet_data_dir)

Config.PrintBanner()
Config.SaveConfig()
