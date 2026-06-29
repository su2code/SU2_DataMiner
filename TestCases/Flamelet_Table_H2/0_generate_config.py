from Common.DataDrivenConfig import Config_FGM
import os

Config = Config_FGM()
Config.SetConfigName("TableGeneration")

# Hydrogen-air flamelets with equivalence ratio between 0.3 and 0.7
Config.SetFuelDefinition(fuel_species=["H2"],fuel_weights=[1.0])
Config.SetOxidizerDefinition(oxidizer_species=["O2", "N2"],
                             oxidizer_weights=[0.23, 0.77])
Config.SetReactionMechanism('h2o2.yaml')

Config.DefineMixtureStatus(False)  # Use equivalence ratio, not mixture fraction

Config.SetMixtureBounds(0.25, 1.25)
Config.SetNpMix(51)
Config.SetUnbTempBounds(250, 500)
Config.SetNpTemp(26)

# not in branch main
Config.SetNpMdot(20)          # burner flames across the mdot range
Config.SetMdotDHTarget(20000.0)   # J/kg target ΔH between flames
Config.SetNpMdotExtra(20)    # synthetic flames linearly interpolated from lowest-mdot burner flame to equilibrium
Config.SetInitialGridLength(0.2)  # Initial flamelet domain length in metres

# Explicitly select which flamelet types to generate.
Config.RunFreeFlames(True)
Config.RunBurnerFlames(True)
Config.RunExtraInterpolatedBurnerFlames(True)
Config.SetSrcInterpExponent(1.5)   # Decay of interpolated flamelets
Config.RunEquilibrium(True)



# Enable preferential diffusion through selecting the "multicomponent" transport model.
Config.SetTransportModel('mixture-averaged')
Config.SetConcatenationFileHeader("LUT_data")

# Setting the Efimov progress variable definition.
#Config.SetProgressVariableDefinition(pv_species=['H2', 'H', 'O2', 'O', 'H2O', 'OH', 'H2O2', 'HO2'],\
#                                     pv_weights=[ 0.0, 0.0,  0.0, 0.0,   1.0,  0.0,    0.0,  0.0])
Config.SetProgressVariableDefinition(pv_species=['H2', 'H', 'O2', 'O', 'H2O', 'OH', 'H2O2', 'HO2'],\
                                     pv_weights=[-7.36, -23.01, -2.04, -4.8, 1.83, -15.31, -57.02, 24.55])

# Preparing flamelet output directory.
flamelet_data_dir = os.getcwd() + "/flamelet_data/"
if not os.path.isdir(flamelet_data_dir):
    os.mkdir(flamelet_data_dir)
Config.SetOutputDir(flamelet_data_dir)

Config.PrintBanner()
Config.SaveConfig()
