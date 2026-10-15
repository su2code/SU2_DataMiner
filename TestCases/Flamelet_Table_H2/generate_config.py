
from Common.DataDrivenConfig import FlameletAIConfig
import os

Config = Config_FGM()
Config.SetConfigName("TableGeneration")

# Hydrogen-air flamelets with equivalence ratio between 0.3 and 1.0
Config.SetFuelDefinition(fuel_species=["H2"],fuel_weights=[1.0])
Config.SetReactionMechanism('h2o2.yaml')
Config.SetMixtureBounds(0.25, 2.0)
Config.SetNpMix(100)
Config.SetUnbTempBounds(300, 800)
Config.SetNpTemp(100)

# Enable preferential diffusion through selecting the "multicomponent" transport model.
Config.SetTransportModel('multicomponent')

# Preferential diffusion formulation. SOURCE_TERM selects model B2, which also switches the mixture
# fraction to the elemental (hydrogen) definition of Eq. (11). BETA_CORRECTION keeps the legacy
# beta-scalar path and Bilger's mixture fraction.
Config.SetPreferentialDiffusionMethod("SOURCE_TERM")

# Major species of FGM model B2 (Schepers & van Oijen, C&F 280 (2025) 114332). Their Eq. (14)
# molecular and thermal flux coefficients are tabulated per species and summed by the CFD solver
# against its own gradients; only the remaining species are gathered into the Eq. (16) closure
# sources S_PV, S_Z and S_h. Leaving this unset puts every species into the closure sources, which
# pre-contracts the summation against the 1D flamelet gradients and loses all curvature and
# cross-diffusion response in the CFD.
Config.SetPDMajorSpecies(["H2", "H2O", "H"])

Config.SetConcatenationFileHeader("LUT_data")

# Setting the progress variable definition.
Config.SetProgressVariableDefinition(pv_species=['H2', 'H', 'O2', 'O', 'H2O', 'OH', 'H2O2', 'HO2'],\
                                     pv_weights=[-2.59, 8.51e-02, -1.10e+00, -3.21e-01, +2.65e+00, -1.91e+00, +8.86e-02, +1.40e+00])

# Preparing flamelet output directory.
flamelet_data_dir = os.getcwd() + "/flamelet_data/"
if not os.path.isdir(flamelet_data_dir):
    os.mkdir(flamelet_data_dir)
Config.SetOutputDir(flamelet_data_dir)

Config.PrintBanner()
Config.SaveConfig()
