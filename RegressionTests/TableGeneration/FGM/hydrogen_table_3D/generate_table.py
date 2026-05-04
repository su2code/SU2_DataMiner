#!/usr/bin/env python3
import os 
from su2dataminer.config import Config_FGM 
from su2dataminer.generate_data import ComputeFlameletData
from su2dataminer.process_data import FlameletConcatenator
from su2dataminer.manifold import SU2TableGenerator_FGM

n_cores = 4
# Create configuration
Config = Config_FGM()
# Hydrogen-air flamelets with equivalence ratio between 0.3 and 0.7
Config.SetFuelDefinition(fuel_species=["H2"],fuel_weights=[1.0])
Config.SetReactionMechanism('h2o2.yaml')
Config.SetMixtureBounds(0.4, 0.7)
Config.SetNpMix(10)
Config.SetUnbTempBounds(300, 400)
Config.SetNpTemp(20)

# Enable preferential diffusion through selecting the "multicomponent" transport model.
Config.SetTransportModel('multicomponent')

Config.SetConcatenationFileHeader("LUT_data")

# Setting the Efimov progress variable definition.
Config.SetProgressVariableDefinition(pv_species=['H2', 'H', 'O2', 'O', 'H2O', 'OH', 'H2O2', 'HO2'],\
                                     pv_weights=[-2.59, 8.51e-02, -1.10e+00, -3.21e-01, +2.65e+00, -1.91e+00, +8.86e-02, +1.40e+00])

# Preparing flamelet output directory.
flamelet_data_dir = os.getcwd() + "/flamelet_data/"
if not os.path.isdir(flamelet_data_dir):
    os.mkdir(flamelet_data_dir)
Config.SetOutputDir(flamelet_data_dir) 
Config.SetConfigName("LUT_test")
Config.SaveConfig()

# Generate flamelet data 
ComputeFlameletData(Config, run_parallel=True, N_processors=n_cores)

# Accumulate flamelet data 
FC = FlameletConcatenator(Config, verbose_level=0)
FC.SetNFlameletNodes(100)
FC.ConcatenateFlameletData()

# Table generator 
Tgen = SU2TableGenerator_FGM(Config)
Tgen.SetNCores(n_cores)
# Manually set mixture fraction limits.
Tgen.SetMixtureFractionLimits(mix_frac_max=0.02, mix_frac_min=0.00939225575395504)
Tgen.SetNTableLevels(10)
Tgen.SetTableVars(["Temperature","ProdRateTot_PV"])
Tgen.SetNnodes_Target(600)
# Visualize the interpolated reaction rate at equivalence ratio 0.5.
cv_target = Config.GetUnburntScalars(equivalence_ratio=0.5, temperature=300.0)
pv_target = cv_target[0]
z_target = cv_target[2]
# # Visualize table level connectivity at equivalence ratio 0.5.
Tgen.InsertMixtureFractionLevel(z_target)
# Visualize table level connectivity at equivalence ratio 0.5.

# # Generate table connectivity and interpolate flamelet data.
Tgen.GenerateTableNodes()

Tgen.WriteTableFile("LUT_test.drg")