# Collect flamelet data into data sets for table generation.
from Common.DataDrivenConfig import Config_FGM
from Data_Processing.collectFlameletData import FlameletConcatenator

Config = Config_FGM("TableGeneration.cfg")

Concat = FlameletConcatenator(Config)

# The (arbitrary, non-equal) equivalence-ratio bounds in the config exist only
# to keep MixtureFraction in the controlling variables (see 0_generate_config.py);
# they don't describe a real premixed sweep here. Without this, non-premixed
# (isPremixed()==False) flamelet data -- i.e. our counterflow flames -- gets
# silently clipped to the narrow mixture-fraction band those bounds imply,
# discarding almost the entire flame structure.
Concat.IgnoreMixtureBounds(True)

# Include NOx reaction rates and heat release in the flamelet data set.
Concat.SetAuxilarySpecies(["H2", "CO2", "H2O", "CO", "NOx"])
Concat.SetLookUpVars(["Heat_Release", "Density", "Y-OH", "X-H2", "X-CO2", "X-H2O", "X-CO"])

# Apply source term and chemical equilibrium data corrections for table generation.
Concat.WriteLUTData(True)

Concat.SetNEquilibriumNodes(1)   # sample N rows from each equilibrium file

# Read and concatenate flamelet data.
Concat.ConcatenateFlameletData()
