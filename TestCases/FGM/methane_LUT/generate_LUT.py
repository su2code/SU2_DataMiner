from su2dataminer.config import Config_FGM
from su2dataminer.manifold import SU2TableGenerator
c = Config_FGM("config_methane.cfg")

tgen = SU2TableGenerator(c, n_near=6, p_fac=2)
tgen.SetBaseCellSize(7e-3)
tgen.VisualizeTableLevel()
tgen.VisualizeTableLevel(var_to_plot="Temperature")
tgen.GenerateTableNodes()
tgen.WriteTableFile()