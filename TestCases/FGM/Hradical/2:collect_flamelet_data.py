from su2dataminer.config import Config_FGM
from su2dataminer.process_data import FlameletConcatenator 

c = Config_FGM("config_Hradical.cfg")

F = FlameletConcatenator(c)
F.SetNFlameletNodes(2**6)
F.ConcatenateFlameletData()
F.CollectBoundaryData()

c.ClearOutputGroups()
c.AddOutputGroup(["Cp", "Beta_Enth_Thermal"])
c.AddOutputGroup(["Beta_Enth","Beta_MixFrac","Beta_ProgVar","MolarWeightMix"])
c.AddOutputGroup(["Conductivity","DiffusionCoefficient"])
c.AddOutputGroup(["ProdRateTot_PV","Heat_Release"])
c.AddOutputGroup(["Temperature","ViscosityDyn"])
c.AddOutputGroup(["Y_dot_net-H"])

c.SaveConfig()
