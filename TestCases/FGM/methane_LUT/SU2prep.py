from su2dataminer.config import Config_FGM
c = Config_FGM("config_methane.cfg")

print(c.GetUnburntScalars(0.5, 300.0))