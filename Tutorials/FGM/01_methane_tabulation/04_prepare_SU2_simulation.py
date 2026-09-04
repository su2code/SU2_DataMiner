import gmsh
import numpy as np
from su2dataminer.config import Config_FGM
from Data_Generation.FlameletSolvers import FreeFlameSolver

# Reactant temperature and equivalence ratio for simulation
Tu = 350
eq_ratio = 0.7

# Run an adiabatic flamelet simulation to estimate adiabatic flame speed and grid size.
config = Config_FGM("methane_tabulation.cfg")
freeflame = FreeFlameSolver(config)
freeflame.setReactantTemperature(Tu)
freeflame.setMixtureStatus(eq_ratio)
freeflame.startSolver()
solution = freeflame.getSolution()

flame_velocity = solution["Velocity"].iloc[0]
grid = np.array(solution["Distance"])
dx = grid[1:]-grid[:-1]
dx_min = min(dx)
T = np.array(solution["Temperature"])
dTdx = (T[1:] - T[:-1])/dx
max_dTdx = max(dTdx)
t_flame = (max(T) - min(T))/max_dTdx

# Generate mesh
gmsh.initialize()
gmsh.model.add("2D flame domain")
gmsh.option.setNumber("General.Verbosity", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints",0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary",0)
gmsh.option.setNumber("Mesh.MeshSizeMax",10*dx_min)
gmsh.option.setNumber("Mesh.Algorithm", 8)
gmsh.option.setNumber("Mesh.SaveAll",0)
gmsh_geo = gmsh.model.geo
gmsh_mesher = gmsh.model.mesh

x_inlet = -4*t_flame
x_outlet = 10*t_flame
w_domain = 2*t_flame
dx_burner = 1.5*t_flame

pts_inlet = [gmsh_geo.addPoint(x_inlet, 0, 0), gmsh_geo.addPoint(x_inlet, w_domain, 0)]

pts_burner_wall = [gmsh_geo.addPoint(-dx_burner, w_domain, 0),\
                   gmsh_geo.addPoint(-dx_burner, 0.5*w_domain, 0),\
                   gmsh_geo.addPoint(0, 0.5*w_domain, 0),\
                   gmsh_geo.addPoint(0, w_domain, 0)]

pts_outflow = [gmsh_geo.addPoint(x_outlet, w_domain, 0), gmsh_geo.addPoint(x_outlet, 0, 0)]

line_inlet = gmsh_geo.addLine(pts_inlet[0],pts_inlet[1])
lines_burner_wall = [gmsh_geo.addLine(pts_burner_wall[0],pts_burner_wall[1]),\
                     gmsh_geo.addLine(pts_burner_wall[1],pts_burner_wall[2]),\
                     gmsh_geo.addLine(pts_burner_wall[2],pts_burner_wall[3])]

line_outlet = gmsh_geo.addLine(pts_outflow[0],pts_outflow[1])
lines_sym = [gmsh_geo.addLine(pts_inlet[1], pts_burner_wall[0]),\
             gmsh_geo.addLine(pts_burner_wall[-1], pts_outflow[0]),\
             gmsh_geo.addLine(pts_outflow[1], pts_inlet[0])]

crvloop = gmsh_geo.addCurveLoop([line_inlet] + lines_sym + lines_burner_wall+[line_outlet],reorient=True)
fluid_plane = gmsh_geo.addPlaneSurface([crvloop])

gmsh_geo.addPhysicalGroup(1, [line_inlet],name="inlet")
gmsh_geo.addPhysicalGroup(1, [line_outlet],name="outlet")
gmsh_geo.addPhysicalGroup(1, lines_burner_wall,name="burner_wall")
gmsh_geo.addPhysicalGroup(1, lines_sym,name="symmetry")
gmsh_geo.addPhysicalGroup(2, [fluid_plane],name="fluid")

gmsh_geo.synchronize()
gmsh_mesher.setRecombine(2, fluid_plane)
gmsh_mesher.generate(2)
gmsh.write("su2mesh.su2")
gmsh.fltk.run()
gmsh.finalize()

# Prepare SU2 configuration file
SU2ConfigOptions = """SOLVER = INC_NAVIER_STOKES
KIND_TURB_MODEL= NONE
RESTART_SOL= NO
INC_DENSITY_MODEL= VARIABLE
INC_ENERGY_EQUATION = YES
INC_VELOCITY_INIT= (__U_INIT__, 0.0, 0.0 )
INC_TEMPERATURE_INIT= __T_U__
THERMODYNAMIC_PRESSURE= 101325
INC_NONDIM= DIMENSIONAL

FLUID_MODEL= FLUID_FLAMELET
VISCOSITY_MODEL= FLAMELET
CONDUCTIVITY_MODEL= FLAMELET
DIFFUSIVITY_MODEL= FLAMELET
KIND_SCALAR_MODEL= FLAMELET
INTERPOLATION_METHOD= LUT
FILENAMES_INTERPOLATOR= (LUT_methane.drg)
PREFERENTIAL_DIFFUSION= NO

FLAME_INIT_METHOD= SPARK
SPARK_INIT= (__X_SPARK__, __Y_SPARK__, 0, __R_SPARK__, 300, 50)

SPARK_REACTION_RATES= (200, 0, 0, 0)

SPECIES_INIT = (__PV_INIT__, __ENTH_INIT__, 0, 0)

% Passive reactants in flamelet problem
SPECIES_CLIPPING= YES
SPECIES_CLIPPING_MAX= 1.0, 1e7, 1.0, 1.0
SPECIES_CLIPPING_MIN= -1.0, -1e7, 0.0, 0.0

CONTROLLING_VARIABLE_NAMES= (ProgressVariable, EnthalpyTot)
CONTROLLING_VARIABLE_SOURCE_NAMES= (ProdRateTot_PV, NULL)
USER_SCALAR_NAMES= (Y-CO,Y-CO2)
USER_SOURCE_NAMES= ( Y_dot_pos-CO,  Y_dot_neg-CO, Y_dot_net-CO2, NULL)

MARKER_INLET_SPECIES = (inlet, __PV_INIT__, __ENTH_INIT__, 0, 0)

MARKER_ISOTHERMAL= (burner_wall, __T_U__)
MARKER_SPECIES_STRONG_BC=(burner_wall)
MARKER_SYM= (symmetry)
INC_INLET_TYPE= VELOCITY_INLET
MARKER_INLET = (inlet, __T_U__, __U_INLET__, 1.0, 0.0, 0.0)
INC_OUTLET_TYPE= PRESSURE_OUTLET
MARKER_OUTLET= (outlet, 0.0)
MARKER_ANALYZE_AVERAGE = AREA

NUM_METHOD_GRAD= GREEN_GAUSS
CFL_NUMBER= 200
CFL_ADAPT= NO
ITER=4000
OUTPUT_WRT_FREQ= 20

LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-4
LINEAR_SOLVER_ITER=20

CONV_NUM_METHOD_FLOW= FDS
CONV_NUM_METHOD_SPECIES= BOUNDED_SCALAR
MUSCL_FLOW= YES
MUSCL_SPECIES= YES
SLOPE_LIMITER_FLOW = NONE
SLOPE_LIMITER_SPECIES= NONE
TIME_DISCRE_FLOW= EULER_IMPLICIT
TIME_DISCRE_SPECIES= EULER_IMPLICIT

CONV_FIELD = RMS_EnthalpyTot
CONV_RESIDUAL_MINVAL= -10
CONV_STARTITER= 20
SCREEN_OUTPUT = INNER_ITER RMS_PRESSURE RMS_ProgressVariable RMS_EnthalpyTot
HISTORY_OUTPUT = WALL_TIME RMS_RES
VOLUME_OUTPUT = SOLUTION, PRIMITIVE, SOURCE, LOOKUP

MESH_FORMAT= SU2
MESH_FILENAME = su2mesh.su2
OUTPUT_FILES = (RESTART,PARAVIEW)
TABULAR_FORMAT = CSV
CONV_FILENAME= history
VOLUME_FILENAME= flow
"""

pv_init, enth_init,_ = config.GetUnburntScalars(eq_ratio, Tu)
SU2ConfigOptions = SU2ConfigOptions.replace("__U_INIT__", "%.3e" % (2*flame_velocity))
SU2ConfigOptions = SU2ConfigOptions.replace("__U_INLET__", "%.3e" % (flame_velocity))
SU2ConfigOptions = SU2ConfigOptions.replace("__T_U__", "%.3e" % Tu)
SU2ConfigOptions = SU2ConfigOptions.replace("__X_SPARK__", "%.3e" % (0.5*x_outlet))
SU2ConfigOptions = SU2ConfigOptions.replace("__Y_SPARK__", "%.3e" % (0.5*w_domain))
SU2ConfigOptions = SU2ConfigOptions.replace("__R_SPARK__", "%.3e" % (0.5*w_domain))
SU2ConfigOptions = SU2ConfigOptions.replace("__PV_INIT__", "%.3e" % pv_init)
SU2ConfigOptions = SU2ConfigOptions.replace("__ENTH_INIT__", "%.3e" % enth_init)

with open("config_SU2.cfg","w+") as fid:
    fid.write(SU2ConfigOptions)
