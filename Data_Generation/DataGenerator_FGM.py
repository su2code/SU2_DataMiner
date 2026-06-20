###############################################################################################
#       #      _____ __  _____      ____        __        __  ____                   #        #
#       #     / ___// / / /__ \    / __ \____ _/ /_____ _/  |/  (_)___  ___  _____   #        #
#       #     \__ \/ / / /__/ /   / / / / __ `/ __/ __ `/ /|_/ / / __ \/ _ \/ ___/   #        #
#       #    ___/ / /_/ // __/   / /_/ / /_/ / /_/ /_/ / /  / / / / / /  __/ /       #        #
#       #   /____/\____//____/  /_____/\__,_/\__/\__,_/_/  /_/_/_/ /_/\___/_/        #        #
#       #                                                                            #        #
###############################################################################################

############################ FILE NAME: DataGenerator_NICFD.py ################################
#=============================================================================================#
# author: Evert Bunschoten                                                                    |
#    :PhD Candidate ,                                                                         |
#    :Flight Power and Propulsion                                                             |
#    :TU Delft,                                                                               |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#  Class for generating fluid data for Flamelet-Generated Manifold data mining operations.    |
#                                                                                             |
# Version: 3.1.0                                                                              |
#                                                                                             |
#=============================================================================================#

#---------------------------------------------------------------------------------------------#
# Importing general packages
#---------------------------------------------------------------------------------------------#
import cantera as ct
import numpy as np
import csv
from os import path, mkdir
from os import listdir as _listdir

from joblib import Parallel, delayed
from threadpoolctl import threadpool_limits

#---------------------------------------------------------------------------------------------#
# Importing DataMiner classes and functions
#---------------------------------------------------------------------------------------------#
from Common.DataDrivenConfig import Config_FGM
from Data_Generation.DataGenerator_Base import DataGenerator_Base
from Common.CommonMethods import ComputeLewisNumber
from Common.Properties import DefaultSettings_FGM,FGMVars

class DataGenerator_Cantera(DataGenerator_Base):
    """Generate flamelet data using Cantera.

    :param Config: Config_FGM class describing the flamelet generation settings.
    :type: Config_FGM
    """
    # Generate flamelet data from Cantera computation.
    _Config:Config_FGM

    # Save directory for computed flamelet data
    __matlab__output_dir:str = "./"

    __fuel_string:str = ''
    __oxidizer_string:str = ''

    __n_flamelets:int = DefaultSettings_FGM.Np_temp       # Number of adiabatic and burner flame computations per mixture fraction
    __n_mdot_extra_flamelets:int = 20                       # Number of interpolation steps for extra interpolated burner-stabilized flamelets
    __T_unburnt_upper:float = DefaultSettings_FGM.T_max   # Highest unburnt reactant temperature
    __T_unburnt_lower:float = DefaultSettings_FGM.T_min   # Lowest unburnt reactant temperature

    __reaction_mechanism:str = DefaultSettings_FGM.reaction_mechanism   # Cantera reaction mechanism
    __transport_model:str = DefaultSettings_FGM.transport_model

    __initial_grid_length:float = 0.2  # Flamelet grid width
    __loglevel:int = 0  # Cantera solver verbosity level (0=silent)
    __initial_grid_Np:int = 50          # Number of initial grid nodes.

    __define_equivalence_ratio:bool = not DefaultSettings_FGM.run_mixture_fraction # Define unburnt mixture via the equivalence ratio
    __unb_mixture_status:list[float] = []

    __translate_to_matlab:bool = False # Save a copy of the flamelet data file in Matlab table generator format

    __run_freeflames:bool = DefaultSettings_FGM.include_freeflames      # Run adiabatic flame computations
    __run_extra_interpolated_burnerflames:bool = True       # Run extra interpolated burner-stabilized flame computations
    __run_burnerflames:bool = DefaultSettings_FGM.include_burnerflames    # Run burner stabilized flame computations
    __run_equilibrium:bool = DefaultSettings_FGM.include_equilibrium    # Run chemical equilibrium computations
    __run_counterflames:bool = DefaultSettings_FGM.include_counterflames   # Run counter-flow diffusion flamelet simulations.
    __counterflow_fixed_strain:bool  = DefaultSettings_FGM.counterflow_fixed_strain  # Fixed-strain mode for counterflow flames.
    __counterflow_strain_rate:float  = DefaultSettings_FGM.counterflow_strain_rate   # Global strain rate [1/s] for fixed-strain mode.

    __u_fuel:float = 1.0       # Fuel stream velocity in counter-flow diffusion flame.
    __u_oxidizer:float = None   # Oxidizer stream velocity in counter-flow diffusion flame.

    __free_flame_refine:dict    = {"ratio": 3, "slope": 0.03, "curve": 0.03, "prune": 0.01}  # Refine criteria for free-flame solver.
    __burner_flame_refine:dict  = {"ratio": 3, "slope": 0.15, "curve": 0.15, "prune": 0.05}  # Refine criteria for burner-stabilized flame solver.
    __counter_flame_refine:dict = {"ratio": 3, "slope": 0.04, "curve": 0.06, "prune": 0.02}  # Refine criteria for counter-flow flame solver.
    __mdot_dH_target:float = 0.0  # Target enthalpy step between successive burner flames (0 = linspace, >0 = adaptive).
    __src_interp_exponent:float = 2.0  # Power-law exponent for source-term decay in extra interpolated burner flamelets.

    def __init__(self, Config:Config_FGM=None):
        """Constructur, load flamelet generation settings from Config_FGM.

        :param Config: Config_FGM containing respective settings.
        :type Config: Config_FGM
        """
        DataGenerator_Base.__init__(self, Config_in=Config)


        if Config is None:
            print("Initializing flamelet generator with default settings")
            self._Config = Config_FGM()
        else:
            print("Initializing flamelet generator from Config_FGM with name " + self._Config.GetConfigName())
            self.__SynchronizeSettings()

        return

    def __SynchronizeSettings(self):
        """Update settings from configuration
        """
        self.__fuel_string = self._Config.GetFuelString()
        self.__oxidizer_string = self._Config.GetOxidizerString()

        self.__reaction_mechanism = self._Config.GetReactionMechanism()
        self.__transport_model = self._Config.GetTransportModel()

        self.gas = ct.Solution(self._Config.GetReactionMechanism())

        self.__n_flamelets = self._Config.GetNpTemp()
        self.__n_mdot_flamelets = self._Config.GetNpMdot()
        self.__n_mdot_extra_flamelets = self._Config.GetNpMdotExtra()
        self.__mdot_dH_target = self._Config.GetMdotDHTarget()
        self.__src_interp_exponent = self._Config.GetSrcInterpExponent()
        self.__initial_grid_length = self._Config.GetInitialGridLength()
        [self.__T_unburnt_lower, self.__T_unburnt_upper] = self._Config.GetUnbTempBounds()

        self.__define_equivalence_ratio = (not self._Config.GetMixtureStatus())
        self.__unb_mixture_status = np.linspace(self._Config.GetMixtureBounds()[0], self._Config.GetMixtureBounds()[1], self._Config.GetNpMix())
        self.__run_freeflames = self._Config.GenerateFreeFlames()
        self.__run_extra_interpolated_burnerflames = self._Config.GenerateExtraInterpolatedBurnerFlames()
        self.__run_burnerflames = self._Config.GenerateBurnerFlames()
        self.__run_equilibrium = self._Config.GenerateEquilibrium()
        self.__run_counterflames = self._Config.GenerateCounterFlames()
        self.__counterflow_fixed_strain = getattr(self._Config, '_Config_FGM__counterflow_fixed_strain', DefaultSettings_FGM.counterflow_fixed_strain)
        self.__counterflow_strain_rate  = getattr(self._Config, '_Config_FGM__counterflow_strain_rate',  DefaultSettings_FGM.counterflow_strain_rate)
        self.__counterflow_fixed_strain = self._Config.GetCounterFlowFixedStrain()
        self.__counterflow_strain_rate  = self._Config.GetCounterFlowStrainRate()


        self.__PrepareOutputDirectories()
        self.__translate_to_matlab = self._Config.WriteMatlabFiles()
        if self.__translate_to_matlab:
            self.__PrepareOutputDirectories_Matlab()
        self._Config.ComputeMixFracConstants()
        self.z_i = self._Config.GetMixtureFractionCoefficients()
        self.c = self._Config.GetMixtureFractionConstant()
        return

    def SetLoglevel(self, loglevel:int=0):
        """Set Cantera solver verbosity level (0=silent)."""
        self.__loglevel = loglevel

    def SetInitialGridLength(self, length:float):
        """Set the initial domain length (in metres) for new flamelet grids."""
        self.__initial_grid_length = length

    def SetFreeFlameRefineCriteria(self, ratio:float=3, slope:float=0.03, curve:float=0.03, prune:float=0.01):
        """Set the grid refinement criteria for the free-flame (adiabatic) solver.

        :param ratio: Maximum ratio of adjacent grid spacings, defaults to 3.
        :param slope: Maximum relative slope of the solution, defaults to 0.03.
        :param curve: Maximum relative curvature of the solution, defaults to 0.03.
        :param prune: Threshold for grid point removal, defaults to 0.01.
        :raises Exception: if any value is not strictly positive.
        """
        if (any(v <= 0 for v in (ratio, slope, curve)) or (prune < 0)):
            raise Exception("All refine criteria values must be strictly positive.")
        self.__free_flame_refine = {"ratio": ratio, "slope": slope, "curve": curve, "prune": prune}
        return

    def SetBurnerFlameRefineCriteria(self, ratio:float=3, slope:float=0.15, curve:float=0.15, prune:float=0.01):
        """Set the grid refinement criteria for the burner-stabilized flame solver.

        :param ratio: Maximum ratio of adjacent grid spacings, defaults to 3.
        :param slope: Maximum relative slope of the solution, defaults to 0.15.
        :param curve: Maximum relative curvature of the solution, defaults to 0.15.
        :param prune: Threshold for grid point removal, defaults to 0.01.
        :raises Exception: if any value is not strictly positive.
        """
        if (any(v <= 0 for v in (ratio, slope, curve)) or (prune < 0)):
            raise Exception("All refine criteria values must be strictly positive.")
        self.__burner_flame_refine = {"ratio": ratio, "slope": slope, "curve": curve, "prune": prune}
        return

    def SetCounterFlameRefineCriteria(self, ratio:float=3, slope:float=0.04, curve:float=0.06, prune:float=0.01):
        """Set the grid refinement criteria for the counter-flow diffusion flame solver.

        :param ratio: Maximum ratio of adjacent grid spacings, defaults to 3.
        :param slope: Maximum relative slope of the solution, defaults to 0.04.
        :param curve: Maximum relative curvature of the solution, defaults to 0.06.
        :param prune: Threshold for grid point removal, defaults to 0.01.
        :raises Exception: if any value is not strictly positive.
        """
        if (any(v <= 0 for v in (ratio, slope, curve)) or (prune < 0)):
            raise Exception("All refine criteria values must be strictly positive.")
        self.__counter_flame_refine = {"ratio": ratio, "slope": slope, "curve": curve, "prune": prune}
        return

    def SetMdotDHTarget(self, dH_target:float):
        """Set the target enthalpy step between successive burner-stabilized flamelets.
        When greater than zero, mdot stepping is adaptive: after each solved flame the next
        Δmdot is rescaled so that |ΔH| matches this target.  When zero, the traditional
        uniform linspace spacing is used.

        :param dH_target: Target |ΔH| between successive burner flames in J/kg. 0 disables adaptive mode.
        :type dH_target: float
        :raises Exception: if the value is negative.
        """
        if dH_target < 0:
            raise Exception("mdot ΔH target must be non-negative (0 = linspace mode).")
        self.__mdot_dH_target = dH_target
        return

    def SetSrcInterpExponent(self, exponent:float=2.0):
        """Set power-law exponent for source-term decay in extra interpolated burner flamelets.
        1 = linear; > 1 = faster decay toward zero at the cold equilibrium end.

        :param exponent: Power-law exponent (>= 1), defaults to 2.
        :type exponent: float
        :raises Exception: If exponent is less than 1.
        """
        if exponent < 1.0:
            raise Exception("Source term interpolation exponent must be >= 1.")
        self.__src_interp_exponent = exponent
        return

    def SetFuelDefinition(self, fuel_species:list[str], fuel_weights:list[float]):
        """Manually define the fuel composition

        :param fuel_species: list of fuel species names.
        :type fuel_species: list[str]
        :param fuel_weights: list of fuel molar fraction weights.
        :type fuel_weights: list[float]
        :raises Exception: if no fuel species are provided.
        :raises Exception: if the number of species does not correspond to the number of weights.
        """
        self._Config.SetFuelDefinition(fuel_species, fuel_weights)
        self.__SynchronizeSettings()

        return

    def SetOxidizerDefinition(self, oxidizer_species:list[str], oxidizer_weights:list[float]):
        """Manually define the oxidizer composition

        :param oxidizer_species: list of oxidizer species names.
        :type oxidizer_species: list[str]
        :param oxidizer_weights: list of oxidizer molar fraction weights.
        :type oxidizer_weights: list[float]
        :raises Exception: if no oxidizer species are provided.
        :raises Exception: if the number of species does not correspond to the number of weights.
        """
        self._Config.SetOxidizerDefinition(oxidizer_species, oxidizer_weights)
        self.__SynchronizeSettings()

        return

    def SetNpTemp(self, n_flamelets_new:int):
        """Set the number of flamelets generated between the minimum and maximum reactant temperature manually.

        :param n_flamelets_new: number of flamelets generated between the minimum and maximum reactant temperature.
        :type n_flamelets_new: int
        :raises Exception: if the provided number is lower than one.
        """
        if n_flamelets_new < 1:
            raise Exception("Number of flamelets should be higher than one.")
        self.__n_flamelets = n_flamelets_new
        return

    def SetUnbTempBounds(self, T_unb_lower:float, T_unb_upper:float):
        """
        Define lower and upper reactant temperature for flamelet data generation.

        :param T_unb_lower: Lower reactant temperature in Kelvin.
        :type T_unb_lower: float
        :param T_unb_upper: Upper reactant temperature in Kelvin.
        :type T_unb_upper: float
        :raise: Exception: if lower temperature value exceeds upper temperature value.

        """

        if (T_unb_lower >= T_unb_upper):
            raise Exception("Lower unburnt temperature bound should be below upper bound.")
        else:
            self.__T_unburnt_upper = T_unb_upper
            self.__T_unburnt_lower = T_unb_lower
        return

    def RunMixtureFraction(self):
        """Define the mixture status as mixture fraction instead of equivalence ratio.
        """
        self.__define_equivalence_ratio = False
        return

    def RunEquivalenceRatio(self):
        """Define the mixture status as equivalence ratio instead of mixture fraction.
        """
        self.__define_equivalence_ratio = True
        return

    def RunFreeFlames(self, input:bool=True):
        """Include adiabatic free-flame data in the manifold.

        :param input: Generate adiabatic free-flame data.
        :type input: bool
        """
        self.__run_freeflames = input
        return

    def RunBurnerFlames(self, input:bool=True):
        """Include burner-stabilized flame data in the manifold.

        :param input: Generate burner-stabilized flamelet data.
        :type input: bool
        """
        self.__run_burnerflames = input
        return

    def RunEquilibrium(self, input:bool=True):
        """Include chemical equilibrium data in the manifold.

        :param input: Generate chemical equilibrium data.
        :type input: bool
        """
        self.__run_equilibrium = input
        return

    def RunCounterFlowFlames(self, input:bool=True):
        """Include counter-flow diffusion flame data in the manifold.

        :param input: Generate counter-flow diffusion flamelet data.
        :type input: bool
        """
        self.__run_counterflames = input
        return

    def RunExtraInterpolatedBurnerFlames(self, input:bool=True):
        """Include extra interpolated burner-stabilized flame data in the manifold.

        :param input: Generate extra interpolated burner-stabilized flamelet data.
        :type input: bool
        """
        self.__run_extra_interpolated_burnerflames = input
        return

    def SetNpMdotExtra(self, n_extra:int):
        """Set the number of interpolation steps for extra interpolated burner-stabilized flamelets.

        :param n_extra: Number of interpolated steps between the lowest-mdot burner flame and equilibrium.
        :type n_extra: int
        :raises Exception: if the provided number is lower than one.
        """
        if n_extra < 1:
            raise Exception("Number of extra interpolated mdot flamelets should be higher than one.")
        self.__n_mdot_extra_flamelets = n_extra
        return

    def SetMixtureValues(self, mixture_values:list[float]):
        """Set the reactant mixture status values manually.

        :param mixture_values: list of equivalence ratio or mixture fraction values.
        :type mixture_values: list[float]
        :raises Exception: If an empty list is provided.
        """
        if len(mixture_values) == 0:
            raise Exception("At least one mixture status value should be provided.")
        self.__unb_mixture_status = []
        for phi in mixture_values:
            self.__unb_mixture_status.append(phi)
        return

    def SetReactionMechanism(self, reaction_mechanism:str):
        """Define the reaction mechanism manually.

        :param reaction_mechanism: name of the reaction mechanism.
        :type reaction_mechanism: str
        """
        self._Config.SetReactionMechanism(reaction_mechanism)
        self.__SynchronizeSettings()
        return

    def SetTransportModel(self, transport_model:str):
        """Overwrite the transport mechanism from the loaded configuration.

        :param transport_model: Cantera transport model.
        :type transport_model: str
        """

        self._Config.SetTransportModel(transport_model)
        self.__SynchronizeSettings()
        return

    def TranslateToMatlab(self):
        """Save a copy of the flamelet data in Matlab TableMaster format.
        """
        self.__translate_to_matlab = True
        return

    def SetOutputDir(self, output_dir_new:str):
        """Define the flamelet data output directory manually.

        :param output_dir_new: Flamelet data output directory.
        :type output_dir_new: str
        :raises Exception: If provided directory doesn't exist.
        """
        self._Config.SetOutputDir(output_dir=output_dir_new)
        self.__SynchronizeSettings()
        return

    def SetMatlabOutputDir(self, output_dir_new):
        self.__matlab__output_dir = output_dir_new
        self.__PrepareOutputDirectories_Matlab()

    def __PrepareOutputDirectories(self):
        """Create sub-directories for the different types of flamelet data.
        """
        if (not path.isdir(self.GetOutputDir()+'/freeflame_data')) and self.__run_freeflames:
            mkdir(self.GetOutputDir()+'/freeflame_data')
        if (not path.isdir(self.GetOutputDir()+'/burnerflame_data')) and self.__run_burnerflames:
            mkdir(self.GetOutputDir()+'/burnerflame_data')
        if (not path.isdir(self.GetOutputDir()+'/equilibrium_data')) and self.__run_equilibrium:
            mkdir(self.GetOutputDir()+'/equilibrium_data')
        if (not path.isdir(self.GetOutputDir()+'/counterflame_data')) and self.__run_counterflames:
            mkdir(self.GetOutputDir()+'/counterflame_data')
        return

    def __PrepareOutputDirectories_Matlab(self):
        if (not path.isdir(self.__matlab__output_dir+'freeflame_data_MATLAB')) and self.__run_freeflames:
            mkdir(self.__matlab__output_dir+'freeflame_data_MATLAB')
        if (not path.isdir(self.__matlab__output_dir+'burnerflame_data_MATLAB')) and self.__run_burnerflames:
            mkdir(self.__matlab__output_dir+'burnerflame_data_MATLAB')
        if (not path.isdir(self.__matlab__output_dir+'equilibrium_data_MATLAB')) and self.__run_equilibrium:
            mkdir(self.__matlab__output_dir+'equilibrium_data_MATLAB')
        if (not path.isdir(self.__matlab__output_dir+'counterflame_data_MATLAB')) and self.__run_counterflames:
            mkdir(self.__matlab__output_dir+'counterflame_data_MATLAB')
        return


    def ComputeFreeFlames(self, mix_status:float, T_ub:float, i_freeflame:int=0, prev_flame:ct.FreeFlame=None):
        """Generate adiabatic free-flamelet data for a specific mixture fraction or equivalence ratio and reactant temperature.

        :param mix_status: Equivalence ratio or mixture fraction value.
        :type mix_status: float
        :param T_ub: Reactant temperature in Kelvin.
        :type T_ub: float
        :param i_freeflame: Solution index, defaults to 0
        :type i_freeflame: int, optional
        :param prev_flame: Converged FreeFlame from the previous temperature step to restart from, defaults to None.
        :type prev_flame: ct.FreeFlame, optional
        :return: Converged FreeFlame object on success, None on failure.
        :rtype: ct.FreeFlame or None
        """
        if self.__define_equivalence_ratio:
            folder_header = "phi"
        else:
            folder_header = "mixfrac"
        # Setting unburnt temperature and pressure
        self.gas.TP = T_ub, ct.one_atm
        # Defining mixture ratio based on equivalence ratio or mixture fraction.
        if self.__define_equivalence_ratio:
            self.gas.set_equivalence_ratio(mix_status, self.__fuel_string, self.__oxidizer_string)
        else:
            self.gas.set_mixture_fraction(mix_status, self.__fuel_string, self.__oxidizer_string)

        if prev_flame is None:
            # First flame: create a new object from scratch.
            initialgrid = np.linspace(0, self.__initial_grid_length, self.__initial_grid_Np)
            freeflame = ct.FreeFlame(self.gas, grid=initialgrid)
            freeflame.set_refine_criteria(**self.__free_flame_refine)
            freeflame.max_grid_points = 2000
            freeflame.transport_model = self.__transport_model
            freeflame.set_initial_guess(locs=[0.0, 0.3, 0.5, 1.0])
        else:
            # Subsequent flames: reuse the previous flame object directly.
            freeflame = prev_flame
            freeflame.inlet.T = T_ub
            freeflame.inlet.Y = self.gas.Y

        # Build output paths early so we can check for existing files.
        if not path.isdir(self.GetOutputDir()+'/freeflame_data/'):
            mkdir(self.GetOutputDir()+'/freeflame_data/')
        if not path.isdir(self.GetOutputDir()+'/freeflame_data/'+folder_header+'_'+str(round(mix_status, 6))):
            mkdir(self.GetOutputDir()+'/freeflame_data/'+folder_header+'_'+str(round(mix_status, 6)))
        freeflame_filename     = "freeflamelet_"+folder_header+str(round(mix_status,6))+"_Tu"+str(round(T_ub, 4))+".csv"
        filename_plus_folder   = self.GetOutputDir()+"/freeflame_data/"+folder_header+'_'+str(round(mix_status, 6)) + "/"+freeflame_filename
        yaml_plus_folder       = filename_plus_folder.replace('.csv', '.yaml')

        # If a YAML exists, restore as warm-start before solving.
        # Always re-solve and overwrite the CSV (use auto=False when restoring).
        warm_started = False
        if path.isfile(yaml_plus_folder):
            try:
                freeflame.restore(yaml_plus_folder, 'solution')
                warm_started = True
                print("Freeflame at %s %.3e, Tu %.3f: warm-start from YAML." % (folder_header, mix_status, T_ub))
            except Exception as exc:
                print("Freeflame at %s %.3e, Tu %.3f: YAML restore failed (%s), re-solving from scratch." % (folder_header, mix_status, T_ub, exc))

        # Try to solve the flamelet solution. If solution diverges, move on to next flamelet.
        try:
            freeflame.solve(loglevel=self.__loglevel, refine_grid=True, auto=(not warm_started and prev_flame is None))

            # Computing mass flow rate for later burner flame evaluation
            self.m_dot_free_flame = freeflame.velocity[0]*freeflame.density[0]

            # Check if mixture is burning
            if np.max(freeflame.T) <= DefaultSettings_FGM.T_threshold:
                print("Flamelet at %s %.3e, Tu %.3f is not burning" % (folder_header, mix_status, T_ub))
                return None

            variables, data_calc = self.__SaveFlameletData(freeflame, self.gas)

            # Directories were already created before the skip-check above.
            freeflame.save(yaml_plus_folder, 'solution', overwrite=True)
            fid = open(filename_plus_folder, 'w+')
            fid.write(variables + "\n")
            csvWriter = csv.writer(fid)
            csvWriter.writerows(data_calc)
            fid.close()

            if self.__translate_to_matlab:
                if not path.isdir(self.__matlab__output_dir+'/freeflame_data_MATLAB/'+folder_header+'_'+str(round(mix_status, 6))):
                        mkdir(self.__matlab__output_dir+'/freeflame_data_MATLAB/'+folder_header+'_'+str(round(mix_status, 6)))
                self.__TranslateToMatlabFile(filename_plus_folder,freeflame_filename, self.__matlab__output_dir + "/freeflame_data_MATLAB/"+folder_header+'_'+str(round(mix_status, 6)) + "/")
            self.last_Y_flamelet = freeflame.Y
            self.last_h_flamelet = freeflame.enthalpy_mass
            self.last_T_flamelet = freeflame.T

            print("Successfull Freeflame simulation at "+folder_header+": "+str(mix_status)+ " T_u: " +str(T_ub) + " Np=%d" % len(freeflame.grid) + " ("+str(i_freeflame+1)+"/"+str(self.__n_flamelets)+")")
            return freeflame

        except:
            print("Unsuccessfull Freeflame simulation at "+folder_header+": "+str(mix_status)+ " T_u: " +str(T_ub) + " ("+str(i_freeflame+1)+"/"+str(self.__n_flamelets)+")")
            return None

    def compute_SingleBurnerFlame(self, mix_status:float, T_burner:float, m_dot:float, prev_flame:ct.BurnerFlame=None):
        """Compute the solution of a single burner-stabilized flamelet.

        :param mix_status: mixture fraction or equivalence ratio.
        :type mix_status: float
        :param T_burner: burner plate temperature
        :type T_burner: float
        :param m_dot: mass flux [kg m/s]
        :type m_dot: float
        :param prev_flame: Converged BurnerFlame object from the previous mass-flux step,
            reused directly to avoid file I/O and object reconstruction. When None the flame
            is solved from scratch using Cantera's auto multi-stage strategy, defaults to None.
        :type prev_flame: ct.BurnerFlame, optional
        :return: converged burner flame object
        :rtype: cantera.BurnerFlame
        """
        if prev_flame is None:
            # First flame: set gas state for the fresh BurnerFlame object.
            # ChemEquil (called internally by auto=True) requires T >= 300 K,
            # so clamp to 300 K here; burner.T is restored to T_burner afterwards.
            T_init = max(T_burner, 300.0)
            self.gas.TP = T_init, DefaultSettings_FGM.pressure
            if self.__define_equivalence_ratio:
                self.gas.set_equivalence_ratio(mix_status, self.__fuel_string, self.__oxidizer_string)
            else:
                self.gas.set_mixture_fraction(mix_status, self.__fuel_string, self.__oxidizer_string)
            initialgrid = np.linspace(0, self.__initial_grid_length, self.__initial_grid_Np)
            burner_flame = ct.BurnerFlame(self.gas, grid=initialgrid)
            burner_flame.burner.T = T_burner
            burner_flame.set_refine_criteria(**self.__burner_flame_refine)
            burner_flame.transport_model = self.__transport_model
        else:
            # Subsequent burner flames: reuse the previous flame object directly.
            # Just update the mass flux — no file I/O, no object reconstruction,
            # and the solver's internal state (Jacobian estimate) carries over.
            burner_flame = prev_flame
        burner_flame.burner.mdot = m_dot

        if self.__define_equivalence_ratio:
            mix_label = "phi"
        else:
            mix_label = "mixfrac"
        print("  Burner flame: %s=%.4f  T_burner=%.1f K  mdot=%.4e kg/m2/s" % (mix_label, mix_status, T_burner, m_dot))

        # For the first flame use auto=True so Cantera's built-in multi-stage strategy
        # (frozen chemistry → reacting → refine) reliably finds the ignited solution.
        # For subsequent flames the previous converged solution is already a great
        # initial guess, so auto=False is sufficient and avoids redundant work.
        use_auto = (prev_flame is None)
        burner_flame.solve(loglevel=self.__loglevel, refine_grid=True, auto=use_auto)

        return burner_flame

    def ComputeBurnerFlames(self, mix_status:float, m_dot:np.ndarray=None, T_burner:float=None, free_flame:ct.FreeFlame=None):
        """Generate burner-stabilized flamelet data for a specific mixture fraction or equivalence ratio and mass flux.

        :param mix_status: Equivalence ratio or mixture fraction value.
        :type mix_status: float
        :param m_dot: Pre-built mass flux array (kg s^{-1} m^{-1}). When None and
            ``__mdot_dH_target > 0`` adaptive ΔH stepping is used; when None and
            ``__mdot_dH_target == 0`` a uniform linspace array is generated automatically.
        :type m_dot: np.ndarray, optional
        :param free_flame: Converged FreeFlame used to seed the first burner flame and set
            its domain length. When None the first burner flame is solved from scratch.
        :type free_flame: ct.FreeFlame, optional
        """
        if self.__define_equivalence_ratio:
            folder_header = "phi"
        else:
            folder_header = "mixfrac"

        if T_burner == None:
            T_burner = self.__T_unburnt_lower
        print("Computing burner flamelets at "+folder_header+": "+str(mix_status)+ " T_burner: " +str(T_burner))
        self.gas.TP = T_burner, ct.one_atm

        if self.__define_equivalence_ratio:
            self.gas.set_equivalence_ratio(mix_status, self.__fuel_string, self.__oxidizer_string)
        else:
            self.gas.set_mixture_fraction(mix_status, self.__fuel_string, self.__oxidizer_string)

        m_dot_start = 0.98 * self.m_dot_free_flame
        m_dot_stop  = 0.001 * self.m_dot_free_flame

        if m_dot is None:
            if self.__mdot_dH_target > 0:
                # Adaptive ΔH stepping: step mdot so that successive flames are
                # separated by approximately __mdot_dH_target J/kg in enthalpy.
                m_dot_iter = None  # sentinel: adaptive path
            else:
                m_dot_iter = np.linspace(m_dot_start, m_dot_stop, self.__n_mdot_flamelets + 1)[:-1]
        else:
            m_dot_iter = m_dot

        prev_burner_flame = None
        H_prev = None
        delta_mdot = (m_dot_start - m_dot_stop) / self.__n_mdot_flamelets  # initial seed step
        m_dot_current = m_dot_start
        i_burnerflame = 0

        def _iterate_mdot():
            """Yield (index, m_dot_value) for either the pre-built or adaptive path."""
            nonlocal m_dot_current, delta_mdot, H_prev, prev_burner_flame
            if m_dot_iter is not None:
                for idx, val in enumerate(m_dot_iter):
                    yield idx, val
            else:
                idx = 0
                while m_dot_current >= m_dot_stop:
                    yield idx, m_dot_current
                    idx += 1

        for i_burnerflame, m_dot_next in _iterate_mdot():
            # Build output paths before solving so we can check for existing files.
            if not path.isdir(self.GetOutputDir()+'/burnerflame_data/'):
                mkdir(self.GetOutputDir()+'/burnerflame_data/')
            if not path.isdir(self.GetOutputDir()+'/burnerflame_data/'+folder_header+'_'+str(round(mix_status, 6))):
                mkdir(self.GetOutputDir()+'/burnerflame_data/'+folder_header+'_'+str(round(mix_status, 6)))
            burnerflame_filename = "burnerflamelet_%s%.6f_mdot%.5f.csv" % (folder_header, mix_status, m_dot_next)
            filename_plus_folder = self.GetOutputDir()+"/burnerflame_data/"+folder_header+'_'+str(round(mix_status, 6)) + "/"+burnerflame_filename
            yaml_plus_folder     = filename_plus_folder.replace('.csv', '.yaml')

            # If a YAML exists, restore as warm-start before solving.
            # Always re-solve and overwrite the CSV (use auto=False when restoring).
            warm_started_burner = False
            if path.isfile(yaml_plus_folder):
                try:
                    # Reconstruct a minimal BurnerFlame shell to restore into.
                    self.gas.TP = T_burner, ct.one_atm
                    if self.__define_equivalence_ratio:
                        self.gas.set_equivalence_ratio(mix_status, self.__fuel_string, self.__oxidizer_string)
                    else:
                        self.gas.set_mixture_fraction(mix_status, self.__fuel_string, self.__oxidizer_string)
                    initialgrid = np.linspace(0, self.__initial_grid_length, self.__initial_grid_Np)
                    burner_flame = ct.BurnerFlame(self.gas, grid=initialgrid)
                    burner_flame.burner.T    = T_burner
                    burner_flame.burner.mdot = m_dot_next
                    burner_flame.set_refine_criteria(**self.__burner_flame_refine)
                    burner_flame.transport_model = self.__transport_model
                    burner_flame.restore(yaml_plus_folder, 'solution')
                    prev_burner_flame = burner_flame
                    warm_started_burner = True
                    print("Burnerflame at %s %.3e, mdot %.2e: warm-start from YAML." % (folder_header, mix_status, m_dot_next))
                except Exception as exc:
                    print("Burnerflame at %s %.3e, mdot %.2e: YAML restore failed (%s), re-solving." % (folder_header, mix_status, m_dot_next, exc))

            try:
                burner_flame = self.compute_SingleBurnerFlame(mix_status, T_burner, m_dot_next,
                                                              prev_burner_flame if warm_started_burner else prev_burner_flame)
                if np.max(burner_flame.T) <= DefaultSettings_FGM.T_threshold:
                    print("Burnerflame at %s %.3e, mdot %.2e is not burning" % (folder_header, mix_status, m_dot_next))
                    if m_dot_iter is None:
                        m_dot_current -= delta_mdot
                    continue

                # Extracting flamelet data
                variables, data_calc = self.__SaveFlameletData(burner_flame, self.gas)

                # Directories and filenames were already built before the skip-check above.
                burner_flame.save(yaml_plus_folder, 'solution', overwrite=True)
                fid = open(filename_plus_folder, 'w+')
                fid.write(variables + "\n")
                csvWriter = csv.writer(fid)
                csvWriter.writerows(data_calc)
                fid.close()

                if self.__translate_to_matlab:
                    if not path.isdir(self.__matlab__output_dir+'/burnerflame_data_MATLAB/'+folder_header+'_'+str(round(mix_status, 6))):
                        mkdir(self.__matlab__output_dir+'/burnerflame_data_MATLAB/'+folder_header+'_'+str(round(mix_status, 6)))
                    self.__TranslateToMatlabFile(filename_plus_folder,burnerflame_filename, self.__matlab__output_dir + "/burnerflame_data_MATLAB/"+folder_header+'_'+str(round(mix_status, 6)) + "/")

                Y_max, Y_min = np.max(burner_flame.Y,axis=1), np.min(burner_flame.Y,axis=1)
                delta_Y_flamelet = Y_max - Y_min
                if max(delta_Y_flamelet) > 1e-5:
                    self.last_Y_flamelet = burner_flame.Y
                    self.last_h_flamelet = burner_flame.enthalpy_mass
                    self.last_T_flamelet = burner_flame.T

                # Adaptive ΔH step control: rescale delta_mdot to target __mdot_dH_target.
                if m_dot_iter is None:
                    H_current = float(burner_flame.enthalpy_mass[0])
                    if H_prev is not None:
                        dH_actual = abs(H_prev - H_current)
                        if dH_actual > 0:
                            scale = float(np.clip(self.__mdot_dH_target / dH_actual, 0.2, 5.0))
                            delta_mdot *= scale
                    H_prev = H_current
                    m_dot_current -= delta_mdot

                prev_burner_flame = burner_flame
                print("Successfull burnerflame simulation at "+folder_header+": "+ str(mix_status)+" mdot: " + str(m_dot_next)+ " ("+str(i_burnerflame+1)+"/"+str(self.__n_mdot_flamelets)+")")

            except:
                print("Unsuccessfull burnerflame simulation at "+folder_header+": "+ str(mix_status)+" mdot: " + str(m_dot_next)+ " ("+str(i_burnerflame+1)+"/"+str(self.__n_mdot_flamelets)+")")
                if m_dot_iter is None:
                    m_dot_current -= delta_mdot

    def ComputeCounterFlowFlames(self, v_fuel:float, v_ox:float, T_ub:float):
        """Generate counter-flow diffusion flamelet data for a given temperature, and reactant velocities.
        Strain rate is gradually increased until extinction in order to distribute data over the progress variable spectrum.

        :param v_fuel: Fuel reactant velocity in meters per second.
        :type v_fuel: float
        :param v_ox: Oxidizer reactant velocity in meters per second.
        :type v_ox: float
        :param T_ub: Reactant temperature in Kelvin.
        :type T_ub: float
        :raises Exception: If either of the velocity values is lower than zero.
        :raises Exception: If the reactant temperature is lower than 200 K.
        """
        if (v_fuel <= 0) or (v_ox <= 0):
            raise Exception("Reactant velocities should be higher than zero.")
        if T_ub < 200:
            raise Exception("Reactant temperature should be higher than 200K.")
        counterflame = ct.CounterflowDiffusionFlame(self.gas, width=self.__initial_grid_length)

        self.gas.set_mixture_fraction(1.0, self.__fuel_string, self.__oxidizer_string)
        self.gas.TP = T_ub, ct.one_atm
        rho_fuel = self.gas.density

        self.gas.set_mixture_fraction(0.0, self.__fuel_string, self.__oxidizer_string)
        self.gas.TP = T_ub, ct.one_atm
        rho_oxidizer = self.gas.density

        counterflame.P = ct.one_atm
        counterflame.fuel_inlet.Y = self.__fuel_string
        counterflame.fuel_inlet.T = T_ub
        counterflame.fuel_inlet.mdot = rho_fuel*v_fuel
        counterflame.oxidizer_inlet.Y = self.__oxidizer_string
        counterflame.oxidizer_inlet.T = T_ub
        counterflame.oxidizer_inlet.mdot = rho_oxidizer*v_ox
        counterflame.set_refine_criteria(**self.__counter_flame_refine)

        counterflame.solve(loglevel=self.__loglevel, auto=True)
        variables, data_calc = self.__SaveFlameletData(counterflame, self.gas)

        counterflame_filename = "counterflamelet_strain_0_Tu"+str(round(T_ub, 4))+".csv"
        if not path.isdir(self.GetOutputDir()+"/counterflame_data"):
            mkdir(self.GetOutputDir()+"/counterflame_data")
        fid = open(self.GetOutputDir()+"/counterflame_data/"+counterflame_filename, 'w+')
        fid.write(variables + "\n")
        csvWriter = csv.writer(fid)
        csvWriter.writerows(data_calc)
        fid.close()
        # Compute counterflow diffusion flames at increasing strain rates at 1 bar
        # The strain rate is assumed to increase by 25% in each step until the flame is
        # extinguished
        strain_factor = 1.25
        # Exponents for the initial solution variation with changes in strain rate
        # Taken from Fiala and Sattelmayer (2014)
        exp_d_a = -0.05
        exp_u_a = 1. / 2.
        exp_V_a = 1.
        exp_lam_a = 2.
        exp_mdot_a = 1. / 2.

        n_iter = 1
        strain_overload = False
        while not strain_overload:
            # Create an initial guess based on the previous solution
            # Update grid
            counterflame.flame.grid *= strain_factor ** exp_d_a
            normalized_grid = counterflame.grid / (counterflame.grid[-1] - counterflame.grid[0])
            # Update mass fluxes
            counterflame.fuel_inlet.mdot *= strain_factor ** exp_mdot_a
            counterflame.oxidizer_inlet.mdot *= strain_factor ** exp_mdot_a
            # Update velocities
            counterflame.set_profile('velocity', normalized_grid,
                        counterflame.velocity * strain_factor ** exp_u_a)
            counterflame.set_profile('spread_rate', normalized_grid,
                        counterflame.spread_rate * strain_factor ** exp_V_a)
            # Update pressure curvature
            counterflame.set_profile('lambda', normalized_grid, counterflame.L * strain_factor ** exp_lam_a)

            try:
                # Try solving the flame
                counterflame.solve(loglevel=self.__loglevel)
                self.last_counterflame_massfraction = counterflame.Y
                variables, data_calc = self.__SaveFlameletData(counterflame, self.gas)

                counterflame_filename = "counterflamelet_strain_"+str(n_iter)+"_Tu"+str(round(T_ub, 4))+".csv"
                fid = open(self.GetOutputDir()+"/counterflame_data/"+counterflame_filename, 'w+')
                fid.write(variables + "\n")
                csvWriter = csv.writer(fid)
                csvWriter.writerows(data_calc)
                fid.close()
                print("Successful Counter-Flow Diffusion Flame at Strain Iteration " + str(n_iter))
            except:
                print("Unsuccessful Counter-Flow Diffusion Flame at Strain Iteration " + str(n_iter))
                strain_overload = True
            n_iter += 1

    def ComputeCounterFlowFlamesFixedStrain(self, T_ub:float):
        """Solve one counterflow diffusion flame at a fixed global strain rate.

        Inlet velocities are recomputed from the momentum-balance condition at
        every temperature level so the global strain rate stays exactly at the
        configured value regardless of density changes.

        A YAML restart file is saved alongside each CSV so that subsequent runs
        can warm-start from the nearest already-solved flame.  If a CSV already
        exists for this temperature it is skipped without re-solving.

        Output is written to counterflame_data/counterflamelet_fixedstrain_Tu{T}.csv.

        :param T_ub: Reactant temperature in Kelvin.
        :type T_ub: float
        :raises Exception: If the reactant temperature is lower than 200 K.
        """
        if T_ub < 200:
            raise Exception("Reactant temperature should be higher than 200 K.")

        out_dir  = self.GetOutputDir() + "/counterflame_data"
        if not path.isdir(out_dir):
            mkdir(out_dir)

        t_tag    = str(round(T_ub, 4))
        csv_path = path.join(out_dir, "counterflamelet_fixedstrain_Tu" + t_tag + ".csv")
        yaml_path = path.join(out_dir, "counterflamelet_fixedstrain_Tu" + t_tag + ".yaml")

        # If a YAML for this exact temperature exists, use it as the restart and
        # re-solve (overwriting the old CSV).  If no YAML exists, warm-start from
        # the nearest available YAML.  If no YAML exists at all, use auto-solve.
        use_auto = True

        a = self.__counterflow_strain_rate
        L = self.__initial_grid_length

        # Momentum-balanced velocities: a_global = (u_fuel + u_ox) / L,
        #   rho_fuel * u_fuel = rho_ox * u_ox
        self.gas.set_mixture_fraction(1.0, self.__fuel_string, self.__oxidizer_string)
        self.gas.TP = T_ub, ct.one_atm
        rho_fuel = self.gas.density_mass
        Y_fuel   = self.gas.Y.copy()

        self.gas.set_mixture_fraction(0.0, self.__fuel_string, self.__oxidizer_string)
        self.gas.TP = T_ub, ct.one_atm
        rho_ox = self.gas.density_mass
        Y_ox   = self.gas.Y.copy()

        u_fuel = a * L / (1.0 + rho_fuel / rho_ox)
        u_ox   = u_fuel * rho_fuel / rho_ox

        flame = ct.CounterflowDiffusionFlame(self.gas, width=L)
        flame.fuel_inlet.mdot     = rho_fuel * u_fuel
        flame.fuel_inlet.Y        = Y_fuel
        flame.fuel_inlet.T        = T_ub
        flame.oxidizer_inlet.mdot = rho_ox * u_ox
        flame.oxidizer_inlet.Y    = Y_ox
        flame.oxidizer_inlet.T    = T_ub
        flame.P = ct.one_atm
        flame.set_refine_criteria(**self.__counter_flame_refine)

        # Warm-start from the nearest already-solved YAML if available.
        existing_yamls = [f for f in _listdir(out_dir)
                          if f.startswith("counterflamelet_fixedstrain_Tu") and f.endswith(".yaml")]
        if existing_yamls:
            def _parse_T(fname):
                try:
                    return float(fname.replace("counterflamelet_fixedstrain_Tu", "").replace(".yaml", ""))
                except ValueError:
                    return None
            candidates = [(abs(_parse_T(f) - T_ub), f) for f in existing_yamls if _parse_T(f) is not None]
            if candidates:
                _, nearest_yaml = min(candidates)
                guess_path = path.join(out_dir, nearest_yaml)
                try:
                    flame.restore(guess_path, 'solution')
                    # Reapply boundary conditions after restore.
                    flame.fuel_inlet.mdot     = rho_fuel * u_fuel
                    flame.fuel_inlet.Y        = Y_fuel
                    flame.fuel_inlet.T        = T_ub
                    flame.oxidizer_inlet.mdot = rho_ox * u_ox
                    flame.oxidizer_inlet.Y    = Y_ox
                    flame.oxidizer_inlet.T    = T_ub
                    use_auto = False
                    print("Counter-flow flame at Tu=%.1f K: warm-start from %s" % (T_ub, nearest_yaml))
                except Exception as exc:
                    print("Counter-flow flame at Tu=%.1f K: warm-start failed (%s), using auto-solve." % (T_ub, exc))

        flame.solve(loglevel=self.__loglevel, auto=use_auto)
        print("Counter-flow flame (fixed strain %.1f 1/s) at Tu=%.1f K: max T = %.1f K, %d grid points"
              % (a, T_ub, np.max(flame.T), len(flame.grid)))

        # Save YAML for future warm-starts before writing the CSV.
        flame.save(yaml_path, 'solution', overwrite=True)

        variables, data_calc = self.__SaveFlameletData(flame, self.gas)

        fid = open(csv_path, 'w+')
        fid.write(variables + "\n")
        csvWriter = csv.writer(fid)
        csvWriter.writerows(data_calc)
        fid.close()
        return

    def ComputeEquilibrium(self, mix_status:float, T_range:np.ndarray[float], burnt:bool=False):
        """Generate chemical equilibrium data for a given mixture status and temperature range.

        :param mix_status: Mixture fraction or equivalence ratio.
        :type mix_status: float
        :param T_range: Reactant or product temperature range.
        :type T_range: np.array[float]
        :param burnt: Compute reaction product properties, defaults to False
        :type burnt: bool, optional
        """
        if self.__define_equivalence_ratio:
            folder_header = "phi"
        else:
            folder_header = "mixfrac"

        gas_eq = ct.Solution(self.__reaction_mechanism)

        if burnt:
            fileHeader = "equilibrium_b_"
        else:
            fileHeader = "equilibrium_ub_"
        if not path.isdir(self.GetOutputDir()+'/equilibrium_data/'):
                        mkdir(self.GetOutputDir()+'/equilibrium_data/')
        if not path.isdir(self.GetOutputDir() + "/equilibrium_data/" + folder_header+"_"+str(round(mix_status,6))):
            mkdir(self.GetOutputDir() + "/equilibrium_data/" + folder_header+"_"+str(round(mix_status,6)))

        is_lean = False
        if self.__define_equivalence_ratio:
            gas_eq.set_equivalence_ratio(mix_status, self.__fuel_string, self.__oxidizer_string)
            if mix_status <= 1.0:
                is_lean = True
        else:
            gas_eq.set_equivalence_ratio(1.0, self.__fuel_string, self.__oxidizer_string)
            z_stoch = gas_eq.mixture_fraction(self.__fuel_string, self.__oxidizer_string)
            if mix_status <= z_stoch:
                is_lean = True
            gas_eq.set_mixture_fraction(mix_status, self.__fuel_string, self.__oxidizer_string)

        gas_eq.TP = max(T_range), ct.one_atm
        H_max = gas_eq.enthalpy_mass
        # In case of reaction products, set the maximum enthalpy to that of the reactants at the maximum temperature.
        if burnt:
            gas_eq.TP = min(T_range), ct.one_atm
            if is_lean:
                gas_eq.equilibrate("TP")
            else:
                gas_eq.equilibrate('HP')
            gas_eq.HP = H_max, ct.one_atm
            T_range = np.linspace(min(T_range), gas_eq.T, len(T_range))

        for i, T in enumerate(T_range):

            gas_eq.TP = T, ct.one_atm

            if i == 0:
                if not path.isdir(self.GetOutputDir()+'/equilibrium_data/'+folder_header+'_'+str(round(mix_status, 6))):
                    mkdir(self.GetOutputDir()+'/equilibrium_data/'+folder_header+'_'+str(round(mix_status, 6)))
                variables, data_calc = self.__SaveFlameletData(gas_eq, self.gas)
                fid = open(self.GetOutputDir()+"/equilibrium_data/"+folder_header+"_"+str(round(mix_status,6))+"/"+ fileHeader +folder_header+"_"+str(round(mix_status,6))+".csv", 'w+')
                fid.write(variables + "\n")
                fid.close()
            else:
                variables, data_calc_2 = self.__SaveFlameletData(gas_eq, self.gas)
                data_calc = np.append(data_calc, data_calc_2, axis=0)

        eq_filename = fileHeader +folder_header+"_"+str(round(mix_status,6))+".csv"
        filename_plus_folder = self.GetOutputDir()+"/equilibrium_data/"+folder_header+"_"+str(round(mix_status,6))+"/"+ eq_filename
        fid = open(filename_plus_folder, 'a+')
        csvWriter = csv.writer(fid)
        csvWriter.writerows(data_calc)
        fid.close()

        if self.__translate_to_matlab:
            if not path.isdir(self.__matlab__output_dir+'/equilibrium_data_MATLAB/'+folder_header+'_'+str(round(mix_status, 6))):
                    mkdir(self.__matlab__output_dir+'/equilibrium_data_MATLAB/'+folder_header+'_'+str(round(mix_status, 6)))
            self.__TranslateToMatlabFile(filename_plus_folder, eq_filename, self.__matlab__output_dir + "/equilibrium_data_MATLAB/"+folder_header+'_'+str(round(mix_status, 6)) + "/")


    def ComputeFlameletsOnMixStatus(self, mix_status:float):
        """Generate flamelet data for a given mixture fraction or equivalence ratio.

        :param mix_status: Mixture fraction or equivalence ratio value.
        :type mix_status: float
        :raises Exception: If mixture status value is below zero.
        """

        if mix_status < 0:
            raise Exception("Mixture status value should be positive.")

        T_unburnt_range = np.linspace(self.__T_unburnt_upper, self.__T_unburnt_lower, self.__n_flamelets)
        # Generate adiabatic freeflame data
        if self.__run_freeflames:
            print("Starting Free Flame simulations ...")
            # Generate and safe adiabatic flamelet data.
            prev_free_flame = None
            for i_freeflame, T_ub in enumerate(T_unburnt_range):
                prev_free_flame = self.ComputeFreeFlames(mix_status=mix_status, T_ub=T_ub, i_freeflame=i_freeflame, prev_flame=prev_free_flame)

        # Generate burner-stabilized flamelet data
        if self.__run_burnerflames:
            print("Starting Burner Stabilized simulations ...")
            # Generate a single freeflamelet solution for reference
            if not self.__run_freeflames:
                prev_free_flame = self.ComputeFreeFlames(mix_status=mix_status, T_ub=self.__T_unburnt_lower, i_freeflame=0)

            # Generate burner-stabilized flamelet data.
            # When __mdot_dH_target > 0 pass m_dot=None to enable adaptive ΔH stepping;
            # otherwise ComputeBurnerFlames generates a uniform linspace internally.
            self.ComputeBurnerFlames(mix_status=mix_status, m_dot=None, free_flame=prev_free_flame)

        # Generate chemical equilibrium data
        if self.__run_equilibrium:

            # Generate unburnt reactants data.
            self.ComputeEquilibrium(mix_status=mix_status,\
                                    T_range=np.linspace(self.__T_unburnt_lower, self.__T_unburnt_upper, 2*self.__n_flamelets),\
                                    burnt=False)

            # Generate reaction products data.
            self.ComputeEquilibrium(mix_status=mix_status,\
                                    T_range=np.linspace(self.__T_unburnt_lower, self.__T_unburnt_upper, 2*self.__n_flamelets),\
                                    burnt=True)

        # Generate extra interpolated burner-stabilized flamelet data.
        if self.__run_extra_interpolated_burnerflames:
            if self.__define_equivalence_ratio:
                folder_header = "phi"
            else:
                folder_header = "mixfrac"

            burnerfolder = self.GetOutputDir() + '/burnerflame_data/' + folder_header + '_' + str(round(mix_status, 6))
            if not path.isdir(burnerfolder):
                return

            # Find the burner flame file with the lowest mass flux.
            burner_files = [f for f in _listdir(burnerfolder)
                            if f.startswith("burnerflamelet_") and f.endswith(".csv")
                            and "_int" not in f]
            if not burner_files:
                return

            def _mdot_from_file(fname):
                try:
                    fpath = burnerfolder + "/" + fname
                    with open(fpath, newline='') as _f:
                        _reader = csv.reader(_f)
                        _header = next(_reader)
                        _row = list(map(float, next(_reader)))
                    _idx_vel = _header.index("Velocity")
                    _idx_rho = _header.index(FGMVars.Density.name)
                    return _row[_idx_vel] * _row[_idx_rho]
                except Exception:
                    return float("inf")

            last_file = min(burner_files, key=_mdot_from_file)
            last_filepath = burnerfolder + "/" + last_file
            print("Interpolating from lowest-mdot burner flame: " + last_file)

            # Read header and all rows from the lowest-mdot burner flame.
            with open(last_filepath, newline='') as f:
                reader = csv.reader(f)
                headerline = next(reader)
                burner_rows = [list(map(float, row)) for row in reader]

            # Read first row of cooled burnt equilibrium as the endpoint.
            eq_filename = "equilibrium_b_" + folder_header + "_" + str(round(mix_status, 6)) + ".csv"
            eq_filepath = path.join(self.GetOutputDir(), "equilibrium_data",
                                    folder_header + "_" + str(round(mix_status, 6)), eq_filename)
            if not path.isfile(eq_filepath):
                print("Equilibrium file not found, skipping interpolated flames: " + eq_filepath)
                return
            eq_data = np.loadtxt(eq_filepath, delimiter=',', skiprows=1, ndmin=2)
            eq_row = eq_data[0].tolist()

            # Identify source-term column indices: Y_dot_net/pos/neg and heat release rate
            # decay faster than linear because Arrhenius kinetics make rates negligible
            # well before the cold equilibrium end is reached.
            _src_prefixes = ("Y_dot_net-", "Y_dot_pos-", "Y_dot_neg-")
            _src_cols = frozenset(
                j for j, h in enumerate(headerline)
                if h.startswith(_src_prefixes) or h == FGMVars.Heat_Release.name
            )

            # Interpolate N_int synthetic flames from the lowest-mdot burner flame to the
            # equilibrium endpoint.  State variables (T, Y, h, …) use linear weights;
            # source-term columns use a power-law weight so they reach zero sooner.
            N_int = self.__n_mdot_extra_flamelets
            exp = self.__src_interp_exponent
            for i in range(N_int):
                ratio = float(i + 1) / float(N_int)
                w_a_lin = 1.0 - ratio                   # linear weight toward burner flame
                w_a_src = (1.0 - ratio) ** exp          # power-law weight for source terms
                w_b_src = 1.0 - w_a_src
                int_filename = burnerfolder + "/" + "burnerflamelet_%s%.6f_int%04d.csv" % (folder_header, mix_status, i)
                with open(int_filename, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(headerline)
                    for burner_row in burner_rows:
                        row_mix = [
                            w_a_src * a + w_b_src * b if j in _src_cols
                            else w_a_lin * a + ratio * b
                            for j, (a, b) in enumerate(zip(burner_row, eq_row))
                        ]
                        writer.writerow(row_mix)

        return

    def ComputeFlamelets(self):
        """Generate and store all flamelet data for the current settings.
        """

        T_unburnt_range = np.linspace(self.__T_unburnt_upper, self.__T_unburnt_lower, self.__n_flamelets)

        # Generate counter-flow diffusion flamelet data
        if self.__run_counterflames:

            if not path.isdir(self.GetOutputDir()+'counterflame_data'):
                mkdir(self.GetOutputDir()+'counterflame_data')
            for T_ub in T_unburnt_range:
                if self.__counterflow_fixed_strain:
                    self.ComputeCounterFlowFlamesFixedStrain(T_ub=T_ub)
                else:
                    self.gas.TP = T_ub, 101325
                    self.gas.set_mixture_fraction(1.0, self.__fuel_string, self.__oxidizer_string)
                    rho_fuel = self.gas.density_mass
                    rhou_fuel = rho_fuel * self.__u_fuel
                    self.gas.set_mixture_fraction(0.0, self.__fuel_string, self.__oxidizer_string)
                    rho_ox = self.gas.density_mass
                    self.__u_oxidizer = rhou_fuel / rho_ox
                    self.ComputeCounterFlowFlames(v_fuel=self.__u_fuel, v_ox=self.__u_oxidizer, T_ub=T_ub)

        # Generate all other flamelet types.
        for mix_status in self.__unb_mixture_status:
            self.ComputeFlameletsOnMixStatus(mix_status)

    def __SaveFlameletData(self,flame, gas:ct.Solution):
        """Save flamelet or chemical equilibrium data in csv file.

        :param flame: Converged Cantera flamelet class.
        :type flame: cantera.FreeFlame, cantera.BurnerFlame, or cantera.CounterFlowDiffusionFlame
        :param gas: Cantera Solution object containing molecular properties of the respective mixture.
        :type gas: cantera.Solution
        :return: Flamelet variables string and data array
        :rtype: str, np.ndarray
        """

        # Check if chemical equilibrium or flamelet data are supplied.
        flame_is_gas = (np.shape(flame.Y) == np.shape(gas.Y))
        molar_weights = np.reshape(gas.molecular_weights, [gas.n_species, 1])

        # Extract species mass and molar fractions, reaction rates, and species specific heat values.
        if flame_is_gas:
            Y = np.reshape(flame.Y, [gas.n_species, 1])
            X = np.reshape(flame.X, [gas.n_species, 1])
            net_reaction_rate = np.zeros(np.shape(Y))#flame.net_production_rates[:,np.newaxis]
            neg_reaction_rate =np.zeros(np.shape(Y))#flame.destruction_rates[:,np.newaxis]
            pos_reaction_rate = np.zeros(np.shape(Y))#net_reaction_rate- neg_reaction_rate
            cp_i = np.reshape(flame.partial_molar_cp/gas.molecular_weights, [gas.n_species, 1])
            enth_i = np.reshape(flame.partial_molar_enthalpies/gas.molecular_weights, [gas.n_species, 1])
            grid = np.zeros([1,1])
            velocity = np.zeros([1,1])
        else:
            Y = flame.Y
            X = flame.X
            net_reaction_rate = flame.net_production_rates
            neg_reaction_rate =flame.destruction_rates
            pos_reaction_rate = flame.net_production_rates + neg_reaction_rate
            cp_i = (flame.partial_molar_cp.T/gas.molecular_weights)
            enth_i = (flame.partial_molar_enthalpies.T/gas.molecular_weights)
            grid= flame.grid
            velocity = flame.velocity[:,np.newaxis]
        Y = Y.T
        try:
            mixture_fraction = flame.mixture_fraction("Bilger")
        except:
            mixture_fraction = np.sum(Y.T * np.reshape(self.z_i, [self.gas.n_species, 1]), axis=0) + self.c

        mean_molar_weights = np.dot(molar_weights.T, X)
        enthalpy = flame.enthalpy_mass

        density = flame.density
        cp = flame.cp_mass
        k = flame.thermal_conductivity

        T = flame.T

        viscosity = flame.viscosity

        Y_dot_net = net_reaction_rate * molar_weights
        Y_dot_pos = pos_reaction_rate * molar_weights
        Y_dot_neg = -neg_reaction_rate * molar_weights / (Y.T+1e-11)  # negated: SU2 uses source_prod + source_cons*Y directly

        Le_i = ComputeLewisNumber(flame)
        if self.__transport_model == "unity-Lewis-number":
            Le_i = Le_i / Le_i

        cp_i = np.reshape(cp_i, np.shape(Y))
        enth_i = np.reshape(enth_i, np.shape(Y))

        Le_i = Le_i.T

        if flame_is_gas:
            Le_i = np.reshape(Le_i, [1, gas.n_species])

        if flame_is_gas:
            heat_rel = 0.0
        else:
            heat_rel = flame.heat_release_rate

        # Define variables and output data array.
        variables = 'Distance,'
        data_matrix = np.reshape(grid, [len(grid), 1])
        variables += 'Velocity,'
        data_matrix = np.append(data_matrix, velocity,axis=1)
        variables += ','.join("Y-"+s for s in gas.species_names)
        data_matrix = np.append(data_matrix, Y,axis=1)
        variables += ',' + ','.join("Y_dot_net-"+s for s in gas.species_names)
        data_matrix = np.append(data_matrix, Y_dot_net.T, axis=1)
        variables += ',' + ','.join("Y_dot_pos-"+s for s in gas.species_names)
        data_matrix = np.append(data_matrix, Y_dot_pos.T, axis=1)
        variables += ',' + ','.join("Y_dot_neg-"+s for s in gas.species_names)
        data_matrix = np.append(data_matrix, Y_dot_neg.T, axis=1)
        variables += ',' + ','.join("Cp-"+s for s in gas.species_names)
        data_matrix = np.append(data_matrix, cp_i, axis=1)
        variables += ',' + ','.join("h-"+s for s in gas.species_names)
        data_matrix = np.append(data_matrix, enth_i, axis=1)
        variables += ',' + ','.join("Le-"+s for s in gas.species_names)
        data_matrix = np.append(data_matrix, Le_i, axis=1)
        variables += ',' + ','.join("X-"+s for s in gas.species_names)
        data_matrix = np.append(data_matrix, X.T, axis=1)

        if flame_is_gas:
            variables += ','+DefaultSettings_FGM.name_enth+','
            data_matrix = np.append(data_matrix, np.array([[enthalpy]]), axis=1)
            variables += DefaultSettings_FGM.name_mixfrac+','
            data_matrix = np.append(data_matrix, np.array([mixture_fraction]), axis=1)
            variables += '%s,' % FGMVars.Temperature.name
            data_matrix = np.append(data_matrix, np.array([[T]]), axis=1)
            variables += '%s,' % FGMVars.Density.name
            data_matrix = np.append(data_matrix, np.array([[density]]), axis=1)
            variables += '%s,' % FGMVars.MolarWeightMix.name
            data_matrix = np.append(data_matrix, mean_molar_weights.T, axis=1)
            variables += '%s,' % FGMVars.Cp.name
            data_matrix = np.append(data_matrix, np.array([[cp]]), axis=1)
            variables += '%s,' % FGMVars.Conductivity.name
            data_matrix = np.append(data_matrix, np.array([[k]]), axis=1)
            variables += '%s,' % FGMVars.ViscosityDyn.name
            data_matrix = np.append(data_matrix, np.array([[viscosity]]), axis=1)
            variables += '%s' % FGMVars.Heat_Release.name
            data_matrix = np.append(data_matrix, np.array([[heat_rel]]), axis=1)
        else:
            variables += ','+DefaultSettings_FGM.name_enth+','
            data_matrix = np.append(data_matrix, np.reshape(enthalpy, [len(enthalpy),1]), axis=1)
            variables += DefaultSettings_FGM.name_mixfrac+','
            data_matrix = np.append(data_matrix, np.reshape(mixture_fraction, [len(mixture_fraction),1]), axis=1)
            variables += '%s,' % FGMVars.Temperature.name
            data_matrix = np.append(data_matrix, np.reshape(T, [len(T), 1]), axis=1)
            variables += '%s,' % FGMVars.Density.name
            data_matrix = np.append(data_matrix, np.reshape(density, [len(density), 1]), axis=1)
            variables += '%s,' % FGMVars.MolarWeightMix.name
            data_matrix = np.append(data_matrix, mean_molar_weights.T, axis=1)
            variables += '%s,' % FGMVars.Cp.name
            data_matrix = np.append(data_matrix, np.reshape(cp, [len(cp), 1]), axis=1)
            variables += '%s,' % FGMVars.Conductivity.name
            data_matrix = np.append(data_matrix, np.reshape(k, [len(k), 1]), axis=1)
            variables += '%s,' % FGMVars.ViscosityDyn.name
            data_matrix = np.append(data_matrix, np.reshape(viscosity, [len(viscosity), 1]), axis=1)
            variables += '%s' % FGMVars.Heat_Release.name
            data_matrix = np.append(data_matrix, np.reshape(heat_rel, [len(heat_rel), 1]), axis=1)

        return variables, data_matrix

    def __TranslateToMatlabFile(self, filename:str, filename_out:str, output_dir:str):
        """Translate default FlameletAI output file to TableMaster compatible file.

        :param filename: default FlameletAI output file name.
        :type filename: str
        :param filename_out: output file name.
        :type filename_out: str
        :param output_dir: folder where to store the translated file.
        :type output_dir: str
        """
        fid = open(filename, "r")
        variables = fid.readline().strip().split(',')
        fid.close()

        data_flamelet = np.loadtxt(filename,delimiter=',',skiprows=1)

        species_in_flamelet = []
        species_molecular_weights = []
        for v in variables:
            if v[:2] == 'Y-':
                species_in_flamelet.append(v[2:])
                species_molecular_weights.append(self.gas.molecular_weights[self.gas.species_index(v[2:])])

        variables_1 = ['Distance',\
            'Temperature',\
            'Density',\
            'Conductivity',\
            'Dynamic_Viscosity',\
            'Cp',\
            'Total_Enthalpy',\
            'Heat_Release',\
            'Mixture_Fraction']

        variables_translated = ['Distance',\
                                'T',\
                                'rho',\
                                'Conductivity',\
                                'ViscosityDyn',\
                                'cp',\
                                'Enthalpy total',\
                                'Heat release rate',\
                                'Mixture Fraction']

        units = ['m',\
                'K', \
                'kg m^-3',\
                'W/m/K',\
                'kg/m/s',\
                'J/kg/K',\
                'J/kg',\
                'W/m^3',\
                '-']

        fid = open(output_dir + "/" + filename_out, 'w+')
        fid.write("Cantera (Bosch edit) flamelet\n\n")
        fid.write("Molecular weights:\n")
        fid.write(",".join(species_in_flamelet) + "\n")
        fid.write(",".join([str(m) for m in species_molecular_weights]) + "\n\n")
        fid.write(",".join([variables_translated[i] + " ("+units[i]+")" for i in range(len(variables_translated))]) + ",")
        fid.write(",".join(["Y-"+s for s in species_in_flamelet]) + ",")
        fid.write(",".join(["ReacRatePos-"+s for s in species_in_flamelet]) + ",")
        fid.write(",".join(["ReacRateNeg-"+s for s in species_in_flamelet]) + ",")
        fid.write(",".join(["cp-"+s for s in species_in_flamelet]) + ",")
        fid.write(",".join(["Enthalpy-"+s for s in species_in_flamelet]) + ",")
        fid.write(",".join(["Le-"+s for s in species_in_flamelet]))

        fid.write('\n\n')
        fid.close()

        idx_vars = [variables.index(v) for v in variables_1]
        idx_massfrac = [variables.index("Y-"+s) for s in species_in_flamelet]
        idx_pos_reacrate = [variables.index("Y_dot_pos-"+s) for s in species_in_flamelet]
        idx_neg_reacrate = [variables.index("Y_dot_neg-"+s) for s in species_in_flamelet]
        idx_cp_sp = [variables.index("Cp-"+s) for s in species_in_flamelet]
        idx_h_sp = [variables.index("h-"+s) for s in species_in_flamelet]
        idx_le_sp = [variables.index("Le-"+s) for s in species_in_flamelet]

        thermophysical_props = data_flamelet[:, [i for i in idx_vars]]
        massfracs = data_flamelet[:, [i for i in idx_massfrac]]
        pos_reacrate = data_flamelet[:, [i for i in idx_pos_reacrate]] / np.array([species_molecular_weights])
        neg_reacrate = data_flamelet[:, [i for i in idx_neg_reacrate]] / np.array([species_molecular_weights])
        cp_sp = data_flamelet[:, [i for i in idx_cp_sp]]
        h_sp = data_flamelet[:, [i for i in idx_h_sp]]
        le_sp = data_flamelet[:, [i for i in idx_le_sp]]

        total_data = np.hstack([thermophysical_props,\
                            massfracs,\
                            pos_reacrate,\
                            neg_reacrate,\
                            cp_sp,\
                            h_sp,le_sp])

        with open(output_dir + "/" + filename_out, "a+") as fid:
            csvWriter = csv.writer(fid)
            csvWriter.writerows(total_data)

def ComputeFlameletData(Config:Config_FGM, run_parallel:bool=False, N_processors:int=2, loglevel:int=0,
                        free_flame_refine:dict=None, burner_flame_refine:dict=None,
                        counter_flame_refine:dict=None):
    """Generate flamelet data according to Config_FGM settings either in serial or parallel.

    :param Config: Config_FGM class containing manifold and flamelet generation settings.
    :type Config: Config_FGM
    :param run_parallel: Generate flamelet data in parallel, defaults to False
    :type run_parallel: bool, optional
    :param N_processors: Number of parallel jobs when generating flamelet data in parallel, defaults to 0
    :type N_processors: int, optional
    :param loglevel: Cantera solver verbosity level (0=silent), defaults to 0
    :type loglevel: int, optional
    :param free_flame_refine: Cantera free-flame refinement criteria dict with keys ratio, slope, curve, prune.
        If None, the DataGenerator_Cantera defaults are used.
    :type free_flame_refine: dict, optional
    :param burner_flame_refine: Cantera burner-flame refinement criteria dict with keys ratio, slope, curve, prune.
        If None, the DataGenerator_Cantera defaults are used.
    :type burner_flame_refine: dict, optional
    :param counter_flame_refine: Cantera counter-flow flame refinement criteria dict with keys ratio, slope, curve, prune.
        If None, the DataGenerator_Cantera defaults are used.
    :type counter_flame_refine: dict, optional
    :raises Exception: If number of processors is set to zero when running in parallel.
    """

    if run_parallel and (N_processors == 0):
        raise Exception("Number of processors should be higher than zero when running in parallel.")

    mix_bounds = Config.GetMixtureBounds()
    Np_unb_mix = Config.GetNpMix()
    Config.gas.TP=300,101325
    Config.gas.set_equivalence_ratio(1.0, Config.GetFuelString(), Config.GetOxidizerString())
    if Config.GetMixtureStatus():
        mix_status_stoch = Config.gas.mixture_fraction(Config.GetFuelString(), Config.GetOxidizerString())
    else:
        mix_status_stoch = Config.gas.equivalence_ratio(Config.GetFuelString(), Config.GetOxidizerString())
    if mix_bounds[0] < mix_status_stoch and mix_bounds[1] > mix_status_stoch:
        mixture_range_lean = np.linspace(mix_bounds[0], mix_status_stoch, int(Np_unb_mix/2))
        mixture_range_rich = np.linspace(mix_status_stoch, mix_bounds[1], int(Np_unb_mix/2)+1)
        mixture_range = np.append(mixture_range_lean, mixture_range_rich[1:])
    else:
        # Equivalence ratios to calculate flamelets for are system inputs
        mixture_range = np.linspace(mix_bounds[0], mix_bounds[1], Np_unb_mix)

    def _make_generator():
        F = DataGenerator_Cantera(Config)
        F.SetLoglevel(loglevel)
        if free_flame_refine is not None:
            F.SetFreeFlameRefineCriteria(**free_flame_refine)
        if burner_flame_refine is not None:
            F.SetBurnerFlameRefineCriteria(**burner_flame_refine)
        if counter_flame_refine is not None:
            F.SetCounterFlameRefineCriteria(**counter_flame_refine)
        return F

    # Set up Cantera flamelet generator object

    def _ComputeFlameletData(mix_input):
        F = _make_generator()
        F.ComputeFlameletsOnMixStatus(mix_input)

    if run_parallel:
        with threadpool_limits(limits=1):
            Parallel(n_jobs=N_processors)(delayed(_ComputeFlameletData)(mix_status) for mix_status in mixture_range)
    else:
        F = _make_generator()
        F.SetMixtureValues(mixture_range)
        F.ComputeFlamelets()

def ComputeBoundaryData(Config:Config_FGM, run_parallel:bool=False, N_processors:int=2):

    def ComputeEquilibriumData(mix_input):
        F = DataGenerator_Cantera(Config)
        F.RunMixtureFraction()
        F.RunEquilibrium(True)
        F.RunFreeFlames(False)
        F.RunBurnerFlames(False)
        F.RunCounterFlowFlames(False)
        F.ComputeFlameletsOnMixStatus(mix_input)


    Np_unb_mix = Config.GetNpMix()
    Config.gas.TP=300,101325
    Config.gas.set_equivalence_ratio(1.0, Config.GetFuelString(), Config.GetOxidizerString())
    mix_status_stoch = Config.gas.mixture_fraction(Config.GetFuelString(), Config.GetOxidizerString())
    mixture_range_lean = np.linspace(0, mix_status_stoch, int(Np_unb_mix/2))
    mixture_range_rich = np.linspace(mix_status_stoch, 1, int(Np_unb_mix/2)+1)
    mixture_range = np.append(mixture_range_lean, mixture_range_rich[1:])

    if run_parallel:
        Parallel(n_jobs=N_processors)(delayed(ComputeEquilibriumData)(mix_status) for mix_status in mixture_range)
    else:
        for z in mixture_range:
            ComputeEquilibriumData(z)