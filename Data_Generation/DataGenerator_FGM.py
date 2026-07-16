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
from typing import Dict
import copy
from os import getcwd

from joblib import Parallel, delayed
from threadpoolctl import threadpool_limits

#---------------------------------------------------------------------------------------------#
# Importing DataMiner classes and functions
#---------------------------------------------------------------------------------------------#
from Common.DataDrivenConfig import Config_FGM
from Data_Generation.DataGenerator_Base import DataGenerator_Base
from Common.Properties import DefaultSettings_FGM

from Data_Generation.FlameletSolvers import FlameletSolverDict, FlameletSolver_Cantera, FreeFlameSolver, BurnerFlameSolver, EquilibriumSolver

class DataGenerator_Cantera(DataGenerator_Base):
    """Generate flamelet data using Cantera.

    :param Config: Config_FGM class describing the flamelet generation settings.
    :type: Config_FGM
    """
    # Generate flamelet data from Cantera computation.
    _Config:Config_FGM

    # Save directory for computed flamelet data
    __matlab__output_dir:str = getcwd()

    __fuel_string:str = ''
    __oxidizer_string:str = ''

    __tag_for_mixture_status:str = "phi"

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
    __flameletSolverDict:Dict[str, FlameletSolver_Cantera]

    __run_freeflames:bool = DefaultSettings_FGM.include_freeflames      # Run adiabatic flame computations
    __run_extra_interpolated_burnerflames:bool = True       # Run extra interpolated burner-stabilized flame computations
    __run_burnerflames:bool = DefaultSettings_FGM.include_burnerflames    # Run burner stabilized flame computations
    __run_equilibrium:bool = DefaultSettings_FGM.include_equilibrium    # Run chemical equilibrium computations
    __run_counterflames:bool = DefaultSettings_FGM.include_counterflames   # Run counter-flow diffusion flamelet simulations.

    __freeFlameSolver:FreeFlameSolver = None
    __burnerFlameSolver:BurnerFlameSolver = None
    __equilibriumSolver:EquilibriumSolver = None


    __freeflame_storage_folder:str = None
    __burnerflame_storage_folder:str = None
    __equilibrium_storage_folder:str = None
    __counterflame_storage_folder:str = None

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
        self.gas = ct.Solution(self._Config.GetReactionMechanism())

        self._Config.ComputeMixFracConstants()

        self.__flameletSolverDict = {}
        for flamelet_type in self._Config.getFlameletTypes():
            self.__flameletSolverDict[flamelet_type] = FlameletSolverDict[flamelet_type](self._Config)

        return

    def SetLoglevel(self, loglevel:int=0):
        """Set Cantera solver verbosity level (0=silent)."""
        self.__loglevel = loglevel
        # Propagate loglevel to all existing flamelet solvers
        for solver in self.__flameletSolverDict.values():
            solver.setCanteraVerbose(loglevel)
            solver.setSolverVerbose(loglevel)

    def SetInitialGridLength(self, length:float):
        """Set the initial domain length (in metres) for new flamelet grids."""
        self.__initial_grid_length = length

    def setRefinementCriteria(self, flamelet_type:str, ratio:float=3, slope:float=0.03, curve:float=0.03, prune:float=0.01):
        if flamelet_type not in self.__flameletSolverDict.keys():
            available = list(self.__flameletSolverDict.keys())
            raise Exception("%s is NOT included in the available flamelet types. Available: %s" % (flamelet_type, available))
        self.__flameletSolverDict[flamelet_type].setGridRefinementCriteria(ratio, slope, curve, prune)
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
        self._Config.SetNpTemp(n_flamelets_new)
        self.__SynchronizeSettings()
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
        self._Config.SetUnbTempBounds(T_unb_lower, T_unb_upper)
        self.__SynchronizeSettings()
        return

    def RunMixtureFraction(self):
        """Define the mixture status as mixture fraction instead of equivalence ratio.
        """
        self._Config.DefineMixtureStatus(True)
        self.__SynchronizeSettings()
        return

    def RunEquivalenceRatio(self):
        """Define the mixture status as equivalence ratio instead of mixture fraction.
        """
        self._Config.DefineMixtureStatus(False)
        self.__SynchronizeSettings()
        return

    def includeFlameletType(self, flamelet_type:str="FREEFLAME"):
        self._Config.includeFlameletType(flamelet_type)
        self.__SynchronizeSettings()

        return

    def excludeFlameletType(self, flamelet_type:str):
        self._Config.excludeFlameletType(flamelet_type)
        self.__SynchronizeSettings()

        return

    def getFlameletTypes(self):
        return self._Config.getFlameletTypes()

    def RunFreeFlames(self, input:bool=True):
        """Include adiabatic free-flame data in the manifold.

        :param input: Generate adiabatic free-flame data.
        :type input: bool
        """
        if input:
            self._Config.includeFlameletType("FREEFLAME")
        else:
            self._Config.excludeFlameletType("FREEFLAME")
        self.__SynchronizeSettings()
        return

    def RunBurnerFlames(self, input:bool=True):
        """Include burner-stabilized flame data in the manifold.

        :param input: Generate burner-stabilized flamelet data.
        :type input: bool
        """
        if input:
            self._Config.includeFlameletType("BURNERFLAME")
        else:
            self._Config.excludeFlameletType("BURNERFLAME")
        self.__SynchronizeSettings()
        return

    def RunEquilibrium(self, input:bool=True):
        """Include chemical equilibrium data in the manifold.

        :param input: Generate chemical equilibrium data.
        :type input: bool
        """
        if input:
            self._Config.includeFlameletType("EQUILIBRIUM")
        else:
            self._Config.excludeFlameletType("EQUILIBRIUM")
        self.__SynchronizeSettings()
        return


    def RunExtraInterpolatedBurnerFlames(self, input:bool=True):
        """Include extra interpolated burner-stabilized flame data in the manifold.

        :param input: Generate extra interpolated burner-stabilized flamelet data.
        :type input: bool
        """
        if input:
            self._Config.includeFlameletType("INT_BURNERFLAME")
        else:
            self._Config.excludeFlameletType("INT_BURNERFLAME")
        self.__SynchronizeSettings()
        return

    def SetNpMdotExtra(self, n_extra:int):
        """Set the number of interpolation steps for extra interpolated burner-stabilized flamelets.

        :param n_extra: Number of interpolated steps between the lowest-mdot burner flame and equilibrium.
        :type n_extra: int
        :raises Exception: if the provided number is lower than one.
        """
        if n_extra < 1:
            raise Exception("Number of extra interpolated mdot flamelets should be higher than one.")
        self._Config.SetNpMdotExtra(n_extra)
        return

    def SetMixtureValues(self, mixture_values:list[float]):
        """Set the reactant mixture status values manually.

        :param mixture_values: list of equivalence ratio or mixture fraction values.
        :type mixture_values: list[float]
        :raises Exception: If an empty list is provided.
        """
        if len(mixture_values) == 0:
            raise Exception("At least one mixture status value should be provided.")
        if any([phi < 0 for phi in mixture_values]):
            raise Exception("Mixture values should be strictly positive")

        self.__unb_mixture_status = [phi for phi in mixture_values]
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

    def SetOutputDir(self, output_dir_new:str):
        """Define the flamelet data output directory manually.

        :param output_dir_new: Flamelet data output directory.
        :type output_dir_new: str
        :raises Exception: If provided directory doesn't exist.
        """
        self._Config.SetOutputDir(output_dir=output_dir_new)
        self.__SynchronizeSettings()
        return

    def computeSingleFlamelet(self, flamelet_type:str, **solver_settings):
        flameletSolver:FlameletSolver_Cantera = self.__flameletSolverDict[flamelet_type]
        flameletSolver.solveAndSaveFor(**solver_settings)
        return

    def computePremixedFlameletsFor(self, mix_status:float):
        """Generate flamelet data for a given mixture fraction or equivalence ratio.

        :param mix_status: Mixture fraction or equivalence ratio value.
        :type mix_status: float
        :raises Exception: If mixture status value is below zero.
        """

        if mix_status < 0:
            raise Exception("Mixture status value should be positive.")

        correct_order = FlameletSolverDict.keys()
        for c in correct_order:
            if c in self.__flameletSolverDict.keys() and c != "COUNTERFLAME":
                flameletSolver = self.__flameletSolverDict[c]
                flameletSolver.retrieveSolverSettings(self.__flameletSolverDict)
                flameletSolver.solveForMixtureStatus(mix_status)

        return

    def computeNonPremixedFlamelets(self):
        if "COUNTERFLAME" in self.__flameletSolverDict.keys():
            self.__flameletSolverDict["COUNTERFLAME"].solveForMixtureStatus(0)
        return

    def ComputeFlamelets(self):
        """Generate and store all flamelet data for the current settings.
        """
        self.computeNonPremixedFlamelets()

        # Generate all other flamelet types.
        for mix_status in self.__unb_mixture_status:
            self.computePremixedFlameletsFor(mix_status)
        return

def ComputeFlameletData(Config:Config_FGM, run_parallel:bool=False, N_processors:int=2, loglevel:int=0,
                        free_flame_refine:dict=None, burner_flame_refine:dict=None, counter_flame_refine:dict=None):
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
    if Config.GetMixtureStatus():
        mix_status_stoch = Config.GetMixtureFractionConstant()
    else:
        mix_status_stoch = 1.0
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
            F.setRefinementCriteria("COUNTERFLAME", **counter_flame_refine)
        return F

    # Set up Cantera flamelet generator object

    def _ComputeFlameletData(mix_input):
        F = _make_generator()
        F.computePremixedFlameletsFor(mix_input)

    if run_parallel:
        F = DataGenerator_Cantera(Config)
        F.computeNonPremixedFlamelets()

        with threadpool_limits(limits=1):
            Parallel(n_jobs=N_processors)(delayed(_ComputeFlameletData)(mix_status) for mix_status in mixture_range)
    else:
        F = _make_generator()
        F.SetMixtureValues(mixture_range)
        F.ComputeFlamelets()

def ComputeBoundaryData(Config:Config_FGM, run_parallel:bool=False, N_processors:int=2):

    Config_local:Config_FGM = copy.copy(Config)
    Config_local.DefineMixtureStatus(True)
    flameletTypes_orig = Config_local.getFlameletTypes().copy()

    def ComputeEquilibriumData(mix_input):
        F = DataGenerator_Cantera(Config_local)
        F.includeFlameletType("EQUILIBRIUM")
        for flamelet_type in FlameletSolverDict.keys():
            if flamelet_type != "EQUILIBRIUM":
                F.excludeFlameletType(flamelet_type)

        F.computePremixedFlameletsFor(mix_input)


    Np_unb_mix = Config_local.GetNpMix()
    mix_status_stoch = Config_local.GetMixtureFractionConstant()
    mixture_range_lean = np.linspace(0, mix_status_stoch, int(Np_unb_mix/2))
    mixture_range_rich = np.linspace(mix_status_stoch, 1, int(Np_unb_mix/2)+1)
    mixture_range = np.append(mixture_range_lean, mixture_range_rich[1:])
    if run_parallel:
        Parallel(n_jobs=N_processors)(delayed(ComputeEquilibriumData)(mix_status) for mix_status in mixture_range)
    else:
        for z in mixture_range:
            ComputeEquilibriumData(z)

    Config.setFlameletTypes(flameletTypes_orig)