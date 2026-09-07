###############################################################################################
#       #      _____ __  _____      ____        __        __  ____                   #        #
#       #     / ___// / / /__ \    / __ \____ _/ /_____ _/  |/  (_)___  ___  _____   #        #
#       #     \__ \/ / / /__/ /   / / / / __ `/ __/ __ `/ /|_/ / / __ \/ _ \/ ___/   #        #
#       #    ___/ / /_/ // __/   / /_/ / /_/ / /_/ /_/ / /  / / / / / /  __/ /       #        #
#       #   /____/\____//____/  /_____/\__,_/\__/\__,_/_/  /_/_/_/ /_/\___/_/        #        #
#       #                                                                            #        #
###############################################################################################

############################## FILE NAME: FlameletSolvers.py ##################################
#=============================================================================================#
# author: Evert Bunschoten                                                                    |
#    :PhD Candidate ,                                                                         |
#    :Flight Power and Propulsion                                                             |
#    :TU Delft,                                                                               |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#  Solvers used for flamelet data generation in SU2 DataMiner.                                |
#                                                                                             |
# Version: 3.2.0                                                                              |
#                                                                                             |
#=============================================================================================#
import cantera as ct
from os import sep, path, mkdir, getcwd
from typing import Dict
import numpy as np
import pandas as pd
from Common.DataDrivenConfig import Config_FGM
from Common.Properties import DefaultSettings_FGM, FGMVars
from Common.CommonMethods import ComputeLewisNumber

class FlameletSolver_Cantera:
    """Base class for Cantera flamelet solutions
    """
    _Config:Config_FGM = None
    _flamelet_type:str = "None"

    _flameletTypeOutputFolder:str = ""
    _plotLabel:str = None
    _flamelet_filename:str
    _is_premixed:bool = True
    _is_scalar:bool = False

    _canteraSolution:ct.Solution = None
    _flameletSolution:ct.FlameBase = None

    _iteration:int = 0
    _n_1D_iterations:int = 10
    _print_iteration:bool =False
    _keep_iterating:bool = True

    _from_restart:bool=False
    _from_file:bool=False
    _flameletSolutionForRestart:ct.FlameBase = None

    
    _T_reactants:float = DefaultSettings_FGM.T_min
    _reactant_mixture_status:float = 0.0
    _pressure:float = ct.one_atm

    # Initial grid and grid refinement options.
    _initial_grid_length:float = 1.8e-2
    _max_grid_length:float = 10.0
    __initial_grid_number_of_points:int = 30
    _initial_grid:np.ndarray[float] = None

    __grid_refinement_ratio:int = 3
    __grid_refinement_slope:float = 0.025
    __grid_refinement_curve:float = 0.025
    __grid_refinement_prune:float = 0.01
    _max_grid_points:int=2000

    __cantera_loglevel:int=0
    __flameletSolverLogLevel:int=1

    _flamelet_is_burning:bool = True
    _converged_solution:bool = True

    _output_filepath:str = getcwd()                 # System file location where flamelet data are saved.
    _thermochemical_solution:pd.DataFrame = None    # Flamelet solution data

    def __init__(self, config_input:Config_FGM):
        """Flamelet solver base class initializer.

        :param config_input: SU2 DataMiner configuration from which the reaction mechanism, transport mechanism, and output directory are retrieved.
        :type config_input: Config_FGM
        """
        self._Config = config_input
        self._canteraSolution = ct.Solution(self._Config.GetReactionMechanism())
        return
    
    def _initializeFlameletSolver(self):
        return
    
    def solveFor(self, **flamelet_solution_settings):
        """Calculate flamelet solution for custom settings.
        """
        self._parseInputSettings(flamelet_solution_settings)
        self.startSolver()
        return
    
    def solveAndSaveFor(self, **flamelet_solution_settings):
        """Calculate flamelet and save solution for custom settings.
        """
        self._parseInputSettings(flamelet_solution_settings)
        self.startSolver()
        self.saveFlameletSolution()
        return
    
    def solveForMixtureStatus(self, val_mixture_status:float, save:bool=True):
        """Calculate flamelet solutions for a specific mixture.

        :param val_mixture_status: mixture fraction or equivalence ratio.
        :type val_mixture_status: float
        :param save: store flamelet solutions, defaults to True
        :type save: bool, optional
        """
        self.setMixtureStatus(val_mixture_status)

        vals_input_settings = self._prepareSettingRange()
        self._iterateFlamelets(vals_input_settings, save)

        self.resetRestart()
        return

    def _iterateFlamelets(self, vals_input_settings, save):
        for i, q in enumerate(vals_input_settings):
            solver_specific_settings = self._writeSolverSettings(i, q)
            self._parseInputSettings(solver_specific_settings)
            self.startSolver()
            if save:
                self.saveFlameletSolution()
            if not self._keep_iterating:
                break
        return 
    
    def saveFlameletSolution(self):
        """Store flamelet solution in appropriately named folder.
        """
        self._prepareStorageFolder()
        self._writeOutput()
        return
    
    def _parseInputSettings(self, flamelet_solution_settings):
        return
    
    def _prepareStorageFolder(self):

        folder_for_flamelet_type = self._createFolderForFlameletType()
        mixture_subfolder = self.__createSubFolderForMixture()

        filepath_for_flamelet_data = sep.join((folder_for_flamelet_type, mixture_subfolder))
        if not path.isdir(filepath_for_flamelet_data):
            mkdir(filepath_for_flamelet_data)
        
        self._output_filepath = filepath_for_flamelet_data
        return
    
    def _createFolderForFlameletType(self):
        storage_folder = sep.join((self._Config.GetOutputDir(), self.getFlameletFolder()))
        if not path.isdir(storage_folder):
            mkdir(storage_folder)
        return storage_folder
    
    def __createSubFolderForMixture(self):
        if self._Config.GetMixtureStatus():
            tag_for_mixture_status="mixfrac"
        else:
            tag_for_mixture_status="phi"

        mixture_subfolder = "%s_%.4f" % (tag_for_mixture_status, self._reactant_mixture_status)
        return mixture_subfolder
    
    def startSolver(self):
        """Initiate flamelet simulation and report status to terminal.
        """
        self._prepareFlameletSimulation()
        self._computeFlameletSolution()
        self._postProcessResults()
        if self.__printToTerminal():
            self._printStatusToTerminal()
        return

    def _prepareFlameletSimulation(self):
        self._prepareReactants()
        self._fromRestart()
        self._solverSpecificPreprocessing()
        self._commonPreprocessing()
        return
    
    def _prepareFlameletSolver(self):
        if not self._from_file:
            if (not self._from_restart):
                self._initializeFlameletSolver()
                self._flameletSolution.max_grid_points = self._max_grid_points
            else:
                self._flameletSolution = self._flameletSolutionForRestart
        
        return
    
    def __printToTerminal(self):
        return self.__flameletSolverLogLevel > 0
    

    def _printStatusToTerminal(self):
        return
    

    def setReactantTemperature(self, Temp_reactants:float=DefaultSettings_FGM.T_min):
        """Specify the temperature of the reactants at the inflow boundary.

        :param Temp_reactants: reactant temperature value in Kelvin, defaults to 300K
        :type Temp_reactants: float, optional
        :raises Exception: if the specified value is not strictly positive.
        """
        if Temp_reactants <= 0:
            raise Exception("Reactant temperature should be strictly positive.")
        self._T_reactants = Temp_reactants
        return
    
    def setPressure(self, val_pressure:float=DefaultSettings_FGM.pressure):
        """Specify the pressure at which flamelet solutions are calculated.

        :param val_pressure: pressure in Pascals, defaults to 101325 Pa.
        :type val_pressure: float, optional
        :raises Exception: when specified value is not strictly positive.
        """
        if val_pressure <= 0:
            raise Exception("Pressure should be strictly positive")
        self._pressure = val_pressure
        return
    
    def setMixtureStatus(self, val_mixture_status:float):
        """Specify the equivalence ratio or mixture fraction for premixed flamelets.

        :param val_mixture_status: equivalence ratio or mixture fraction, depending on the mixture status definition in the config.
        :type val_mixture_status: float
        :raises Exception: if the mixture fraction is higher than one.
        :raises Exception: if the mixture status is not strictly positive.
        """
        if self._Config.GetMixtureStatus():
            if val_mixture_status > 1.0:
                raise Exception("Mixture fraction should be between zero and one.")
        if val_mixture_status < 0.0:
            raise Exception("Mixture status should be strictly positive.")
        
        self._reactant_mixture_status = val_mixture_status
        return
    
    def getReactantTemperature(self):
        return self._T_reactants
    
    def getMixtureStatus(self):
        return self._reactant_mixture_status
    
    def getFlameletFolder(self):
        return self._flameletTypeOutputFolder
    
    def getFlameletType(self):
        return self._flamelet_type
    
    def getPlotLabel(self):
        return self._plotLabel
    
    def setInitialGrid(self, initial_grid_length_in_meters:float=1.8e-2, number_of_nodes:int=100):
        """Specify the settings for the initial grid used to calculate the flamelet solution.

        :param initial_grid_length_in_meters: maximum grid length, defaults to 1.8e-2
        :type initial_grid_length_in_meters: float, optional
        :param number_of_nodes: number of nodes in the initial grid, defaults to 100
        :type number_of_nodes: int, optional
        :raises Exception: if the initial grid length is negative or exceeds the maximum grid length.
        :raises Exception: if the specified number of nodes is lower than 10.
        """
        if initial_grid_length_in_meters < 0 or initial_grid_length_in_meters >= self._max_grid_length:
            raise Exception("Initial grid length should be between 0 and %.1f." % self._max_grid_length)
        
        if number_of_nodes < 10:
            raise Exception("The initial grid should contain at least 10 nodes.")

        self._initial_grid_length = initial_grid_length_in_meters
        self.__initial_grid_number_of_points = number_of_nodes
        return
    
    def setGridRefinementCriteria(self, ratio:int=2, slope:float=0.025, curve:float=0.025, prune=0.01):
        """Specify the settings for refining the mesh of the flamelet simulation. See https://cantera.org/dev/reference/onedim/grid-refinement.html for more details.

        :param ratio: maximum ratio of grid length between adjacent elements, defaults to 2
        :type ratio: int, optional
        :param slope: maximum derivative in the solution between adjacent elements, defaults to 0.025
        :type slope: float, optional
        :param curve: maximum curvature in the solution between adjacent elements, defaults to 0.025
        :type curve: float, optional
        :param prune: remove grid nodes where curvature and slope falls below the specified value, defaults to 0.01
        :type prune: float, optional
        """
        self.__grid_refinement_ratio = ratio
        self.__grid_refinement_curve = curve
        self.__grid_refinement_prune = prune
        self.__grid_refinement_slope = slope
        return
    
    def getGridRefinementCriteria(self):
        return self.__grid_refinement_ratio, self.__grid_refinement_slope, self.__grid_refinement_curve, self.__grid_refinement_prune
    
    def setCanteraVerbose(self, verbose_level:int=0):
        """Specify verbosity of Cantera solution process.

        :param verbose_level: Cantera verbosity level, defaults to 0
        :type verbose_level: int, optional
        """
        self.__cantera_loglevel = verbose_level
        return
    
    def setSolverVerbose(self, verbose_level:int=1):
        """Specify the verbosity of the FlameletSolver solution process.

        :param verbose_level: FlameletSolver verbosity level, defaults to 1
        :type verbose_level: int, optional
        """
        self.__flameletSolverLogLevel = verbose_level
        return
    
    def _prepareSettingRange(self):
        self._print_iteration=True
        return
    
    def setInputVariable(self, val_input:float):
        return
    
    def getFlameletFileName(self, val_input:float):
        return ""
    
    def retrieveSolverSettings(self, solvers):
        return
    
    
    def _writeSolverSettings(self, solution_index:int, setting_1D:float):
        self._iteration = solution_index
        return
    
    
    
    def _writeOutput(self):
        if self.isConverged() and self.isBurning():
            self.saveSolution()
            self._flameletSolutionForRestart = self._flameletSolution
        return
    
    
    
    def _prepareReactants(self):
        if self._Config.GetMixtureStatus():
            self._canteraSolution.set_mixture_fraction(self._reactant_mixture_status, self._Config.GetFuelString(), self._Config.GetOxidizerString())
        else:
            self._canteraSolution.set_equivalence_ratio(self._reactant_mixture_status, self._Config.GetFuelString(), self._Config.GetOxidizerString())
        self._canteraSolution.TP = self._T_reactants, self._pressure
        return
    
    def _fromRestart(self):
        if self._flameletSolutionForRestart:
            self._from_restart = True
        else:
            self._from_restart = False
        return
    
    def _solverSpecificPreprocessing(self):
        self._setInitialGrid()
        self._prepareFlameletSolver()
        return
    
    def _setInitialGrid(self):
        self._initial_grid = np.linspace(0, self._initial_grid_length, self.__initial_grid_number_of_points)
        return
    
    def _commonPreprocessing(self):
        self._flameletSolution.set_refine_criteria(slope=self.__grid_refinement_slope,ratio=self.__grid_refinement_ratio,curve=self.__grid_refinement_curve,prune=self.__grid_refinement_prune)
        self._flameletSolution.transport_model = self._Config.GetTransportModel()
        return
    
    def _computeFlameletSolution(self):
        try:
            automatic_grid_refinement = (not self._from_restart)

            self._flameletSolution.solve(loglevel=self.__cantera_loglevel, refine_grid=True, auto=automatic_grid_refinement)
            self._converged_solution = True

            no_ignition = np.max(self._flameletSolution.T) <= DefaultSettings_FGM.T_threshold
            domain_too_long = (max(self._flameletSolution.grid) - min(self._flameletSolution.grid)) > self._max_grid_length
            too_few_nodes = len(self._flameletSolution.grid) < 10
            if no_ignition or domain_too_long or too_few_nodes:
                self._flamelet_is_burning = False
            else:
                self._flamelet_is_burning = True
        except:
            self._converged_solution = False
            self._flamelet_is_burning = False
        return
    
    def _postProcessResults(self):
        if self.isConverged():
            self._extractSolutionDataForOutput()
            self._flameletSolutionForRestart = self._flameletSolution

        if self._from_file:
            self._from_file = False
        return
    
    def getFlameletSolution(self):
        """Retrieve Cantera oneDim solution

        :return: Cantera flamelet solution object.
        :rtype: cantera.oneDim
        """
        return self._flameletSolution
    
    def isConverged(self):
        return self._converged_solution

    def isBurning(self):
        return self._flamelet_is_burning
    
    def getSolution(self):
        """Retrieve thermochemical state data extracted from flamelet solution.

        :return: flamelet solution data.
        :rtype: pandas.DataFrame
        """
        return self._thermochemical_solution
    
    def _extractSolutionDataForOutput(self):
        
        self._thermochemical_solution = pd.DataFrame()

        # Flamelet solution variables
        self._extractFlameletDiscretization()

        solution_1D = (np.asarray(self._flameletSolution.T).ndim == 1)

        # Species mass fractions data
        mass_fractions = self._flameletSolution.Y
        species_labels = ["Y-%s" % sp for sp in self._canteraSolution.species_names]
        self._concatenateThermoChemicalData(species_labels, mass_fractions)
        
        self._collectSourceTermData()

        if solution_1D:
            molecular_weights = self._canteraSolution.molecular_weights[:,np.newaxis]
        else:
            molecular_weights = self._canteraSolution.molecular_weights
        species_specific_heat = self._flameletSolution.partial_molar_cp / molecular_weights
        species_specific_enthalpy = self._flameletSolution.partial_molar_enthalpies / molecular_weights

        species_cp_labels = ["%s-%s" % (FGMVars.Cp.name, sp) for sp in self._canteraSolution.species_names]
        species_h_labels = ["h-%s" % (sp) for sp in self._canteraSolution.species_names]
        self._concatenateThermoChemicalData(species_cp_labels, species_specific_heat)
        self._concatenateThermoChemicalData(species_h_labels, species_specific_enthalpy)
        
        Le_i = ComputeLewisNumber(self._flameletSolution)
        if self._Config.GetTransportModel() == "unity-Lewis-number":
            Le_i = np.ones(Le_i.shape)

        Le_labels = ["Le-%s" % sp for sp in self._canteraSolution.species_names]
        self._concatenateThermoChemicalData(Le_labels, Le_i)

        # Enthalpy and mixture fraction
        total_enthalpy = self._flameletSolution.enthalpy_mass

        mixture_fraction_species_coefficients = self._Config.GetMixtureFractionCoefficients()
        if solution_1D:
            mixture_fraction_species_coefficients = mixture_fraction_species_coefficients[:,np.newaxis]
        mixture_fraction_offset = self._Config.GetMixtureFractionConstant()
        mixture_fraction = mixture_fraction_offset + np.sum(mixture_fraction_species_coefficients * mass_fractions,axis=0)

        controlvar_df = pd.DataFrame()
        controlvar_df[FGMVars.EnthalpyTot.name] = np.atleast_1d(total_enthalpy)
        controlvar_df[FGMVars.MixtureFraction.name] = np.atleast_1d(mixture_fraction)

        self._thermochemical_solution = pd.concat((self._thermochemical_solution, controlvar_df),axis=1)

        thermodynamic_df = pd.DataFrame()
        temperature = self._flameletSolution.T
        thermodynamic_df[FGMVars.Temperature.name] = np.atleast_1d(temperature)

        density = self._flameletSolution.density
        thermodynamic_df[FGMVars.Density.name] = np.atleast_1d(density)

        molar_fractions = self._flameletSolution.X
        mean_molecular_weight = np.sum(molecular_weights * molar_fractions, axis=0)
        thermodynamic_df[FGMVars.MolarWeightMix.name] = np.atleast_1d(mean_molecular_weight)

        specific_heat_cp = self._flameletSolution.cp_mass
        thermodynamic_df[FGMVars.Cp.name] = np.atleast_1d(specific_heat_cp)

        conductivity = self._flameletSolution.thermal_conductivity
        thermodynamic_df[FGMVars.Conductivity.name] = np.atleast_1d(conductivity)

        dynamic_viscosity = self._flameletSolution.viscosity
        thermodynamic_df[FGMVars.ViscosityDyn.name] = np.atleast_1d(dynamic_viscosity)

        if solution_1D:
            heat_release = self._flameletSolution.heat_release_rate
        else:
            heat_release = 0.0
        thermodynamic_df[FGMVars.Heat_Release.name] = np.atleast_1d(heat_release)

        self._thermochemical_solution = pd.concat((self._thermochemical_solution, thermodynamic_df),axis=1)

        self._writeInflowSettings()
        
        return

    def _collectSourceTermData(self):
        mass_fractions = self._flameletSolution.Y
        molecular_weights = self._canteraSolution.molecular_weights[:,np.newaxis]

        net_reaction_rate = self._flameletSolution.net_production_rates
        neg_reaction_rate = self._flameletSolution.destruction_rates
        pos_reaction_rate = net_reaction_rate - neg_reaction_rate
        Y_dot_net = net_reaction_rate * molecular_weights
        Y_dot_pos = pos_reaction_rate * molecular_weights
        Y_dot_neg = -neg_reaction_rate * molecular_weights / (mass_fractions+1e-11)

        net_rates_labels = ["Y_dot_net-%s" % sp for sp in self._canteraSolution.species_names]
        pos_rates_labels = ["Y_dot_pos-%s" % sp for sp in self._canteraSolution.species_names]
        neg_rates_labels = ["Y_dot_neg-%s" % sp for sp in self._canteraSolution.species_names]

        self._concatenateThermoChemicalData(net_rates_labels, Y_dot_net)
        self._concatenateThermoChemicalData(pos_rates_labels, Y_dot_pos)
        self._concatenateThermoChemicalData(neg_rates_labels, Y_dot_neg)
        return
    
    def _concatenateThermoChemicalData(self, labels:list[str], data:np.ndarray[float]):
        solution_1D = (np.asarray(self._flameletSolution.T).ndim == 1)
        if solution_1D:
            data_for_output = data.transpose()
        else:
            data_for_output = np.atleast_2d(data)
        output_df = pd.DataFrame()
        output_df[labels] = data_for_output
        self._thermochemical_solution = pd.concat((self._thermochemical_solution, output_df),axis=1)
        return
    
    def _writeInflowSettings(self):
        gas = ct.Solution(self._Config.GetReactionMechanism())
        if self._Config.GetMixtureStatus():
            reactant_mixture_fraction = self._reactant_mixture_status
            gas.set_mixture_fraction(self._reactant_mixture_status, self._Config.GetFuelString(), self._Config.GetOxidizerString())
            reactant_equivalence_ratio = gas.equivalence_ratio(self._Config.GetFuelString(), self._Config.GetOxidizerString())
        else:
            reactant_equivalence_ratio = self._reactant_mixture_status
            gas.set_equivalence_ratio(self._reactant_mixture_status, self._Config.GetFuelString(), self._Config.GetOxidizerString())
            reactant_mixture_fraction = gas.mixture_fraction(self._Config.GetFuelString(), self._Config.GetOxidizerString())
        
        inflow_df = pd.DataFrame()
        n_grid_points = len(self._thermochemical_solution[FGMVars.Temperature.name])
        inflow_df["ReactantMixtureFraction"] = reactant_mixture_fraction*np.ones(n_grid_points)
        inflow_df["ReactantEquivalenceRatio"] = reactant_equivalence_ratio*np.ones(n_grid_points)
        inflow_df["ReactantTemperature"] = self._T_reactants*np.ones(n_grid_points)
        self._thermochemical_solution = pd.concat((self._thermochemical_solution, inflow_df),axis=1)
        return
    
    def _extractFlameletDiscretization(self):
        grid= self._flameletSolution.grid
        self._thermochemical_solution["Distance"] = grid
        velocity = self._flameletSolution.velocity
        self._thermochemical_solution["Velocity"] = velocity
        return
    
    def saveSolution(self):
        flamelet_filename = self.getFlameletFileName()
        filename_plus_folder = sep.join((self._output_filepath, flamelet_filename))
        self._thermochemical_solution.to_csv(filename_plus_folder+".csv",index=False)
        return filename_plus_folder+".csv"
    
    def resetRestart(self):
        self._flameletSolutionForRestart = None
        self._from_restart = False
        return
    
    def loadSolution(self, flameletFileName:str):
        """Load flamelet solution data from csv file and initialize flamelet simulation from loaded data.

        :param flameletFileName: solution file path name.
        :type flameletFileName: str
        """
        self._from_file = True
        self._thermochemical_solution = pd.read_csv(flameletFileName)
        self._initial_grid = self._thermochemical_solution["Distance"]

        if self._Config.GetMixtureStatus():
            mixtureStatus = self._thermochemical_solution["ReactantMixtureFraction"][0]
        else:
            mixtureStatus = self._thermochemical_solution["ReactantEquivalenceRatio"][0]
        self.setMixtureStatus(mixtureStatus)
        reactantTemperature = self._thermochemical_solution["ReactantTemperature"][0]
        self.setReactantTemperature(reactantTemperature)
        self._prepareReactants()
        self._initializeFlameletSolver()

        initialGuessData = self._prepareInitialGuessData()
        try:
            self._flameletSolution.set_initial_guess(data=initialGuessData)
        except:
            print("Initializing the flamelet solution from data frame will be included in the upcoming Cantera version.")

        return
    
    def _prepareInitialGuessData(self):
        initialGuessData = pd.DataFrame()
        initialGuessData["temperature"] = self._thermochemical_solution[FGMVars.Temperature.name]
        initialGuessData["grid"] = self._thermochemical_solution[FGMVars.Distance.name]
        initialGuessData["density"] = self._thermochemical_solution[FGMVars.Density.name]
        initialGuessData["velocity"] = self._thermochemical_solution[FGMVars.Velocity.name]
        for sp in self._canteraSolution.species_names:
            initialGuessData["Y_%s" % sp] = self._thermochemical_solution["Y-%s" % sp]
        return initialGuessData
    
    def isPremixed(self):
        return self._is_premixed
    
    def isScalar(self):
        return self._is_scalar
    
    def getMassFractions(self):
        """Retrieve mass fraction data from flamelet solution.

        :return: data frame with species mass fractions.
        :rtype: pandas.DataFrame
        """
        Y = self._flameletSolution.Y
        species_names = self._canteraSolution.species_names
        struct = {sp : y for sp, y in zip(species_names, Y)}
        return pd.DataFrame(data=struct)
     
    
class FreeFlameSolver(FlameletSolver_Cantera):
    """Solver class for adiabatic free flamelets
    """
    __mass_flow_rate:float = 0.0

    def __init__(self, config_input:Config_FGM):
        """Adiabatic flamelet solver class. The reactant temperature and mixture status of the inflow boundary can be specified. The reaction mechanism and transport mechanism are retrieved from the configuration.

        :param config_input: SU2 DataMiner configuration.
        :type config_input: Config_FGM
        """
        FlameletSolver_Cantera.__init__(self, config_input)
        self._flameletTypeOutputFolder = "freeflame_data"
        self._flamelet_type = "Freeflame"
        self._plotLabel = "Adiabatic free flame"
        self._is_premixed = True
        self._is_scalar = False
        self._n_1D_iterations = self._Config.GetNpTemp()
        self.setGridRefinementCriteria(ratio=3, slope=0.03, curve=0.03, prune=0.01)
        self.setInitialGrid(0.2, 50)
        self._initializeFlameletSolver()
        return
    
    def _initializeFlameletSolver(self):
        self._flameletSolution = ct.FreeFlame(self._canteraSolution, grid=self._initial_grid)
        return
    
    def _solverSpecificPreprocessing(self):
        super()._solverSpecificPreprocessing()
        if not self._from_restart and not self._from_file:
            self._flameletSolution.set_initial_guess(locs=[0.0, 0.3, 0.5, 1.0])
        self._flameletSolution.inlet.T = self._T_reactants
        return
    
    def _prepareSettingRange(self):
        super()._prepareSettingRange()
        Tu_bounds = self._Config.GetUnbTempBounds()
        Tu_range = np.linspace(Tu_bounds[1], Tu_bounds[0], self._n_1D_iterations)
        return Tu_range
    
    def _parseInputSettings(self, flamelet_solution_settings):
        if "mixture_status" in flamelet_solution_settings.keys():
            self.setMixtureStatus(flamelet_solution_settings["mixture_status"])
        if "reactant_temperature" in flamelet_solution_settings.keys():
            self.setReactantTemperature(flamelet_solution_settings["reactant_temperature"])
        return super()._parseInputSettings(flamelet_solution_settings)
    
    def _writeSolverSettings(self, solution_index:int, setting_1D:float):
        super()._writeSolverSettings(solution_index, setting_1D)
        freeflame_settings = {"reactant_temperature":setting_1D, "solution_index":solution_index}
        return freeflame_settings
    
    def setInputVariable(self, val_input:float):
        self.setReactantTemperature(val_input)
        return
    
    def getFlameletFileName(self):
        if self._Config.GetMixtureStatus():
            tag_for_mixture_status="mixfrac"
        else:
            tag_for_mixture_status="phi"
        flamelet_filename = "%s_%s%.4f_Tu%.1f" % (self._flamelet_type, \
                                                        tag_for_mixture_status, \
                                                        self._reactant_mixture_status, \
                                                        self._T_reactants)
        return flamelet_filename
    
    def _printStatusToTerminal(self):
        if self._Config.GetMixtureStatus():
            tag_for_mixture_status="Z"
        else:
            tag_for_mixture_status="phi"
        if not self.isConverged():
            outp_message = "%s simulation at %s=%.3f, Tu=%.2f did not converge " % (self._flamelet_type,\
                                                                                  tag_for_mixture_status,\
                                                                                  self._reactant_mixture_status,\
                                                                                  self._T_reactants)
        elif not self.isBurning():
            outp_message = "%s at %s=%.3f, Tu=%.2f did not ignite " % (self._flamelet_type,\
                                                                    tag_for_mixture_status,\
                                                                    self._reactant_mixture_status,\
                                                                    self._T_reactants)
        else:
            outp_message = "Successful %s simulation at %s=%.3f, Tu=%.2f, Np=%i " % (self._flamelet_type,\
                                                                            tag_for_mixture_status,\
                                                                            self._reactant_mixture_status,\
                                                                            self._T_reactants,\
                                                                            self._thermochemical_solution.shape[0])
        
        if self._print_iteration:
            outp_message += "(%i/%i)" % (self._iteration+1, self._n_1D_iterations)
        print(outp_message)
        return
    
    def _postProcessResults(self):
        super()._postProcessResults()
        if self.isBurning() and self.isConverged():
            self.__mass_flow_rate = self._thermochemical_solution["Velocity"][0] * self._thermochemical_solution["Density"][0]
        return
    
    def getMassFlowRate(self):
        """Retrieve adiabatic mass flow rate.

        :return: mass flow rate evaluated at the inflow boundary.
        :rtype: float
        """
        return self.__mass_flow_rate
    
class BurnerFlameSolver(FlameletSolver_Cantera):
    """Solver class for burner-stabilized flamelets
    """
    __adiabatic_massflow:float = 0.0
    __val_massflow:float
    __val_massflow_enthalpy:float
    __delta_enth:float = None
    __delta_massflow:float = None
    __enth_prev:float = None

    def __init__(self, config_input:Config_FGM):
        FlameletSolver_Cantera.__init__(self, config_input)
        self._flameletTypeOutputFolder = "burnerflame_data"
        self._flamelet_type = "Burnerflame"
        self._plotLabel = "Burner-stabilized flame"
        self._is_premixed = True
        self._is_scalar = False
        deltaEnth_mdot = self._Config.GetMdotDHTarget()
        if deltaEnth_mdot > 0:
            self.__delta_enth = deltaEnth_mdot
        
        self._n_1D_iterations = self._Config.GetNpMdot()
        self.setReactantTemperature(self._Config.GetUnbTempBounds()[0])
        self.setGridRefinementCriteria(ratio=3, slope=0.15, curve=0.15, prune=0.05)
        self._initializeFlameletSolver()
        return
    
    def setReactantMassFlow(self, val_massflow_inlet:float):
        """Specify the mass flow rate at the inflow boundary.

        :param val_massflow_inlet: mass flux in kg /m /s
        :type val_massflow_inlet: float
        :raises Exception: if the specified value is negative.
        """
        if val_massflow_inlet < 0:
            raise Exception("Mass flow rate should be strictly positive.")
        self.__val_massflow = val_massflow_inlet
        return

    def setAdiabaticMassFlow(self, val_ad_massflow:float):
        self.__adiabatic_massflow = val_ad_massflow
        return 
    
    def getReactantMassFlow(self):
        return self.__val_massflow

    def _solverSpecificPreprocessing(self):
        super()._solverSpecificPreprocessing()
        self._flameletSolution.burner.mdot = self.__val_massflow
        self._flameletSolution.burner.T = self._T_reactants
        return
    
    def _initializeFlameletSolver(self):
        self._flameletSolution = ct.BurnerFlame(self._canteraSolution, self._initial_grid)

        return

    def setTargetEnthalpySpacing(self, val_delta_enth:float):
        """Specify the change in total enthalpy between adjacent flamelet solutions targeted by the solver.

        :param val_delta_enth: target value for the total enthalpy offset evalauted at the inflow boundary in Joules per kilo gram.
        :type val_delta_enth: float
        :raises Exception: if the specified value is equal to zero.
        """
        if val_delta_enth == 0:
            raise Exception("Enthalpy difference between adjacent flamelets should be higher than zero.")
        self.__delta_enth = abs(val_delta_enth)
        return 
    
    def retrieveSolverSettings(self, solvers:Dict[str, FlameletSolver_Cantera]):
        """Retrieve adiabatic mass flow rate from the adiabatic flamelet solver.

        :param solvers: dictionary of flamelet solvers for the current workflow
        :type solvers: Dict[str, FlameletSolver_Cantera]
        :raises Exception: if mass flow rate can not be retrieved.
        """
        if "FREEFLAME" in solvers.keys():
            freeflame_solver:FreeFlameSolver = solvers["FREEFLAME"]
            self.__adiabatic_massflow = freeflame_solver.getMassFlowRate()
        else:
            freeflame_solver:FreeFlameSolver = FreeFlameSolver(self._Config)
            freeflame_solver.setMixtureStatus(self._reactant_mixture_status)
            freeflame_solver.setReactantTemperature(self._Config.GetUnbTempBounds()[0])
            freeflame_solver.startSolver()
            if freeflame_solver.isBurning() and freeflame_solver.isConverged():
                self.__adiabatic_massflow = freeflame_solver.getMassFlowRate()
            else:
                raise Exception("Unable to calculate adiabatic mass flow rate")
        return
    
    def _loadSolverSpecificData(self):
        u = self._thermochemical_solution[FGMVars.Velocity.name][0]
        rho = self._thermochemical_solution[FGMVars.Density.name][0]
        self.setReactantMassFlow(u * rho)
        return
    
    def _prepareSettingRange(self):
        super()._prepareSettingRange()
        mdot_max = 0.98 * self.__adiabatic_massflow
        mdot_min = 0.001 * self.__adiabatic_massflow
        if self.__iterate_enthalpy():
            self.__delta_massflow = (mdot_max - mdot_min)/self._n_1D_iterations
            self.__val_massflow_enthalpy = mdot_max
        m_dot_range = np.linspace(mdot_max, mdot_min, self._n_1D_iterations+1)[:-1]
        return m_dot_range

    def _iterateFlamelets(self, vals_input_settings, save):
        if self.__iterate_enthalpy():
            i=0
            while self._keep_iterating:
                solver_specific_settings = self._writeSolverSettings(i, 0)
                self._parseInputSettings(solver_specific_settings)
                self.startSolver()
                if save:
                    self.saveFlameletSolution()
                i += 1
        else:
            super()._iterateFlamelets(vals_input_settings, save)
        return 
    
    def _postProcessResults(self):
        super()._postProcessResults()
        if self.__iterate_enthalpy():
            self.__update_mass_flow_rate()
        return
    
    def __iterate_enthalpy(self):
        return self.__delta_enth is not None

    def __update_mass_flow_rate(self):
        enth_current = self._thermochemical_solution[FGMVars.EnthalpyTot.name][0]
        if self.__enth_prev is not None:
            delta_enth = abs(enth_current - self.__enth_prev)
            scale_mdot = np.clip(self.__delta_enth/delta_enth, 0.2, 5.0)
            self.__delta_massflow *= scale_mdot

        self.__enth_prev = self._thermochemical_solution[FGMVars.EnthalpyTot.name][0]
        self.__val_massflow_enthalpy -= self.__delta_massflow
        self._keep_iterating = (self.__val_massflow_enthalpy > 0.001*self.__adiabatic_massflow)
        
        return 
    
    def setInputVariable(self, val_input:float):
        self.setReactantMassFlow(val_input)
        return
    
    def _parseInputSettings(self, flamelet_solution_settings):
        if "mixture_status" in flamelet_solution_settings.keys():
            self.setMixtureStatus(flamelet_solution_settings["mixture_status"])
        if "mdot" in flamelet_solution_settings.keys():
            self.setReactantMassFlow(flamelet_solution_settings["mdot"])
        if "burner_temperature" in flamelet_solution_settings.keys():
            self.setReactantTemperature(flamelet_solution_settings["burner_temperature"])

        if self.__iterate_enthalpy():
            self.setReactantMassFlow(self.__val_massflow_enthalpy)

        return super()._parseInputSettings(flamelet_solution_settings)
    
    def _writeSolverSettings(self, solution_index:int, setting_1D:float):
        super()._writeSolverSettings(solution_index, setting_1D)
        freeflame_settings = {"mdot":setting_1D}
        return freeflame_settings
    
    
    
    def getFlameletFileName(self):
        if self._Config.GetMixtureStatus():
            tag_for_mixture_status="mixfrac"
        else:
            tag_for_mixture_status="phi"
        flamelet_filename = "%s_%s%.4f_mdot%.3f" % (self._flamelet_type,\
                                                        tag_for_mixture_status,\
                                                        self._reactant_mixture_status,\
                                                        self.__val_massflow)
        return flamelet_filename
    
    def _printStatusToTerminal(self):
        if self._Config.GetMixtureStatus():
            tag_for_mixture_status="Z"
        else:
            tag_for_mixture_status="phi"
        if not self.isConverged():
            output_message = "%s simulation at %s=%.3f, mdot=%.3f kg/m2/s did not converge " % (self._flamelet_type,\
                                                                                  tag_for_mixture_status,\
                                                                                  self._reactant_mixture_status,\
                                                                                  self.__val_massflow)
        elif not self.isBurning():
            output_message = "%s at %s=%.3f, mdot=%.3f kg/m2/s did not ignite " % (self._flamelet_type,\
                                                                    tag_for_mixture_status,\
                                                                    self._reactant_mixture_status,\
                                                                    self.__val_massflow)
        else:
            output_message = "Successful %s simulation at %s=%.3f, mdot=%.3f kg/m2/s, Np=%i " % (self._flamelet_type,\
                                                                            tag_for_mixture_status,\
                                                                            self._reactant_mixture_status,\
                                                                            self.__val_massflow,\
                                                                            self._thermochemical_solution.shape[0])
        if self._print_iteration:
            output_message += "(%i/%i)" % (self._iteration+1, self._n_1D_iterations)

        print(output_message)
        return
    
    def loadSolution(self, flameletFileName):
        super().loadSolution(flameletFileName)
        u = self._thermochemical_solution[FGMVars.Velocity.name][0]
        rho = self._thermochemical_solution[FGMVars.Density.name][0]
        self.__val_massflow =  u*rho
        return

class EquilibriumSolver(FlameletSolver_Cantera):
    """Solver class for chemical equilibrium data
    """
    __is_reaction_products:bool = False
    __is_lean:bool = False
    __accumulated_solution:pd.DataFrame = pd.DataFrame()
    __solution_reactants:pd.DataFrame = pd.DataFrame()
    __solution_products:pd.DataFrame = pd.DataFrame()

    def __init__(self, config_input:Config_FGM):
        FlameletSolver_Cantera.__init__(self, config_input)
        self._flameletTypeOutputFolder = "equilibrium_data"
        self._flamelet_type = "Equilibrium"
        self._plotLabel = "Chemical equilibrium data"
        self._is_premixed = True
        self._is_scalar = True
        self._n_1D_iterations = self._Config.GetNpTemp()
        self._flameletFileExtension = "csv"
        return
    
    def solveForMixtureStatus(self, val_mixture_status:float, save:bool=True):
        self.__is_reaction_products = False
        self.__accumulated_solution = pd.DataFrame()
        super().solveForMixtureStatus(val_mixture_status, save)
        self.__solution_reactants = self.__accumulated_solution

        self.__is_reaction_products = True
        self.__accumulated_solution = pd.DataFrame()
        super().solveForMixtureStatus(val_mixture_status, save)
        self.__solution_products = self.__accumulated_solution

        self.resetRestart()
        return
    
    def _solverSpecificPreprocessing(self):
        self._flameletSolution = self._canteraSolution
        return
    
    def _initializeFlameletSolver(self):
        self._flameletSolution = self._canteraSolution
        return
    
    def _parseInputSettings(self, flamelet_solution_settings):
        if "mixture_status" in flamelet_solution_settings.keys():
            self.setMixtureStatus(flamelet_solution_settings["mixture_status"])
        if "reactant_temperature" in flamelet_solution_settings.keys():
            self.setReactantTemperature(flamelet_solution_settings["reactant_temperature"])
        return super()._parseInputSettings(flamelet_solution_settings)
   
    
    def _writeSolverSettings(self, solution_index:int, setting_1D:float):
        super()._writeSolverSettings(solution_index, setting_1D)
        freeflame_settings = {"reactant_temperature":setting_1D}
        return freeflame_settings
    

    def _prepareSettingRange(self):
        super()._prepareSettingRange()
        Tu_bounds = self._Config.GetUnbTempBounds()
        if self.__is_reaction_products:
            T_min = Tu_bounds[0]
            self._canteraSolution.TP = Tu_bounds[1], ct.one_atm
            enthalpy_max = self._canteraSolution.enthalpy_mass

            self._canteraSolution.TP = T_min, ct.one_atm
            if self.__is_lean:
                self._canteraSolution.equilibrate("TP")
            else:
                self._canteraSolution.equilibrate('HP')
            self._canteraSolution.HP = enthalpy_max, ct.one_atm
            T_max = self._canteraSolution.T
        else:
            T_min = Tu_bounds[0]
            T_max = Tu_bounds[1]

        T_range = np.linspace(T_min, T_max, self._n_1D_iterations)
        return T_range
    
    def setInputVariable(self, val_input:float):
        self.setReactantTemperature(val_input)
        return
    
    def _computeFlameletSolution(self):
        
        try:
            self._canteraSolution.TP = self._T_reactants, ct.one_atm
            self._converged_solution = True
        except:
            self._converged_solution = False
        return
    
    def _extractFlameletDiscretization(self):
        self._thermochemical_solution["Distance"] = np.zeros(1)
        self._thermochemical_solution["Velocity"] = np.zeros(1)
        return
    
    def concatenateSolution(self, other_solution:pd.DataFrame):
        concatenated_data = pd.concat((other_solution, self._thermochemical_solution),axis=0)
        return concatenated_data
    
    def reactionProducts(self, is_reaction_products:bool=True):
        self.__is_reaction_products = is_reaction_products
        return
    
    def isReactionProducts(self):
        return self.__is_reaction_products
    
    def _commonPreprocessing(self):
        return
    
    def getFlameletFileName(self):
        if self._Config.GetMixtureStatus():
            tag_for_mixture_status="mixfrac"
        else:
            tag_for_mixture_status="phi"
        
        if self.__is_reaction_products:
            file_header = "Products"
        else:
            file_header = "Reactants"
        flamelet_filename = "%s_%s%.4f" % (file_header, \
                                                tag_for_mixture_status, \
                                                self._reactant_mixture_status)
        return flamelet_filename
    
    def _postProcessResults(self):
        super()._postProcessResults()
        if self.isConverged():
            if self.__accumulated_solution.empty:
                self.__accumulated_solution = self._thermochemical_solution
            else:
                self.__accumulated_solution = pd.concat((self.__accumulated_solution, self._thermochemical_solution),axis=0)
        return
    
    def _writeOutput(self):
        if self.isConverged():
            self.saveSolution()
        return
    
    def saveSolution(self):
        flamelet_filename = self.getFlameletFileName()
        filename_plus_folder = sep.join((self._output_filepath, flamelet_filename))
        self.__accumulated_solution.to_csv("%s.%s" % (filename_plus_folder, self._flameletFileExtension),index=False)
        return
    
    def loadSolution(self, flameletFileName:str):
        solution_fileName = flameletFileName.split(sep)[-1]
        self.__is_reaction_products = ("Products" in solution_fileName)
        self.__accumulated_solution = pd.read_csv(flameletFileName)
        if self.__is_reaction_products:
            self.__solution_products = self.__accumulated_solution
        else:
            self.__solution_reactants = self.__accumulated_solution
        
        val_mixfrac = np.clip(self.__accumulated_solution[FGMVars.MixtureFraction.name].iloc[-1], 0.0, 1.0)
        self._canteraSolution.set_mixture_fraction(val_mixfrac, self._Config.GetFuelString(), self._Config.GetOxidizerString())
        if self._Config.DefineMixtureStatus():
            self._reactant_mixture_status = self._canteraSolution.mixture_fraction(self._Config.GetFuelString(), self._Config.GetOxidizerString())
        else:
            self._reactant_mixture_status = self._canteraSolution.equivalence_ratio(self._Config.GetFuelString(), self._Config.GetOxidizerString())
        
        self.setReactantTemperature(self.__accumulated_solution[FGMVars.Temperature.name][0])
        return
    
    def _fromRestart(self):

        return

    def _prepareReactants(self):
        super()._prepareReactants()
        if self.__is_reaction_products:
            Tu_bounds = self._Config.GetUnbTempBounds()
            T_min = Tu_bounds[0]
            self._canteraSolution.TP = T_min, ct.one_atm
            if self.__is_lean:
                self._canteraSolution.equilibrate("TP")
            else:
                self._canteraSolution.equilibrate("HP")
        
        self._canteraSolution.TP = self._T_reactants, self._pressure
        return
    
    def setMixtureStatus(self, val_mixture_status:float):
        super().setMixtureStatus(val_mixture_status)
        self.__is_lean = False
        if self._Config.GetMixtureStatus():
            z_stoch = self._Config.GetMixtureFractionConstant()
            if val_mixture_status<= z_stoch:
                self.__is_lean = True
        else:
            if val_mixture_status <= 1.0:
                self.__is_lean = True
        return

    def getSolution(self):
        return self.__accumulated_solution
    
    def getReactionProductData(self):
        return self.__solution_products
    def getReactantData(self):
        return self.__solution_reactants
    
    def _collectSourceTermData(self):
        Y_dot_net = np.zeros(self._canteraSolution.n_species)
        Y_dot_pos = np.zeros(self._canteraSolution.n_species)
        Y_dot_neg = np.zeros(self._canteraSolution.n_species)

        net_rates_labels = ["Y_dot_net-%s" % sp for sp in self._canteraSolution.species_names]
        pos_rates_labels = ["Y_dot_pos-%s" % sp for sp in self._canteraSolution.species_names]
        neg_rates_labels = ["Y_dot_neg-%s" % sp for sp in self._canteraSolution.species_names]

        self._concatenateThermoChemicalData(net_rates_labels, Y_dot_net)
        self._concatenateThermoChemicalData(pos_rates_labels, Y_dot_pos)
        self._concatenateThermoChemicalData(neg_rates_labels, Y_dot_neg)
        return
    
class CooledFlameInterpolator(FlameletSolver_Cantera):
    """Class for interpolated data between burner-stabilized flamelet and chemical equilibrium data
    """
    __burnerFlameSolution:pd.DataFrame = None
    __equilibriumSolution:pd.DataFrame = None

    def __init__(self, config:Config_FGM):
        super().__init__(config)
        self._flameletTypeOutputFolder = "interpolated_burnerflame_data"
        self._flamelet_type = "BurnerflameInt"
        self._plotLabel = "Interpolated burner flame"
        self._is_premixed = True
        self._is_scalar = True
        self._n_1D_iterations = self._Config.GetNpMdotExtra()
        self._flameletFileExtension = "csv"
        return
    
    def setEquilibriumData(self, eq_data:pd.DataFrame):
        self.__equilibriumSolution = eq_data
        return
    
    def setBurnerFlameData(self, burner_data:pd.DataFrame):
        self.__burnerFlameSolution = burner_data
        return
    
    def _prepareSettingRange(self):
        super()._prepareSettingRange()
        iter_values = [i for i in range(self._n_1D_iterations)]
        return iter_values
    
    def setInputVariable(self, val_input:float):
        self._iteration = val_input
        return
    
    def _computeFlameletSolution(self):
        ratio = float(self._iteration + 1) / float(self._n_1D_iterations)
        w_a_lin = 1.0 - ratio
        w_a_src = (1.0 - ratio) ** self._Config.GetSrcInterpExponent()
        w_b_src = 1.0 - w_a_src

        h_eq = self.__equilibriumSolution[FGMVars.EnthalpyTot.name]
        iMin = np.argmin(h_eq)

        exponential_interpolation = w_a_src * self.__burnerFlameSolution.values + w_b_src * self.__equilibriumSolution.values[iMin]

        linear_interpolation = w_a_lin * self.__burnerFlameSolution.values + ratio * self.__equilibriumSolution.values[iMin]

        self._thermochemical_solution = pd.DataFrame()
        for iVar, var in enumerate(self.__burnerFlameSolution.keys()):
            if (var==FGMVars.Heat_Release.name) or "Y_dot" in var:
                self._thermochemical_solution[var] = exponential_interpolation[:, iVar]
            else:
                self._thermochemical_solution[var] = linear_interpolation[:, iVar]
        return
    
    def _writeOutput(self):
        if self.isConverged():
            self.saveSolution()
        return
    
    def saveSolution(self):
        flamelet_filename = self.getFlameletFileName()
        filename_plus_folder = sep.join((self._output_filepath, flamelet_filename))
        self._thermochemical_solution.to_csv("%s.%s" % (filename_plus_folder, self._flameletFileExtension),index=False)
        return
    
    def _extractSolutionDataForOutput(self):
        return
    
    def getFlameletFileName(self):
        if self._Config.GetMixtureStatus():
            tag_for_mixture_status="mixfrac"
        else:
            tag_for_mixture_status="phi"
        flamelet_filename = "%s_%s%.4f_iter%i" % (self._flamelet_type, \
                                                        tag_for_mixture_status, \
                                                        self._reactant_mixture_status, \
                                                        self._iteration)
        return flamelet_filename
    
    def retrieveSolverSettings(self, solvers:Dict[str, FlameletSolver_Cantera]):
        burnerflameSolver = solvers["BURNERFLAME"]
        equilibriumSolver = solvers["EQUILIBRIUM"]
        self.setBurnerFlameData(burnerflameSolver.getSolution())
        self.setEquilibriumData(equilibriumSolver.getReactionProductData())
        return
    
    def _printStatusToTerminal(self):
        message_out = "Interpolated between burner-stabilized and equilibrium data "
        if self._print_iteration:
            message_out += "(%i/%i)" % (self._iteration+1, self._n_1D_iterations)
        print(message_out)
        return
    
    def _solverSpecificPreprocessing(self):
        return
    def _commonPreprocessing(self):
        return

    def _postProcessResults(self):
        return
    
    def _writeSolverSettings(self, solution_index:int, setting_1D:float):
        super()._writeSolverSettings(solution_index, setting_1D)
        return {"iteration":setting_1D}

    def loadSolution(self, flameletFileName:str):
        self._thermochemical_solution = pd.read_csv(flameletFileName)
        massfractions = [self._thermochemical_solution["Y-%s" % sp][0] for sp in self._canteraSolution.species_names]
        self._canteraSolution.Y = massfractions
        self._canteraSolution.TP = self._thermochemical_solution[FGMVars.Temperature.name][0], ct.one_atm
        self._canteraSolution.equilibrate("TP")

        if self._Config.DefineMixtureStatus():
            self._reactant_mixture_status = self._canteraSolution.mixture_fraction(self._Config.GetFuelString(), self._Config.GetOxidizerString())
        else:
            self._reactant_mixture_status = self._canteraSolution.equivalence_ratio(self._Config.GetFuelString(), self._Config.GetOxidizerString())
        self.setReactantTemperature(self._thermochemical_solution[FGMVars.Temperature.name][0])
        
        return

class CounterFlowDiffusionFlameSolver(FlameletSolver_Cantera):
    """Solver class for strained counter-flow diffusion flamelets
    """
    __strain_rate:float = 1.0
    __fuel_density:float = None
    __fuel_velocity:float = None
    __oxidizer_density:float = None
    __oxidizer_velocity:float = None

    def __init__(self, config_input):
        super().__init__(config_input)
        self._flameletTypeOutputFolder = "counterflame_data_n"
        self._flamelet_type = "CounterFlowDiffusionFlame"
        self._plotLabel = "Counter-flow diffusion flame"
        self._is_premixed = False
        self._is_scalar = False
        self.setInitialGrid(2e-1, 300)
        self.setGridRefinementCriteria(ratio=3, slope=0.04, curve=0.06, prune=0.02)
        return
    
    def _initializeFlameletSolver(self):
        self._flameletSolution = ct.CounterflowDiffusionFlame(self._canteraSolution, grid=self._initial_grid)
        return
    
    def _prepareFlameletSimulation(self):
        return super()._prepareFlameletSimulation()
    
    def _prepareReactants(self):
        super()._prepareReactants()
        self._canteraSolution.set_mixture_fraction(0.0, self._Config.GetFuelString(), self._Config.GetOxidizerString())
        self.__oxidizer_density = self._canteraSolution.density
        self._canteraSolution.set_mixture_fraction(1.0, self._Config.GetFuelString(), self._Config.GetOxidizerString())
        self.__fuel_density = self._canteraSolution.density

        L = self._initial_grid_length
        self.__fuel_velocity = self.__strain_rate * L / (1 + self.__fuel_density/self.__oxidizer_density)

        self.__oxidizer_velocity = self.__fuel_velocity * self.__fuel_density / self.__oxidizer_density
        return
    
    def _solverSpecificPreprocessing(self):
        super()._solverSpecificPreprocessing()
        self._flameletSolution.P = ct.one_atm

        self._flameletSolution.fuel_inlet.T = self._T_reactants
        self._flameletSolution.fuel_inlet.mdot = self.__fuel_density * self.__fuel_velocity
        self._flameletSolution.fuel_inlet.Y = self._Config.GetFuelString()

        self._flameletSolution.oxidizer_inlet.T = self._T_reactants
        self._flameletSolution.oxidizer_inlet.mdot = self.__oxidizer_density * self.__oxidizer_velocity
        self._flameletSolution.oxidizer_inlet.Y = self._Config.GetOxidizerString()
        return
    
    def setStrainRate(self, val_strain_rate:float=1.0):
        if val_strain_rate <= 0.0:
            raise Exception("Strain rate should be strictly positive")
        self.__strain_rate = val_strain_rate
        return
    
    def getFlameletFileName(self):
        flamelet_filename = "%s_strain%.3e_Tu%.1f" % (self._flamelet_type, \
                                                        self.__strain_rate, \
                                                        self._T_reactants)
        return flamelet_filename
    
    def _printStatusToTerminal(self):
        if not self.isConverged():
            output_message = "%s simulation at strain=%.2f s^-1 did not converge " % (self._flamelet_type,\
                                                                                  self.__strain_rate)
        elif not self.isBurning():
            output_message = "%s simulation at strain=%.2f s^-1 did not ignite " % (self._flamelet_type,\
                                                                                  self.__strain_rate)
        else:
            output_message = "Successful %s simulation at Tu=%.1f strain=%.2f s^-1, Np=%i " % (self._flamelet_type,\
                                                                            self._T_reactants,\
                                                                            self.__strain_rate,\
                                                                            self._thermochemical_solution.shape[0])
        if self._print_iteration:
            output_message += "(%i/%i)" % (self._iteration+1, self._n_1D_iterations)

        print(output_message)
        return
    
    def _writeInflowSettings(self):
        super()._writeInflowSettings()
        self._thermochemical_solution["StrainRate"] = self.__strain_rate
        self._thermochemical_solution["SpreadRate"] = self._flameletSolution.spread_rate
        return
    
    def loadSolution(self, flameletFileName:str):
        super().loadSolution(flameletFileName)
        grid_length = self._thermochemical_solution["Distance"].iloc[-1] - self._thermochemical_solution["Distance"].iloc[0]
        n_nodes = self._thermochemical_solution.shape[0]
        self.setInitialGrid(grid_length, n_nodes)
        self.__strain_rate = self._thermochemical_solution["StrainRate"][0]
        return
    
    def _prepareInitialGuessData(self):
        initialGuessData = super()._prepareInitialGuessData()
        initialGuessData["spreadRate"] = self._thermochemical_solution["SpreadRate"]
        return initialGuessData
    
    def _parseInputSettings(self, flamelet_solution_settings):
        if "strain_rate" in flamelet_solution_settings.keys():
            self.setStrainRate(flamelet_solution_settings["strain_rate"])
        if "reactant_temperature" in flamelet_solution_settings.keys():
            self.setReactantTemperature(flamelet_solution_settings["reactant_temperature"])
        return super()._parseInputSettings(flamelet_solution_settings)
   
    
    def _writeSolverSettings(self, solution_index:int, setting_1D:float):
        super()._writeSolverSettings(solution_index, setting_1D)
        freeflame_settings = {"reactant_temperature":setting_1D}
        return freeflame_settings
    

    def _prepareSettingRange(self):
        super()._prepareSettingRange()
        Tu_bounds = self._Config.GetUnbTempBounds()
        Tu_range = np.linspace(Tu_bounds[0], Tu_bounds[1], self._n_1D_iterations)
        return Tu_range
    
    def setInputVariable(self, val_input:float):
        self.setReactantTemperature(val_input)
        return
    
    def _prepareStorageFolder(self):

        folder_for_flamelet_type = self._createFolderForFlameletType()
        mixture_subfolder = self.__createSubFolderForStrain()

        filepath_for_flamelet_data = sep.join((folder_for_flamelet_type, mixture_subfolder))
        if not path.isdir(filepath_for_flamelet_data):
            mkdir(filepath_for_flamelet_data)
        
        self._output_filepath = filepath_for_flamelet_data
        return
    
    def __createSubFolderForStrain(self):
        strain_subfolder = "strain_%.3e" % (self.__strain_rate)
        return strain_subfolder
    
FlameletSolverDict:dict = {"FREEFLAME" : FreeFlameSolver,\
                           "BURNERFLAME" : BurnerFlameSolver,\
                           "EQUILIBRIUM" : EquilibriumSolver,\
                            "INT_BURNERFLAME" : CooledFlameInterpolator,\
                            "COUNTERFLAME" : CounterFlowDiffusionFlameSolver}