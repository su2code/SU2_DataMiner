###############################################################################################
#       #      _____ __  _____      ____        __        __  ____                   #        #
#       #     / ___// / / /__ \    / __ \____ _/ /_____ _/  |/  (_)___  ___  _____   #        #
#       #     \__ \/ / / /__/ /   / / / / __ `/ __/ __ `/ /|_/ / / __ \/ _ \/ ___/   #        #
#       #    ___/ / /_/ // __/   / /_/ / /_/ / /_/ /_/ / /  / / / / / /  __/ /       #        #
#       #   /____/\____//____/  /_____/\__,_/\__/\__,_/_/  /_/_/_/ /_/\___/_/        #        #
#       #                                                                            #        #
###############################################################################################

############################ FILE NAME: collectFlameletData.py ################################
#=============================================================================================#
# author: Evert Bunschoten                                                                    |
#    :PhD Candidate ,                                                                         |
#    :Flight Power and Propulsion                                                             |
#    :TU Delft,                                                                               |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#  Class for reading flamelet data files, extracting relevant data, and generating a          |
#  homogeneous distribution of flamelet data along the progress variable, enthalpy, and       |
#  mixture fraction direction.                                                                |
#                                                                                             |
# Version: 3.1.0                                                                              |
#                                                                                             |
#=============================================================================================#

import numpy as np
from copy import copy
import pandas as pd
from tqdm import tqdm
np.random.seed(0)
from random import sample
import os
import matplotlib.pyplot as plt
from collections.abc import Callable
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()['color']

from Common.DataDrivenConfig import Config_FGM
from Common.Properties import DefaultSettings_FGM, FGMVars, FGMPlotSymbols
from Data_Generation.FlameletSolvers import FlameletSolverDict, FlameletSolver_Cantera

class FlameletConcatenator:
    """Read, regularize, and concatenate flamelet data for MLP training or LUT generation.

    """
    __Config:Config_FGM = None # FlameletAI configuration for current workflow.

    __nFlamelets:int = 0
    __nFlameletDataPoints:int = 0

    __Np_per_flamelet:int = 2**DefaultSettings_FGM.batch_size_exponent          # Number of data points to extract per flamelet.
    __Np_equilibrium:int = 1            # Number of rows to read from each equilibrium file (1 = only the cold boundary point).
    __custom_resolution:bool = False    # Overwrite average number of data points per flamelet with a specified value.

    __mfrac_skip:int = 1        # Number of mixture status folder to skip while concatenating flamelet data.
    __flameletSolutionIndex:int = 0

    __ignore_mixture_bounds:bool = False
    __mix_status_max:float = DefaultSettings_FGM.eq_ratio_min     # Minimum mixture status value above which to collect flamelet data.
    __mix_status_min:float = DefaultSettings_FGM.eq_ratio_max    # Maximum mixture status value below which to collect flamelet data.

    __boundary_file_header:str = DefaultSettings_FGM.boundary_file_header

    __controlling_variables:list[str] = DefaultSettings_FGM.controlling_variables
    __N_control_vars:int = len(DefaultSettings_FGM.controlling_variables)

    # Thermodynamic data to search for in flamelet data.
    __TD_train_vars = [FGMVars.Temperature.name, \
                       FGMVars.MolarWeightMix.name, \
                       FGMVars.DiffusionCoefficient.name, \
                       FGMVars.Conductivity.name, \
                       FGMVars.ViscosityDyn.name, \
                       FGMVars.Cp.name]

    __TD_flamelet_data:np.ndarray = None

    # Differential diffusion data to search for in flamelet data.
    __PD_train_vars = [FGMVars.Beta_ProgVar.name, \
                       FGMVars.Beta_Enth_Thermal.name, \
                       FGMVars.Beta_Enth.name, \
                       FGMVars.Beta_MixFrac.name]

    __PD_flamelet_data:np.ndarray = None

    __flamelet_ID_vars = ['FlameletID']
    __flamelet_ID:np.ndarray = None

    # Passive species names for which to save production and consumption terms.
    __Species_in_FGM = []

    # Passive look-up terms to include in the manifold.
    __LookUp_vars = []
    __LookUp_flamelet_data:np.ndarray = None

    # Progress variable source term name.
    __PPV_train_vars = [FGMVars.ProdRateTot_PV.name]
    __Sources_vars = [__PPV_train_vars[0]]
    for s in __Species_in_FGM:
        __Sources_vars.append("Y_dot_pos-"+s)
        __Sources_vars.append("Y_dot_neg-"+s)
        __Sources_vars.append("Y_dot_net-"+s)
        __Sources_vars.append("Y-"+s)

    __Sources_flamelet_data = None

    __CV_flamelet_data:np.ndarray = None    # Flamelet controlling variables

    __write_LUT_data:bool = False   # Apply source term and equilibrium data corrections for table data preparations.

    __verbose:int=1

    def __init__(self, Config:Config_FGM,verbose_level:int=1):
        """Class constructor

        :param Config: loaded Config_FGM class for the current workflow.
        :type Config: Config_FGM
        """
        self.__verbose = verbose_level
        if self.__verbose >0:
            print("Loading flameletAI configuration " + Config.GetConfigName())
        self.__Config = copy(Config)
        self.__SynchronizeSettings()

        return

    def __SynchronizeSettings(self):
        self.__Np_per_flamelet = self.__Config.GetNpConcatenation()
        if self.__Np_per_flamelet is not None:
            self.__custom_resolution = True

        [self.__mix_status_min, self.__mix_status_max] = self.__Config.GetMixtureBounds()

        self.SetAuxilarySpecies(self.__Config.GetPassiveSpecies())
        self.SetLookUpVars(self.__Config.GetLookUpVariables())

        self.__controlling_variables = []
        for c in self.__Config.GetControllingVariables():
            self.__controlling_variables.append(c)
        self.__N_control_vars = len(self.__controlling_variables)

        return

    def ConcatenateFlameletData(self):
        """Read flamelets and concatenate relevant flamelet data in the appropriate resolution.
        """

        self.__SizeDataArrays()

        self.__extractFlameletData()

        self.__WriteOutputFiles()
        return

    def CollectBoundaryData(self):
        self.IgnoreMixtureBounds(True)
        self.__Config.setFlameletTypes(["EQUILIBRIUM"])
        self.__Config.SetConcatenationFileHeader(self.__boundary_file_header)
        self.ConcatenateFlameletData()
        return

    def __SizeDataArrays(self):
        """Size the output data arrays according to the number of flamelets and manifold resolution.
        """

        self.__printMsg("Counting number of flamelet solutions...")
        self.__countNumberofFlameletDataPoints()
        self.__printMsg("Done.")

        if not self.__custom_resolution:
            self.__Np_per_flamelet = int(self.__nFlameletDataPoints / self.__nFlamelets)

        self.__printMsg("Number of data-points per flamelet: %i " % self.__Np_per_flamelet)

        self.__CV_flamelet_data = np.zeros([self.__nFlamelets * self.__Np_per_flamelet, self.__N_control_vars])
        self.__TD_flamelet_data = np.zeros([self.__nFlamelets * self.__Np_per_flamelet, len(self.__TD_train_vars)])
        if self.__Config.PreferentialDiffusion():
            self.__PD_flamelet_data = np.zeros([self.__nFlamelets * self.__Np_per_flamelet, len(self.__PD_train_vars)])
        self.__Sources_flamelet_data = np.zeros([self.__nFlamelets * self.__Np_per_flamelet, 1 + 4 * len(self.__Species_in_FGM)])
        self.__LookUp_flamelet_data = np.zeros([self.__nFlamelets * self.__Np_per_flamelet, len(self.__LookUp_vars)])

        return

    def __countNumberofFlameletDataPoints(self):
        self.__nFlameletDataPoints = 0
        self.__nFlamelets = 0
        self.__loopOverFlamelets(self.__incrementNumberOfFlameletData)
        return

    def __loopOverFlamelets(self, taskPerFlamelet:Callable):
        flameletTypes = self.__Config.getFlameletTypes()
        for flameletType in flameletTypes:
            flameletSolver:FlameletSolver_Cantera = FlameletSolverDict[flameletType](self.__Config)
            storageFolder = os.sep.join((self.__Config.GetOutputDir(), flameletSolver.getFlameletFolder()))
            mixtures = os.listdir(storageFolder)
            mixtures.sort()
            self.__printMsg("Processing %s data..." % flameletSolver.getFlameletType())
            for phi in mixtures[::self.__mfrac_skip]:
                flameletsFolder = os.sep.join((storageFolder, phi))
                flameletFilesForMixture = os.listdir(flameletsFolder)
                flameletFilesForMixture.sort()
                for flameletFile in flameletFilesForMixture:
                    flameletFilePath = os.sep.join((flameletsFolder, flameletFile))
                    flameletSolver.loadSolution(flameletFilePath)
                    if self.__isWithinMixtureBounds(flameletSolver):
                        taskPerFlamelet(flameletSolver)
            self.__printMsg("Done.")
        return

    def __isWithinMixtureBounds(self, flameletSolver:FlameletSolver_Cantera):
        if self.__ignore_mixture_bounds or not flameletSolver.isPremixed():
            return True
        else:
            mixture_status = flameletSolver.getMixtureStatus()
            margin = 1e-2
            return (mixture_status-margin <= self.__mix_status_max) and (mixture_status+margin >= self.__mix_status_min)

    def __incrementNumberOfFlameletData(self, flameletSolver:FlameletSolver_Cantera):
        thermochemical_solution = flameletSolver.getThermoChemicalData()
        self.__nFlameletDataPoints += thermochemical_solution.shape[0]
        self.__nFlamelets += 1
        return

    def __extractFlameletData(self):
        self.__printMsg("Extracting thermochemical data from manifold...")
        self.__flameletSolutionIndex = 0
        self.__loopOverFlamelets(self.__interpolateAlongFlamelet)
        self.__printMsg("Done.")
        return

    def __interpolateAlongFlamelet(self, flameletSolution:FlameletSolver_Cantera):
        solutionData = self.__retrieveFlameletSolution(flameletSolution)
        flameletIsBurning = True
        if not flameletSolution.isScalar():
            T_flamelet = solutionData[FGMVars.Temperature.name]
            if np.max(T_flamelet) < DefaultSettings_FGM.T_threshold:
                flameletIsBurning = False

        if flameletIsBurning:
            tracingVariable = self.__calculateTracingVariable(solutionData)
            if self.__reactionProductsForLUT(flameletSolution):
                flameletIsValid = True
            else:
                flameletIsValid = max(tracingVariable) > 0

            if flameletIsValid:
                CV_data = self.__retrieveControlVariables(solutionData)
                TD_data = self.__retrieveThermoPhyiscalData(solutionData)
                LookUp_data = self.__retrievePassiveLookUpData(solutionData)
                Sources_data = self.__retrieveSourceTerms(solutionData)

                if flameletSolution.isScalar() or not flameletSolution.isPremixed():
                    tracingVariable_query = np.linspace(0, 1.0, self.__Np_per_flamelet)
                else:
                    tracingVariable_query = 0.5 - 0.5*np.cos(np.linspace(0, np.pi, self.__Np_per_flamelet))

                CV_data_interpolated = self.__interpolateAlongTracingVariable(tracingVariable, tracingVariable_query, CV_data)
                TD_data_interpolated = self.__interpolateAlongTracingVariable(tracingVariable, tracingVariable_query, TD_data)
                LookUp_data_interpolated = self.__interpolateAlongTracingVariable(tracingVariable, tracingVariable_query, LookUp_data)
                Sources_data_interpolated = self.__interpolateAlongTracingVariable(tracingVariable, tracingVariable_query, Sources_data)

                startIndex = self.__Np_per_flamelet*self.__flameletSolutionIndex
                stopIndex = self.__Np_per_flamelet*(self.__flameletSolutionIndex+1)
                self.__CV_flamelet_data[startIndex:stopIndex] = CV_data_interpolated
                self.__TD_flamelet_data[startIndex:stopIndex] = TD_data_interpolated
                self.__LookUp_flamelet_data[startIndex:stopIndex] = LookUp_data_interpolated
                self.__Sources_flamelet_data[startIndex:stopIndex] = Sources_data_interpolated
                if self.__Config.PreferentialDiffusion():
                    PD_data = self.__retrievePreferentialDiffusionScalars(solutionData)
                    PD_data_interpolated = self.__interpolateAlongTracingVariable(tracingVariable, tracingVariable_query, PD_data)
                    self.__PD_flamelet_data[startIndex:stopIndex] = PD_data_interpolated
        self.__flameletSolutionIndex += 1
        return

    def __retrieveFlameletSolution(self, flameletSolver:FlameletSolver_Cantera):
        solutionData = flameletSolver.getThermoChemicalData()

        if self.__reactionProductsForLUT(flameletSolver):
            solutionData = solutionData.iloc[:self.__Np_equilibrium,:]

        if not flameletSolver.isPremixed() and not self.__ignore_mixture_bounds:
            solutionData = self.__clipNonPremixedFlameletToMixtureBounds(solutionData)
        return solutionData

    def __clipNonPremixedFlameletToMixtureBounds(self, solutionData:pd.DataFrame):
        if self.__Config.GetMixtureStatus():
            mixfrac_upper = self.__mix_status_max
            mixfrac_lower = self.__mix_status_min
        else:
            phi_upper = self.__mix_status_max
            phi_lower = self.__mix_status_min
            self.__Config.gas.set_equivalence_ratio(phi_upper, self.__Config.GetFuelString(), self.__Config.GetOxidizerString())
            mixfrac_upper = self.__Config.gas.mixture_fraction(self.__Config.GetFuelString(), self.__Config.GetOxidizerString())
            self.__Config.gas.set_equivalence_ratio(phi_lower, self.__Config.GetFuelString(), self.__Config.GetOxidizerString())
            mixfrac_lower = self.__Config.gas.mixture_fraction(self.__Config.GetFuelString(), self.__Config.GetOxidizerString())
        mixfrac_solution = solutionData[FGMVars.MixtureFraction.name]
        within_bounds = np.logical_and(mixfrac_solution >= mixfrac_lower, mixfrac_solution <= mixfrac_upper)
        solutionData_out = solutionData.iloc[within_bounds]
        return solutionData_out

    def __WriteOutputFiles(self):
        """Collect all flamelet data arrays, split into train, test, and validation portions, and write to appropriately named files.
        """
        self.__printMsg("Writing accumulated flamelet data files...")

        flameletData = self.__accumulateDataFromFlamelets()

        flameletDataFull = self.__filterFlameletData(flameletData)

        trainData, testData, validationData = self.__splitFlameletDataSet(flameletDataFull)

        filename_full = os.sep.join((self.__Config.GetOutputDir(), self.__Config.GetConcatenationFileHeader()+"_full.csv"))
        flameletDataFull.to_csv(filename_full, index=False)

        filename_train = os.sep.join((self.__Config.GetOutputDir(), self.__Config.GetConcatenationFileHeader()+"_train.csv"))
        trainData.to_csv(filename_train, index=False)

        filename_test = os.sep.join((self.__Config.GetOutputDir(), self.__Config.GetConcatenationFileHeader()+"_test.csv"))
        testData.to_csv(filename_test, index=False)

        filename_validation = os.sep.join((self.__Config.GetOutputDir(), self.__Config.GetConcatenationFileHeader()+"_val.csv"))
        validationData.to_csv(filename_validation, index=False)

        self.__printMsg("Done.")
        return

    def __accumulateDataFromFlamelets(self):
        outputData = pd.DataFrame()
        outputData[self.__Config.GetControllingVariables()] = self.__CV_flamelet_data
        outputData[self.__TD_train_vars] = self.__TD_flamelet_data
        outputData[self.__LookUp_vars] = self.__LookUp_flamelet_data
        outputData[self.__Sources_vars] = self.__Sources_flamelet_data
        if self.__Config.PreferentialDiffusion():
            outputData[self.__PD_train_vars] = self.__PD_flamelet_data
        return outputData

    def __filterFlameletData(self, flameletData:pd.DataFrame):
        uniqueData = flameletData.drop_duplicates()
        noNans = uniqueData.dropna()
        noZeros = noNans.loc[~(noNans == 0).all(axis=1)].reset_index(drop=True)

        return noZeros

    def __splitFlameletDataSet(self, flameletDataFull:pd.DataFrame):
        shuffledData = flameletDataFull.sample(frac=1).reset_index(drop=True)
        Np_total = shuffledData.shape[0]
        Np_train = int(self.__Config.GetTrainFraction()*Np_total)
        Np_test = int(self.__Config.GetTestFraction()*Np_total)

        trainData = shuffledData.iloc[:Np_train,:]
        testData = shuffledData.iloc[Np_train:Np_train+Np_test,:]
        validationData = shuffledData.iloc[Np_train+Np_test:,:]
        return trainData, testData, validationData

    def IgnoreMixtureBounds(self, ignore_bounds:bool=False):
        self.__ignore_mixture_bounds = ignore_bounds
        return

    def WriteLUTData(self, write_LUT_data:bool=False):
        """Apply corrections to chemical equilibrium data and source terms in order to ensure boundary
        correctness for table generation.
        """
        self.__write_LUT_data = write_LUT_data
        return

    def SetNFlameletNodes(self, Np_per_flamelet:int):
        """Manually define the number of data points per flamelet to be included in the manifold.

        :param Np_per_flamelet: number of data points to be interpolated from each flamelet.
        :type Np_per_flamelet: int
        :raises Exception: if the number of points is lower than two.
        """
        if Np_per_flamelet < 2:
            raise Exception("Number of data points per flamelet should be higher than two.")
        self.__Np_per_flamelet = Np_per_flamelet
        self.__custom_resolution = True
        return

    def GetNFlameletNodes(self):
        return self.__Np_per_flamelet

    def SetNEquilibriumNodes(self, Np_equilibrium:int):
        """Define the number of rows to sample from each equilibrium file.

        A value of 1 (default) reads only the fully-cooled boundary point.
        Higher values sample the full equilibrium curve from cold to adiabatic,
        giving the table or MLP denser coverage of the equilibrium arm.

        :param Np_equilibrium: number of points per equilibrium file (>= 1).
        :type Np_equilibrium: int
        :raises Exception: if the value is lower than one.
        """
        if Np_equilibrium < 1:
            raise Exception("Number of equilibrium nodes must be at least 1.")
        self.__Np_equilibrium = min(Np_equilibrium, self.__Config.GetNpTemp())
        return

    def SetMixStep(self, skip_mixtures:int):
        """Skip a number of mixture status values when reading flamelet data to reduce the concatenated file size.

        :param skip_mixtures: step size in mixture status
        :type skip_mixtures: int
        :raises Exception: if the provided step size is lower than one.
        """
        if skip_mixtures < 1:
            raise Exception("Mixture step size should be higher than one.")
        self.__mfrac_skip = skip_mixtures
        return

    def SetMixStatusBounds(self, mix_status_low:float, mix_status_high:float):
        """Define the mixture status bounds between which to read flamelet data.

        :param mix_status_low: lower mixture status value.
        :type mix_status_low: float
        :param mix_status_high: upper mixture status value.
        :type mix_status_high: float
        :raises Exception: if the lower mixture status value is higher than the upper mixture status value.
        """
        if mix_status_low >= mix_status_high:
            raise Exception("Lower mixture status should be lower than upper value.")
        self.__mix_status_min = mix_status_low
        self.__mix_status_max = mix_status_high
        return

    def SetAuxilarySpecies(self, species_input:list[str]):
        """Define the passive species names for which to collect source terms.

        :param input: list of species names.
        :type input: list[str]
        """
        self.__Config.SetPassiveSpecies(species_input)

        self.__Species_in_FGM = []
        for s in species_input:
            self.__Species_in_FGM.append(s)
        self.__Sources_vars = [self.__PPV_train_vars[0]]
        for s in species_input:
            self.__Sources_vars.append("Y_dot_pos-"+s)
            self.__Sources_vars.append("Y_dot_neg-"+s)
            self.__Sources_vars.append("Y_dot_net-"+s)
            self.__Sources_vars.append("Y-"+s)
        return

    def SetControllingVariables(self, controlling_variables:list[str]=DefaultSettings_FGM.controlling_variables):
        self.__Config.SetControllingVariables(controlling_variables)
        self.__SynchronizeSettings()
        return

    def IncludeFlameletType(self, flameletType:str, include:bool=True):
        if include:
            self.__Config.includeFlameletType(flameletType)
        else:
            self.__Config.excludeFlameletType(flameletType)
        return

    def SetLookUpVars(self, input:list[str]):
        """Define passive look-up variables to be included in the manifold data.

        :param input: list of passive look-up variables.
        :type input: list[str]
        """
        self.__Config.SetLookUpVariables(input)
        self.__LookUp_vars = []
        for var in input:
            self.__LookUp_vars.append(var)
        return

    def SetFlameletDir(self, input:str):
        """Manually define the directory where the flamelet data is stored.

        :param input: path to flamelet data directory.
        :type input: str
        :raises Exception: if the provided directory doesn't exist.
        """
        self.__Config.SetOutputDir(input)
        return

    def SetOutputFileName(self, input:str):
        """Define the manifold output file header.

        :param input: manifold file header.
        :type input: str
        """
        self.__Config.SetConcatenationFileHeader(input)
        return

    def SetBoundaryFileName(self, input:str):
        self.__boundary_file_header = input
        return

    def SetTrainFraction(self, input:float=DefaultSettings_FGM.train_fraction):
        """Define the fraction of concatenated flamelet data to be used for training MLP's.

        :param input: train data fraction. Should be between zero and one.
        :type input: float
        :raises Exception: if the provided fraction is lower than zero or higher than one.
        """
        self.__Config.SetTrainFraction(input)
        return

    def SetTestFraction(self, input:float=DefaultSettings_FGM.test_fraction):
        """Define the fraction of concatenated flamelet data to be used for accuracy testing after training MLP's.

        :param input: test data fraction. Should be between zero and one.
        :type input: float
        :raises Exception: if the provided fraction is lower than zero or higher than one.
        """
        self.__Config.SetTestFraction(input)
        return


    def __reactionProductsForLUT(self, flameletSolver:FlameletSolver_Cantera):
        if self.__write_LUT_data and flameletSolver.getFlameletType()=="Equilibrium":
            return flameletSolver.isReactionProducts()
        else:
            return False

    def __calculateTracingVariable(self, solutionData:pd.DataFrame):
        controlVariables = self.__retrieveControlVariables(solutionData)

        cv_max, cv_min = np.max(controlVariables,axis=0), np.min(controlVariables,axis=0)
        scaledControlVariables = (controlVariables - cv_min)/(cv_max - cv_min + 1e-10)
        controlVariableIncrement = scaledControlVariables[1:] - scaledControlVariables[:-1]

        tracingVariableIncrement = np.linalg.norm(controlVariableIncrement,axis=1)
        tracingVariable = np.hstack((0, np.cumsum(tracingVariableIncrement)))
        tracingVariable_scaled = tracingVariable / (np.max(tracingVariable)+1e-10)
        return tracingVariable_scaled

    def __retrieveControlVariables(self, solutionData:pd.DataFrame):
        controlVariables = np.zeros([solutionData.shape[0], len(self.__Config.GetControllingVariables())])
        for iCv, cv in enumerate(self.__Config.GetControllingVariables()):
            if cv==DefaultSettings_FGM.name_pv:
                controlVariables[:, iCv] = self.__Config.ComputeProgressVariable(list(solutionData.keys()), solutionData.values)
            else:
                controlVariables[:, iCv] = solutionData[cv]
        return controlVariables

    def __retrieveThermoPhyiscalData(self, solutionData:pd.DataFrame):
        TD_data = np.zeros([solutionData.shape[0], len(self.__TD_train_vars)])
        for iVar_TD, TD_var in enumerate(self.__TD_train_vars):
            if TD_var == FGMVars.DiffusionCoefficient.name:
                conductivity = solutionData[FGMVars.Conductivity.name]
                cp = solutionData[FGMVars.Cp.name]
                density = solutionData[FGMVars.Density.name]
                TD_data[:, iVar_TD] = conductivity / (cp * density)
            else:
                TD_data[:, iVar_TD] = solutionData[TD_var]
        return TD_data

    def __retrievePassiveLookUpData(self, solutionData:pd.DataFrame):
        LookUp_data = np.zeros([solutionData.shape[0], len(self.__LookUp_vars)])
        for iVar_LookUp, LookUp_var in enumerate(self.__LookUp_vars):
            LookUp_data[:, iVar_LookUp] = solutionData[LookUp_var]
        return LookUp_data

    def __retrievePreferentialDiffusionScalars(self, solutionData:pd.DataFrame):
        vars = list(solutionData.keys())
        beta_pv_flamelet, beta_h1_flamelet, beta_h2_flamelet, beta_z_flamelet = self.__Config.ComputeBetaTerms(vars, solutionData.values)
        PD_data = np.zeros([solutionData.shape[0], len(self.__PD_train_vars)])
        PD_data[:, 0] = beta_pv_flamelet
        PD_data[:, 1] = beta_h1_flamelet
        PD_data[:, 2] = beta_h2_flamelet
        PD_data[:, 3] = beta_z_flamelet
        return PD_data

    def __retrieveSourceTerms(self, solutionData:pd.DataFrame):
        nP_flamelet = solutionData.shape[0]
        species_mass_fraction = np.zeros([nP_flamelet, len(self.__Species_in_FGM)])
        species_production_rate = np.zeros([nP_flamelet, len(self.__Species_in_FGM)])
        species_destruction_rate = np.zeros([nP_flamelet, len(self.__Species_in_FGM)])
        species_net_rate = np.zeros([nP_flamelet, len(self.__Species_in_FGM)])
        for iSp, Sp in enumerate(self.__Species_in_FGM):
            if Sp == "NOx":
                species_production_rate[:, iSp] = np.zeros(nP_flamelet)
                species_destruction_rate[:, iSp] = np.zeros(nP_flamelet)
                species_net_rate[:, iSp] = np.zeros(nP_flamelet)
                species_mass_fraction[:, iSp] = np.zeros(nP_flamelet)
                for NOsp in ["NO2","NO","N2O"]:
                    species_production_rate[:, iSp] += solutionData["Y_dot_pos-%s" % NOsp]
                    species_destruction_rate[:, iSp] += solutionData["Y_dot_neg-%s" % NOsp]
                    species_net_rate[:, iSp] += solutionData["Y_dot_net-%s" % NOsp]
                    species_mass_fraction[:, iSp] += solutionData["Y-%s" % NOsp]
            else:
                species_mass_fraction[:, iSp] = solutionData["Y-%s" % Sp]
                species_production_rate[:, iSp] = solutionData["Y_dot_pos-%s" % Sp]
                species_destruction_rate[:, iSp] = solutionData["Y_dot_neg-%s" % Sp]
                species_net_rate[:, iSp] = solutionData["Y_dot_net-%s" % Sp]

        Sources_data = np.zeros([nP_flamelet, 1 + 4 * len(self.__Species_in_FGM)])
        ppv_flamelet = self.__Config.ComputeProgressVariable_Source(list(solutionData.keys()), solutionData.values)
        Sources_data[:, 0] = ppv_flamelet

        for iSp in range(len(self.__Species_in_FGM)):
            Sources_data[:, 1 + 4*iSp] = species_production_rate[:, iSp]
            Sources_data[:, 1 + 4*iSp + 1] = species_destruction_rate[:, iSp]
            Sources_data[:, 1 + 4*iSp + 2] = species_net_rate[:, iSp]
            Sources_data[:, 1 + 4*iSp + 3] = species_mass_fraction[:, iSp]

        sourceterm_zero_line_numbers = np.zeros(nP_flamelet, dtype=bool)
        sourceterm_zero_line_numbers[0]  = True
        sourceterm_zero_line_numbers[-1] = True

        if self.__write_LUT_data:
            # Extend zeroing to a 2 % temperature-margin band at both ends.
            temp_margin = 2e-2
            T_flamelet = solutionData[FGMVars.Temperature.name]
            T_max, T_min = np.max(T_flamelet), np.min(T_flamelet)
            deltaT = temp_margin * (T_max - T_min)
            sourceterm_zero_line_numbers = np.logical_or(
                sourceterm_zero_line_numbers,
                np.logical_or((T_flamelet - T_min) < deltaT,
                                (T_max - T_flamelet) < deltaT))

        # Only zero the PV source term at the flamelet boundaries.
        # Species production rates and mass fractions are NOT zeroed
        # because (a) Y-{species} is non-zero at PV_max and (b) species
        # rates at PV_max are governed by Cantera's equilibrium values.
        Sources_data[sourceterm_zero_line_numbers, 0] = 0.0
        return Sources_data

    def __interpolateAlongTracingVariable(self, tracingVariableData:np.ndarray[float], tracingVariableQuery:np.ndarray[float], flameletData:np.ndarray[float]):
        flameletData_Sampled = np.zeros([self.__Np_per_flamelet, np.shape(flameletData)[1]])
        for i in range(np.shape(flameletData)[1]):
            flameletData_Sampled[:, i] = np.interp(tracingVariableQuery, tracingVariableData, flameletData[:, i])
        return flameletData_Sampled

    def __printMsg(self, msg:str):
        if self.__verbose > 0:
            print(msg)
        return

class GroupOutputs:
    """Class which groups flamelet data variables into MLP outputs based on their affinity.
    """

    __Config:Config_FGM = None    # FlameletAI configuration for the current problem.
    __controlling_variables:list[str] = DefaultSettings_FGM.controlling_variables
    __vars_to_exclude:list[str] = DefaultSettings_FGM.controlling_variables + ["FlameletID"]   # Variables to exclude from grouping; controlling variables by default.
    __flamelet_variables:list[str]  # Flamelet data variable names.

    __free_variables:list[str]      # Flamelet variables considered for grouping.
    __flamelet_data_filepath:str    # File path where flamelet data collection file is located.
    __flamelet_data:np.ndarray      # Concatenated flamelet data.
    __correlation_matrix:np.ndarray # Cross-correlation values between flamelet data variables.
    __iVar_remove:list[int]         # Variable indices to exclude from data set.

    __theta_threshold:float = 0.7   # Affinity threshold above which groups are accepted.

    __group_leaders_orig:list[str] = [] # Lead variables forced to represent separate groups.
    __n_groups:list[int]                # Number of groups in each combination.
    __group_variables:list[list[str]]   # Variables in each group.
    __group_affinity:list[list[float]]  # Minimum affinity for each combination of groups.

    # Combinations of variables for FGM evaluation.
    __evaluations_TD:list[str] = [FGMVars.Temperature.name, FGMVars.ViscosityDyn.name, FGMVars.MolarWeightMix.name, FGMVars.Cp.name, FGMVars.Conductivity.name, FGMVars.DiffusionCoefficient.name]
    __evaluations_PD:list[str] = [FGMVars.Beta_ProgVar.name, FGMVars.Beta_Enth_Thermal.name, FGMVars.Beta_Enth.name, FGMVars.Beta_MixFrac.name]
    __evaluations_Sources:list[str] = [FGMVars.ProdRateTot_PV.name]

    __most_interesting_groups:list[list[list[str]]] = []    # Combinations of groups with highest affinity for a certain group count.

    __best_group:int = 0
    def __init__(self, Config_in:Config_FGM):
        """Class constructor, load flamelet data.

        :param Config_in: FlameletAI configuration for the current problem.
        :type Config_in: Config_FGM
        """
        self.__Config = Config_in
        self.__flamelet_data_filepath = self.__Config.GetOutputDir()+"/"+self.__Config.GetConcatenationFileHeader()+"_full.csv"

        self.__controlling_variables = self.__Config.GetControllingVariables()
        self.__vars_to_exclude = []
        for var in self.__controlling_variables:
            self.__vars_to_exclude.append(var)
        self.__vars_to_exclude.append("FlameletID")
        self.__FilterVariables(self.__Config.GetControllingVariables() + ["FlameletID"])
        return

    def SetFlameletDataFile(self, filepath_in:str):
        """Define a custom flamelet data file for which to compute output groups.

        :param filepath_in: file path name of flamelet data file.
        :type filepath_in: str
        :raises Exception: if specified path does not exist on current hardware.
        """
        if not os.path.isdir(filepath_in):
            raise Exception("Supplied flamelet data file does not exist.")
        self.__flamelet_data_filepath = filepath_in
        self.__FilterVariables(self.__vars_to_exclude)

    def SetControllingVariables(self, control_vars:list[str]):
        """Define controlling variable names to always be excluded from flamelet data grouping.

        :param control_vars: list with controlling variable names
        :type control_vars: list[str]
        :raises Exception: if any of the specified variable names is not present in flamelet data set.
        """
        for var in control_vars:
            if var not in self.__flamelet_variables:
                raise Exception("Controlling variable " + var + " not present in flamelet data set.")

        vars_originally_excluded = self.__vars_to_exclude[len(self.__controlling_variables):]
        self.__controlling_variables = []
        self.__vars_to_exclude = []
        for var in control_vars:
            self.__controlling_variables.append(var)
            self.__vars_to_exclude.append(var)
        for var in vars_originally_excluded:
            self.__vars_to_exclude.append(var)

    def ExcludeVariables(self, vars_to_exclude:list[str]):
        """Add variables to be excluded from the output grouping.

        :param vars_to_exclude: list with variable names to be omitted from grouping.
        :type vars_to_exclude: list[str]
        :raises Exception: if any of the specified variables is not present in flamelet data set.
        """
        self.__vars_to_exclude = [c for c in self.__controlling_variables]
        self.__vars_to_exclude.append("FlameletID")
        for var in vars_to_exclude:
            if var not in self.__flamelet_variables:
                raise Exception("Variable "+var+" not present in flamelet data set.")
            self.__vars_to_exclude.append(var)
        return

    def SetAffinityThreshold(self, val_threshold:float=0.7):
        """Specify the threshold value for affinity below which groups are not considered.

        :param val_threshold: affinity threshold value. Should be between zero and one.
        :type val_threshold: float
        :raises Exception: if threshold value is not within range.
        """
        if val_threshold <= 0 or val_threshold >= 1:
            raise Exception("Threshold value should be between zero and one.")
        self.__theta_threshold = val_threshold

    def SetGroupLeaders(self, group_leaders_in:list[str]):
        """Specify a set of variables which are forced into separate groups.

        :param group_leaders_in: list of group leading variables.
        :type group_leaders_in: list[str]
        :raises Exception: if any of the variables is not present in the flamelet data set.
        """
        for g in group_leaders_in:
            if g not in self.__flamelet_variables:
                raise Exception("Variable " + g + " not present in flamelet data set.")
        self.__group_leaders_orig = []
        for g in group_leaders_in:
            self.__group_leaders_orig.append(g)

    def __FilterVariables(self, vars_to_remove:list[str]):
        with open(self.__flamelet_data_filepath, 'r') as fid:
            flamelet_variables = fid.readline().strip().split(',')
            self.__flamelet_variables = flamelet_variables
        self.__free_variables = []
        for var in self.__flamelet_variables:
            self.__free_variables.append(var)

        self.__iVar_remove = []
        for var in vars_to_remove:
            if var not in flamelet_variables:
                raise Exception("Variable " + var + " not present in flamelet data.")
            self.__iVar_remove.append(flamelet_variables.index(var))
            self.__free_variables.remove(var)

        self.__LoadFlameletData()
        self.__GenerateCorrelationMatrix()

        self.__correlation_matrix = np.delete(self.__correlation_matrix, self.__iVar_remove,0)
        self.__correlation_matrix = np.delete(self.__correlation_matrix, self.__iVar_remove,1)

    def __LoadFlameletData(self):
        self.__flamelet_data = np.loadtxt(self.__flamelet_data_filepath, delimiter=',',skiprows=1)

    def __GenerateCorrelationMatrix(self):
        self.__correlation_matrix = np.corrcoef(self.__flamelet_data.T)

    def __UpdateGroupLeaders(self, group_leaders_in:list[str]):
        group_variables = []
        group_indices = []
        group_affinity = []
        free_var_indices = [i for i in range(len(self.__correlation_matrix))]
        free_vars = [var for var in self.__free_variables]
        for g in group_leaders_in:
            group_indices.append([self.__free_variables.index(g)])
            group_affinity.append([1])
            group_variables.append([g])
            free_var_indices.remove(self.__free_variables.index(g))
            free_vars.remove(g)
        return group_variables, group_indices, group_affinity, free_var_indices, free_vars

    def __AffinityFunction(self, group_indices:list[int], iVar:int):
        theta = 1
        for k in group_indices:
                theta *= np.abs(self.__correlation_matrix[iVar,iVar])*np.abs(self.__correlation_matrix[iVar, k])
        return theta

    def EvaluateGroups(self):
        """Perform affinity evaluation and generate combinations of groups with a minimum affinity beyond the threshold value.
        """
        self.__FilterVariables(self.__vars_to_exclude)

        self.__group_affinity = []
        self.__group_variables = []
        self.__n_groups = []

        # Specify initial groups according to group leaders.
        group_variables, group_indices, group_affinity, _, free_vars_orig = self.__UpdateGroupLeaders(self.__group_leaders_orig)

        # Repeat 1000 times to come up with plenty of potential groups.
        for _ in tqdm(range(10000)):
            repeat = True

            while repeat:
                # Randomly select species from list of remaining species to act as additional group leaders.
                new_group_vars = sample(free_vars_orig, np.random.randint(1, len(free_vars_orig)))

                group_leaders = [g for g in self.__group_leaders_orig] + [g for g in new_group_vars]

                # Randomly select a species and add to an appropriate group by computing maximum affinity with that group.
                group_variables, group_indices, group_affinity, _, free_vars = self.__UpdateGroupLeaders(group_leaders)

                n_free_vars = len(free_vars)
                while n_free_vars > 0:
                    var_sample = sample(free_vars, 1)[0]
                    iVar = self.__free_variables.index(var_sample)
                    affinity_groups = []
                    for iGroup in range(len(group_leaders)):
                        theta = self.__AffinityFunction(group_indices=group_indices[iGroup], iVar=iVar)
                        affinity_groups.append(theta)
                    best_group_index = np.argmax(affinity_groups)
                    group_variables[best_group_index].append(var_sample)
                    group_indices[best_group_index].append(iVar)
                    group_affinity[best_group_index].append(max(affinity_groups))
                    free_vars.remove(var_sample)
                    n_free_vars -= 1
                repeat = False
                for g in group_affinity:
                    if min(g) < self.__theta_threshold:
                        repeat = True
            min_affinity = 1
            for g in group_affinity:
                min_affinity = min(min_affinity, min(g))
            n_groups = len(group_leaders)
            self.__n_groups.append(n_groups)
            for igroup, g in enumerate(group_variables):
                group_variables[igroup] = sorted(g)
            self.__group_variables.append(sorted(group_variables))
            self.__group_affinity.append(min_affinity)

        self.PostProcessGroups()
        return

    def __ComputeNumberofEvaluations(self, group_variables:list[list[str]]):
        n_networks_eval = 0

        for g in group_variables:
            this_group_TD = False
            this_group_PD = False
            this_group_sources = False
            for var in self.__evaluations_TD:
                if var in g:
                    this_group_TD = True
            if self.__Config.PreferentialDiffusion():
                for var in self.__evaluations_PD:
                    if var in g:
                        this_group_PD = True
            for var in self.__evaluations_Sources:
                if var in g:
                    this_group_sources = True

            if this_group_TD:
                n_networks_eval += 1
            if this_group_PD:
                n_networks_eval += 1
            if this_group_sources:
                n_networks_eval += 1
        return n_networks_eval

    def PostProcessGroups(self):
        """Extract the combinations of variables with the highest affinity and fewest number of network evaluations.
        Groups with most potential are visualized in a figure.
        """
        min_group = min(self.__n_groups)
        max_group = max(self.__n_groups)

        unique_groups = range(min_group, max_group + 1)

        n_network_evals = []
        interesting_groups = []
        self.__most_interesting_groups = []
        for j,i in enumerate(unique_groups):
            same_number_of_groups = np.argwhere(np.array(self.__n_groups) == i)[:,0]
            affinities_combinations = np.array(self.__group_affinity)[same_number_of_groups]
            iMax_affinity = np.argmax(affinities_combinations)
            interesting_group = self.__group_variables[same_number_of_groups[iMax_affinity]]
            n_network_evals.append(self.__ComputeNumberofEvaluations(interesting_group))
            interesting_groups.append(interesting_group)
            self.__most_interesting_groups.append(interesting_group)

        group_fewest_evaluations = np.argmin(np.array(n_network_evals))
        print("Output combinations with fewest number of network evaluations:")
        for iGroup, g in enumerate(interesting_groups[group_fewest_evaluations]):
            print("Output group " + str(iGroup)+ ": [" + ",".join("\"" + s + "\"" for s in g) + "]")
        self.__best_group = group_fewest_evaluations

    def GetInterestingGroup(self, iGroup:int=-1):
        """Get the group or groups with the highest efficiency.

        :param iGroup: combination index for which to display the output groups. If none provided, all combinations are returned.
        :type iGroup: int, optional
        :raises Exception: if specified index exceeds number of combinations.
        :return: list of output groups or list of combinations.
        :rtype: list[str]
        """
        if iGroup == -1:
            return self.__most_interesting_groups[self.__best_group]
        else:
            if iGroup >= len(self.__most_interesting_groups):
                raise Exception("Index exceeds number of best combinations")
            return self.__most_interesting_groups[iGroup]

    def PlotCorrelationMatrix(self, combination_index:int=-1):
        """Plots cross-correlation matrix between filtered flamelet data.

        :param combination_index: variable combination for which to plot output groups, defaults to 0
        :type combination_index: int, optional
        :raises Exception: if index exceeds number of combinations.
        """
        if combination_index == -1:
            combination_index = self.__best_group
        if combination_index >= len(self.__most_interesting_groups):
            raise Exception("Index exceeds number of best combinations")

        N=len(self.__most_interesting_groups[combination_index])
        plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.cubehelix(np.linspace(0,1,N+1)))
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        fig = plt.figure(figsize=[10,10])
        ax = plt.axes()
        ax.matshow(np.abs(self.__correlation_matrix), cmap="gray")
        for i in range(len(self.__free_variables)):
            for j in range(len(self.__free_variables)):
                ax.text(i, j, "%.2f" % (np.abs(self.__correlation_matrix[i,j])),\
                        fontsize=12,\
                        horizontalalignment='center',\
                        verticalalignment='center',\
                        color ='k' if np.abs(self.__correlation_matrix[i,j])>0.5 else 'white')

        for iGroup, g in enumerate(self.__most_interesting_groups[combination_index]):
            color = colors[iGroup]
            for iVar, v in enumerate(g):
                if iVar == 0:
                    ax.plot(self.__free_variables.index(v), self.__free_variables.index(v), 's',markerfacecolor='none',color=color,markersize=36, markeredgewidth=5,label="Group "+str(iGroup+1))

                else:
                    ax.plot(self.__free_variables.index(v), self.__free_variables.index(g[0]), 's',markerfacecolor='none',color=color,markersize=36, markeredgewidth=5)
                ax.text(self.__free_variables.index(v), len(self.__correlation_matrix), "%i" % (iGroup+1),fontsize=20,\
                        horizontalalignment='center',\
                        verticalalignment='center')
        ax.text(-0.5, len(self.__correlation_matrix), r"$J_\mathrm{group}$", fontsize=20,horizontalalignment='right',verticalalignment='center')
        ax.set_xticks(range(len(self.__free_variables)))
        ax.set_yticks(range(len(self.__free_variables)))


        ax.set_xticklabels([FGMPlotSymbols[q] for q in self.__free_variables])
        ax.set_yticklabels([FGMPlotSymbols[q] for q in self.__free_variables])
        ax.tick_params(axis='x',labelrotation=90)
        ax.tick_params(which='both',labelsize=18)
        fig.savefig(self.__Config.GetOutputDir()+"/Group_correlation_matrix.pdf",format='pdf',bbox_inches='tight')
        plt.tight_layout()
        plt.show()

        return

    def UpdateConfig(self, combination_index:int=-1):
        """Update the output groups in the FlameletAI configuration

        :param combination_index: group combination index to store in config, defaults to -1
        :type combination_index: int, optional
        """

        # By default, select the combination with the fewest number of function
        # evaluations.
        if combination_index == -1:
            combination_index = self.__best_group

        # Clear output groups present in the configuration.
        self.__Config.ClearOutputGroups()

        # Add flamelet variable groups to configuration.
        for group in self.__most_interesting_groups[combination_index]:
            self.__Config.AddOutputGroup(group)

        # Save configuration.
        self.__Config.SaveConfig()
        return