###############################################################################################
#       #      _____ __  _____      ____        __        __  ____                   #        #
#       #     / ___// / / /__ \    / __ \____ _/ /_____ _/  |/  (_)___  ___  _____   #        #
#       #     \__ \/ / / /__/ /   / / / / __ `/ __/ __ `/ /|_/ / / __ \/ _ \/ ___/   #        #
#       #    ___/ / /_/ // __/   / /_/ / /_/ / /_/ /_/ / /  / / / / / /  __/ /       #        #
#       #   /____/\____//____/  /_____/\__,_/\__/\__,_/_/  /_/_/_/ /_/\___/_/        #        #
#       #                                                                            #        #
###############################################################################################

############################# FILE NAME: LUTGenerators.py #####################################
#=============================================================================================#
# author: Evert Bunschoten                                                                    |
#    :PhD Candidate ,                                                                         |
#    :Flight Power and Propulsion                                                             |
#    :TU Delft,                                                                               |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#   Table generator classes for generating SU2-supported tables for FGM and NICFD problems    |
# Version: 3.2.0                                                                              |
#                                                                                             |
#=============================================================================================#

import numpy as np
import pandas as pd
import os
from Common.Properties import EntropicVars,DefaultSettings_FGM
from su2dataminer.generate_data import DataGenerator_CoolProp
from Common.DataDrivenConfig import Config_NICFD,Config_FGM
from Manifold_Generation.LUT.LUTGenerator_Base import SU2TableGenerator_Base
from Manifold_Generation.LUT.MeshTools import MeshThermodynamicPlane
class SU2TableGenerator_NICFD(SU2TableGenerator_Base):
    _Config:Config_NICFD = None
    __datagenerator:DataGenerator_CoolProp = None
    __thermodynamic_data:pd.DataFrame = None

    def __init__(self, config_in:Config_NICFD):
        self._state_quantities = [q.name for q in EntropicVars]
        super().__init__(config_in)
        self.__datagenerator = DataGenerator_CoolProp(self._Config)
        self.__datagenerator.SetFDStepSizes(3e-3,3e-3)
        self.__defaultTableVariables()
        return
    
    def __defaultTableVariables(self):
        vars_to_exclude = [EntropicVars.N_STATE_VARS.name]
        if not self._Config.TwoPhase():
            vars_to_exclude.append(EntropicVars.VaporQuality.name)
        if not self._Config.CalcTransportProperties():
            vars_to_exclude.append(EntropicVars.Conductivity.name)
            vars_to_exclude.append(EntropicVars.ViscosityDyn.name)
            vars_to_exclude.append(EntropicVars.VaporQuality.name)
        for var in EntropicVars:
            if var.name not in vars_to_exclude:
                self._table_vars.append(var.name)
        return
    
    def setTableVars(self, table_vars_in:list[str]):
        if self._Config.TwoPhase() and EntropicVars.VaporQuality.name in table_vars_in:
            print("Table generator not configured for two-phase, ignoring vapor quality from table data.")
            table_vars_in.remove(EntropicVars.VaporQuality.name)

        if not self._Config.CalcTransportProperties():
            if EntropicVars.Conductivity.name in table_vars_in:
                print("Table generator not configured for transport properties, ignoring conductivity data")
                table_vars_in.remove(EntropicVars.Conductivity.name)
            if EntropicVars.ViscosityDyn.name in table_vars_in:
                print("Table generator not configured for transport properties, ignoring viscosity data")
                table_vars_in.remove(EntropicVars.ViscosityDyn.name)

        valid_vars = True
        for v in table_vars_in:
            found_var = False
            for q in EntropicVars:
                if v.lower() == q.name.lower():
                    found_var = True
                    self._table_vars.append(q.name)
            if not found_var:
                print("Error, \"%s\" is not supported by SU2 DataMiner" % v)
                valid_vars = False
        if not valid_vars:
            raise Exception("Some specified thermophysical variables are not supported.")
        return super().setTableVars(table_vars_in)
    
    def _checkIfVariableIsValid(self, var_to_check:str):
        if var_to_check in self._state_quantities:
            return True
        else:
            return False
    
    def _passRefinementOptions(self, mesher:MeshThermodynamicPlane):
        if self._Config.TwoPhase():
            saturation_curve_pts_scaled = self.__calculateSaturationCurvePoints()
            mesher.setSaturationCurvePoints(saturation_curve_pts_scaled)

        return super()._passRefinementOptions(mesher)

    def __calculateSaturationCurvePoints(self):
        self.__datagenerator.GenerateSaturationCurveInterpolator()
        n_samples = 5000
        rho_saturation_curve = self.__datagenerator.ComputeSaturationCurve(N_samples=n_samples)[:,0]
        rho_min, rho_max = np.min(rho_saturation_curve), np.max(rho_saturation_curve)
        rho_saturation_curve = np.linspace(rho_min, rho_max, n_samples)
        e_saturation_curve = self.__datagenerator.GetSaturationCurveStaticEnergy(rho_saturation_curve)
        rhoe_saturation_curve = np.column_stack((rho_saturation_curve, e_saturation_curve))
        saturation_curve_pts_scaled = self._scaler_controlling_variables.transform(rhoe_saturation_curve)
        return saturation_curve_pts_scaled
    
    def _initiateMesher(self):
        return MeshThermodynamicPlane()
    
    def _getFluidDataForInterpolator(self):
        self.__datagenerator.PreprocessData()
        self.__datagenerator.ComputeData()
        state_data_pointcloud, valid_mask = self.__datagenerator.GetStateData()
        
        state_dataFrame = pd.DataFrame()
        for var in EntropicVars:
            if var.value != EntropicVars.N_STATE_VARS.value:
                state_dataFrame[var.name] = state_data_pointcloud[valid_mask, var.value]

        self.__thermodynamic_data = state_dataFrame

        return state_dataFrame
    
    def _createPointCloudForTableLevel(self, levelValue:float):

        rhoe_pointcloud = np.column_stack((self.__thermodynamic_data[EntropicVars.Density.name], self.__thermodynamic_data[EntropicVars.Energy.name]))
        const_z = np.zeros(len(rhoe_pointcloud))
        rhoe_pointcloud_scaled = self._scaler_controlling_variables.transform(rhoe_pointcloud)
        cv_pointcloud = np.column_stack((rhoe_pointcloud_scaled, const_z))
        return cv_pointcloud
    
    def _calculateTableStateData(self, cv_table_nodes:np.ndarray[float]):
        rhoe_table_nodes = self._scaler_controlling_variables.inverse_transform(cv_table_nodes[:, :2])
        state_data_out = np.zeros([len(rhoe_table_nodes), EntropicVars.N_STATE_VARS.value])
        for i, rhoe in enumerate(rhoe_table_nodes):
            try:
                self.__datagenerator.UpdateFluid(rhoe[0], rhoe[1])
                state_data, correct_phase = self.__datagenerator.GetStateVector()
                if correct_phase:
                    state_data_out[i] = state_data
            except:
                pass

        return state_data_out

    def _writeAdditionalInfoToTable(self, fid):
        fid.write("Fluid:\n")
        fid.write("%s\n" % self._Config.GetFluid())
        fid.write("Equation of state:\n")
        fid.write("%s\n" % self._Config.GetEquationOfState())
        if self._Config.TwoPhase():
            fid.write("Table contains two-phase data\n\n")
        else:
            fid.write("Table constains single-phase data\n\n")
        return
    
class SU2TableGenerator_FGM(SU2TableGenerator_Base):
    _Config:Config_FGM = None
    __refineGradients:bool = False
    _refine_for_gradients_of:list[str] = []
    _gradient_refinement_factors:list[float] = []
    _gradient_norm_factors:list[float] = []

    __refine_equilibrium:bool = False
    __equilibrium_refinement_factor:float = 0.5
    __margin_equilibrium:float = 2e-2

    def __init__(self, config_in:Config_FGM):
        super().__init__(config_in)
        self._getFluidDataForInterpolator()
        return
    
    def applyRefinementForGradientOf(self, varname:str, coef:float=0.5):
        """Scale table refinement based on the gradients of thermochemical properties.

        :param varname: quantity for which to evaluate gradients.
        :type varname: str
        :param coef: refinement coefficient, defaults to 0.5
        :type coef: float, optional
        :raises Exception: if quantity is not found in table variables.
        :raises Exception: if refinement coefficient value is negative.
        """
        if varname not in self._state_quantities:
            raise Exception("%s is not in the list of available thermophysical state variables" % varname)
        if coef <= 0:
            raise Exception("Refinement coeffcient should be positive")
        self._refine_for_gradients_of.append(varname)
        self._gradient_refinement_factors.append(coef)
        self.__refineGradients = True
        return
    
    def refineEquilibrium(self, coef:float=0.5, margin:float=2e-2):
        """Apply refinement to equilibrium areas of each table level.

        :param coef: refinement coefficient, defaults to 0.5
        :type coef: float, optional
        :raises Exception: if refinement coefficient value is negative.
        """
        if coef <= 0:
            raise Exception("Refinement coefficient should be positive")
        self.__refine_equilibrium = True
        self.__equilibrium_refinement_factor = coef
        self.__margin_equilibrium = margin
        return
    
    def _getFluidDataForInterpolator(self):
        flamelet_data_filename = os.sep.join((self._Config.GetOutputDir(), self._Config.GetConcatenationFileHeader()+"_full.csv"))
        flameletDataPointCloud = pd.read_csv(flamelet_data_filename)
        self._state_quantities = list(flameletDataPointCloud.keys())
        self._table_vars = self._state_quantities.copy()
        return flameletDataPointCloud
    
    def _processTableLevels(self):
        super()._processTableLevels()
        if self._is2D():
            mixture_status = self._Config.GetMixtureBounds()[0]
            if self._Config.DefineMixtureStatus():
                mixture_fraction = mixture_status
            else:
                self._Config.gas.set_equivalence_ratio(mixture_status, self._Config.GetFuelString(),self._Config.GetOxidizerString())
                mixture_fraction = self._Config.gas.mixture_fraction(self._Config.GetFuelString(),self._Config.GetOxidizerString())
            self._table_levels[0] = mixture_fraction
        return
    
    def _createPointCloudForTableLevel(self, levelValue:float):
        self._Config.gas.set_mixture_fraction(levelValue, self._Config.GetFuelString(),self._Config.GetOxidizerString())
        self._Config.gas.TP=self._Config.GetUnbTempBounds()[0],DefaultSettings_FGM.pressure
        h_min_unb = self._Config.gas.enthalpy_mass

        # Compute reactant progress variable for the current mixture fraction.
        pv_unb = self._Config.ComputeProgressVariable(variables=None, flamelet_data=None, Y_flamelet=self._Config.gas.Y[:,np.newaxis])[0]

        # Define maximum enthalpy as the reactant enthalpy at the maximum reactant temperature.
        self._Config.gas.TP=self._Config.GetUnbTempBounds()[1],DefaultSettings_FGM.pressure
        h_max = self._Config.gas.enthalpy_mass

        # Equilibrate at constant enthalpy to get product progress variable value.
        self._Config.gas.equilibrate("TP")
        pv_b = self._Config.ComputeProgressVariable(variables=None, flamelet_data=None, Y_flamelet=self._Config.gas.Y[:,np.newaxis])[0]

        # Define minimum enthalpy as the product enthalpy cooled to minimum reactant temperature.
        self._Config.gas.TP=self._Config.GetUnbTempBounds()[0],DefaultSettings_FGM.pressure
        h_min = self._Config.gas.enthalpy_mass

        # Define 2D grid between minimum and maximum progress variable and total enthalpy
        delta_pv_unb = 1e-2*(pv_b - pv_unb)
        pv_range = np.linspace(pv_unb-delta_pv_unb, pv_b, 100)
        h_range = np.linspace(h_min, h_max, 100)
        xgrid, ygrid = np.meshgrid(pv_range, h_range)
        zgrid = levelValue*np.ones(np.shape(xgrid))

        # 2: Locate nodes that are above the burner-stabilized enthalpy line
        CV_grid_init = np.vstack((xgrid.flatten(), ygrid.flatten(), zgrid.flatten())).transpose()
        pv_grid = CV_grid_init[:,0]
        h_grid = CV_grid_init[:,1]

        h_limit = ((h_min_unb - h_min) * pv_grid + (h_min*pv_unb - h_min_unb*pv_b))/(pv_unb - pv_b)
        idx_keep = h_grid >= h_limit

        cv_pointcloud = CV_grid_init[idx_keep]
        if self._is2D():
            cv_pointcloud_scaled = self._scaler_controlling_variables.transform(cv_pointcloud[:,:2])
            zcoords = np.zeros([len(cv_pointcloud), 1])
            pointCloudNodes = np.column_stack((cv_pointcloud_scaled, zcoords))
        else:
            cv_pointcloud_scaled = self._scaler_controlling_variables.transform(cv_pointcloud)
            pointCloudNodes = cv_pointcloud_scaled
        
        if self.__refineGradients:
            self.__computeGradientNormFactors(cv_pointcloud_scaled)

        return pointCloudNodes
    
    def __computeGradientNormFactors(self, cv_vals:np.ndarray[float]):
        gradients = self._fluid_data_interpolator.Jacobian(cv_vals, varnames=self._refine_for_gradients_of)
        grads_mag = np.linalg.norm(gradients,axis=0)
        self._gradient_norm_factors = np.max(grads_mag,axis=0)
        return
    
    def _refinelocation(self, x:float,y:float,z:float):
        ref_factor = super()._refinelocation(x,y,z)
        
        if self.__refineGradients:
            ref_factor_grads = self.__applyGradientRefinement(x,y,z)
            ref_factor = min(ref_factor, ref_factor_grads)
        
        if self.__refine_equilibrium:
            ref_factor_equilibrium = self.__applyEquilibriumRefinement(x,y,z)
            ref_factor = min(ref_factor, ref_factor_equilibrium)

        return ref_factor
    
    def __applyGradientRefinement(self, x:float, y:float, z:float):
        cv_input = np.array([x,y,z])
        jac = self._fluid_data_interpolator.Jacobian(cv_input[:self._nDim_table], varnames=self._refine_for_gradients_of)
        jac_mag = np.linalg.norm(jac,axis=0)
        jac_mag_norm = jac_mag / self._gradient_norm_factors
        ref_factor = 1.0
        for i, f in enumerate(self._gradient_refinement_factors):
            ref_factor = min(ref_factor, max(f, 1.0 - (1 - f)*jac_mag_norm[i]))
        return ref_factor
    
    def __applyEquilibriumRefinement(self, x:float, y:float, z:float):
        ref_factor = 1.0
        if self._is2D():
            if x <= self.__margin_equilibrium or x >= (1-self.__margin_equilibrium):
                ref_factor = self.__equilibrium_refinement_factor
        return ref_factor
    
    def _writeAdditionalInfoToTable(self, fid):
        fid.write("Fuel:\n")
        fid.write(",".join(["%.2f:%s" % (w, sp) for w, sp in zip(self._Config.GetFuelWeights(), self._Config.GetFuelDefinition())]))
        fid.write("\nOxidizer:\n")
        fid.write(",".join(["%.2f:%s" % (w, sp) for w, sp in zip(self._Config.GetOxidizerWeights(), self._Config.GetOxidizerDefinition())]))
        fid.write("\nProgress variable:\n")
        fid.write("+".join(["(%+.3e*%s)" % (w, sp) for w, sp in zip(self._Config.GetProgressVariableWeights(), self._Config.GetProgressVariableSpecies())]))
        fid.write("\n\n")
        return