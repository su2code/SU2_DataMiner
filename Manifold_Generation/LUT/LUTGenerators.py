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
from Common.Properties import EntropicVars,DefaultSettings_FGM,FGMVars,PDClosureSourceName
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

    __hull_edge_tolerance:float = 1e-3

    # Dead bands, in normalized reaction progress, at which the preferential diffusion closure
    # sources are switched off at the unburnt and burnt ends of each level (see SetSourceDeadBands).
    __deadband_unburnt:float = 0.005
    __deadband_burnt:float = 0.005
    __deadband_nodes:list = None   # per level, the nodes the dead bands switched off (kept by the projection)

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
        """Apply refinement to the frozen progress variable edges of each table level.

        Both edges are covered: the fresh reactants and the equilibrium products. The margin is a
        fraction of the progress variable range of the level concerned, so the refined band
        follows the equilibrium edge across the levels of a 3D table. To refine only the burnt
        side, or a band wider than a margin around both edges, use
        :meth:`applyRefinementWithinNormalized` instead.

        :param coef: refinement coefficient, defaults to 0.5
        :type coef: float, optional
        :param margin: width of the refined band as a fraction of the level's progress variable range, defaults to 2e-2
        :type margin: float, optional
        :raises Exception: if refinement coefficient value is negative.
        :raises Exception: if the margin does not lie between zero and one half.
        """
        if coef <= 0:
            raise Exception("Refinement coefficient should be positive")
        if margin <= 0 or margin >= 0.5:
            raise Exception("Equilibrium refinement margin should lie between zero and one half")
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
        """Refine the two chemically frozen progress variable edges of the level being meshed.

        The margin is measured against the progress variable range of that level rather than
        against the globally scaled coordinate. On a 3D table those differ: the burnt edge is the
        equilibrium progress variable at the level's mixture fraction, which only reaches the
        global maximum on the level where equilibrium is most advanced. Comparing the raw scaled
        coordinate to the margin, as this did while the criterion was restricted to 2D tables,
        would refine the burnt edge of that one level and miss it everywhere else.
        """
        ref_factor = 1.0
        pv_normalized = self._normalizedLevelCoord(0, np.array([x, y, z]))
        if pv_normalized is None:
            return ref_factor
        if pv_normalized <= self.__margin_equilibrium or pv_normalized >= (1 - self.__margin_equilibrium):
            ref_factor = self.__equilibrium_refinement_factor
        return ref_factor
    
    def __closureSourceNames(self):
        """Closure source column names, one per controlling variable ("Res_<cv>")."""
        return [PDClosureSourceName(cv) for cv in self._Config.GetControllingVariables()]

    def __closureSourceVars(self):
        """Names of the tabulated preferential diffusion closure sources actually present."""
        if not self._is3D():
            return []
        if len(self._data_in_table) == 0 or self._data_in_table[0] is None:
            return []
        columns = list(self._data_in_table[0].keys())
        return [v for v in self.__closureSourceNames() if v in columns]

    def _varsExcludedFromSmoothing(self):
        """Leave the preferential diffusion closure sources out of the smoothing pass.

        These sources are sharp, sign-changing dipoles across the flame front whose defining
        property is that they integrate to zero along each flamelet. An RBF smoother spanning
        a hundred neighbours flattens that structure and destroys the cancellation: measured
        on the H2 manifold, smoothing raises the pre-projection conservation defect of S_Z by
        roughly a factor of six. Their noise is not what smoothing is meant to address.
        """
        excluded = super()._varsExcludedFromSmoothing()
        excluded += self.__closureSourceNames()
        return excluded

    def _unsmoothedNodes(self, level_index:int):
        """Leave the padding below the reactant progress variable out of the smoothing.

        Each level's mesh domain reaches a little below the reactant progress variable. No
        flamelet data lie there, so the values are extrapolated from the fresh mixture, and every
        smoothing neighbour of a node there lies inwards, towards the reacting mixture. The
        smoother then lifts the fresh edge towards the flame rather than denoising it: measured
        on the H2 manifold, fresh reactants read 4 K too hot at their correct enthalpy, so a CFD
        inlet specified by temperature started 5 kJ/kg low in enthalpy. Fitting these nodes
        exactly leaves the rest of the level smoothed as before.
        """
        if not self._is3D():
            return super()._unsmoothedNodes(level_index)
        pv_unb, _ = self.__reactionProgressBounds(self._table_levels[level_index])
        return np.where(self._table_nodes[level_index][:, 0] < pv_unb)[0]

    def _hullEdgeAnchors(self, level_index:int):
        """Anchor the chemically frozen progress-variable edges of a table level.

        A level is meshed over progress variable from just below the reactant value at its
        mixture fraction up to the equilibrium value. Both extremes are frozen: the fresh
        mixture has not started reacting and the equilibrium products have nothing left to
        react. Every tabulated rate is therefore exactly zero on both edges, whatever the
        enthalpy, and the composition on the burnt edge is the equilibrium composition the mesh
        domain itself was built from -- so anchoring it keeps the node's composition and its
        progress-variable coordinate describing the same state.

        The reactant edge sits a little below the reactant progress variable (the mesh domain is
        padded there), so its composition is an extrapolation rather than the fresh mixture and
        is left alone; only its rates are restated.

        :param level_index: table level index.
        :type level_index: int
        :return: {variable name: (node indices, values)}.
        :rtype: dict[str, tuple[np.ndarray[int], np.ndarray[float]]] | None
        """
        hull_ids = self._table_hullnodes[level_index]
        if hull_ids is None or len(hull_ids) == 0:
            return None
        hull_ids = np.asarray(hull_ids, dtype=int)

        pv_nodes = self._table_nodes[level_index][:, 0]
        pv_span = np.ptp(pv_nodes)
        if pv_span <= 0.0:
            return None
        tolerance = self.__hull_edge_tolerance * pv_span

        idx_burnt = hull_ids[pv_nodes[hull_ids] >= pv_nodes.max() - tolerance]
        idx_fresh = hull_ids[pv_nodes[hull_ids] <= pv_nodes.min() + tolerance]
        idx_frozen = np.unique(np.concatenate((idx_burnt, idx_fresh)))
        if len(idx_frozen) == 0:
            return None

        anchors = {}
        for var in self.__frozenStateVars():
            anchors[var] = (idx_frozen, np.zeros(len(idx_frozen)))

        # Equilibrium composition of this level, reproduced exactly as _createPointCloudForTableLevel
        # derives the burnt progress variable that bounds the mesh domain.
        if len(idx_burnt) > 0:
            self._Config.gas.set_mixture_fraction(self._table_levels[level_index], \
                                                  self._Config.GetFuelString(), self._Config.GetOxidizerString())
            self._Config.gas.TP = self._Config.GetUnbTempBounds()[1], DefaultSettings_FGM.pressure
            self._Config.gas.equilibrate("TP")
            Y_burnt = self._Config.gas.Y
            for sp in self._Config.gas.species_names:
                var = "Y-%s" % sp
                if var in self._data_in_table[level_index]:
                    anchors[var] = (idx_burnt, np.full(len(idx_burnt), Y_burnt[self._Config.gas.species_index(sp)]))

        return anchors

    def __frozenStateVars(self):
        """Tabulated quantities that vanish wherever the chemistry is frozen."""
        columns = list(self._data_in_table[0].keys())
        frozen = [FGMVars.ProdRateTot_PV.name, FGMVars.Heat_Release.name]
        frozen += [var for var in columns if "Y_dot" in var]
        return [var for var in frozen if var in columns]

    def setHullEdgeTolerance(self, tolerance:float=1e-3):
        """Width of the band, as a fraction of a level's progress variable range, that counts as
        a hull edge when anchoring.

        :param tolerance: band width as a fraction of the progress variable range, defaults to 1e-3
        :type tolerance: float, optional
        :raises Exception: if the tolerance is not positive.
        """
        if tolerance <= 0:
            raise Exception("Hull edge tolerance should be positive")
        self.__hull_edge_tolerance = tolerance
        return

    def SetSourceDeadBands(self, c_unburnt:float=0.005, c_burnt:float=0.005):
        """Width of the dead bands, in normalized reaction progress, at which the preferential
        diffusion closure sources are switched off at the unburnt and burnt ends of each level.

        The closure sources are extracted from the ends of the 1D flamelets, where the diffusion
        fluxes are being differenced against a vanishing flame structure, so the two extremes of
        each table level are the least trustworthy part of the field. The reference OpenFOAM
        implementation trims them with a hard band on its scaled progress variable; measured on
        its manifold that band discards 0.53% of the total enthalpy closure magnitude and 1.35%
        of the mixture fraction closure, which the defaults here bracket (1.7% and 0.8% on the
        same manifold) while cutting on the level's own normalized progress variable instead of a
        table-specific constant.

        Set both to zero to disable. The trim is applied before the conservation projection, so
        the sources still integrate to zero over each flamelet afterwards -- which is exactly why
        it belongs here rather than in the CFD solver, where trimming the tails of a
        divergence-shaped source would leave a net spurious production term behind.

        :param c_unburnt: dead band at the unburnt end, as a fraction of the progress variable
            range between the level's reactant and equilibrium values, defaults to 0.005. Table
            nodes below the reactant value (the padding of the level's mesh domain) are always in
            the band. A state in the band only reads exactly zero if every vertex of its triangle
            is in the band too, so the band should be wider than the node spacing at its edge.
        :type c_unburnt: float, optional
        :param c_burnt: dead band at the burnt end, defaults to 0.005
        :type c_burnt: float, optional
        :raises Exception: if the bands are negative or together span the whole level.
        """
        if c_unburnt < 0.0 or c_burnt < 0.0:
            raise Exception("Source dead bands should not be negative.")
        if c_unburnt + c_burnt >= 1.0:
            raise Exception("Source dead bands should not span the entire progress variable range.")
        self.__deadband_unburnt = c_unburnt
        self.__deadband_burnt = c_burnt
        return

    def ApplyClosureSourceDeadBands(self):
        """Zero the preferential diffusion closure sources within the configured dead bands.

        The normalized progress variable runs from the level's reactant value to its equilibrium
        value, so the trim follows the manifold as both limits move with mixture fraction. It is
        deliberately NOT measured against the meshed extent: the mesh domain is padded below the
        reactant value, and a band measured from the lowest node then lies entirely inside that
        padding and never reaches the fresh reactants it is meant to trim.
        """
        self.__deadband_nodes = None
        if self.__deadband_unburnt <= 0.0 and self.__deadband_burnt <= 0.0:
            return
        source_vars = self.__closureSourceVars()
        if len(source_vars) == 0:
            return

        n_zeroed = n_total = 0
        self.__deadband_nodes = [np.zeros(len(nodes), dtype=bool) for nodes in self._table_nodes]
        for iLevel in range(self._N_table_levels):
            pv_nodes = self._table_nodes[iLevel][:, 0]
            pv_unb, pv_b = self.__reactionProgressBounds(self._table_levels[iLevel])
            if pv_b <= pv_unb:
                continue
            pv_normalized = (pv_nodes - pv_unb) / (pv_b - pv_unb)
            in_band = (pv_normalized < self.__deadband_unburnt) | \
                      (pv_normalized > 1.0 - self.__deadband_burnt)
            for var in source_vars:
                # np.where rather than in-place masking: the frame hands back a read-only view
                # under copy-on-write, and assigning into it raises.
                self._data_in_table[iLevel][var] = np.where(
                    in_band, 0.0, self._data_in_table[iLevel][var].values)
            self.__deadband_nodes[iLevel] = in_band
            n_zeroed += int(np.sum(in_band))
            n_total += len(in_band)

        print("Closure source dead bands (unburnt=%.3f, burnt=%.3f): %.1f%% of nodes zeroed."
              % (self.__deadband_unburnt, self.__deadband_burnt, 100.0 * n_zeroed / max(n_total, 1)))
        return

    def __reactionProgressBounds(self, levelValue:float):
        """Reactant and equilibrium progress variable of a table level, derived exactly as
        _createPointCloudForTableLevel derives the bounds of the level's mesh domain."""
        gas = self._Config.gas
        gas.set_mixture_fraction(levelValue, self._Config.GetFuelString(), self._Config.GetOxidizerString())
        gas.TP = self._Config.GetUnbTempBounds()[0], DefaultSettings_FGM.pressure
        pv_unb = self._Config.ComputeProgressVariable(variables=None, flamelet_data=None, Y_flamelet=gas.Y[:, np.newaxis])[0]
        gas.TP = self._Config.GetUnbTempBounds()[1], DefaultSettings_FGM.pressure
        gas.equilibrate("TP")
        pv_b = self._Config.ComputeProgressVariable(variables=None, flamelet_data=None, Y_flamelet=gas.Y[:, np.newaxis])[0]
        return pv_unb, pv_b

    def EnforceClosureSourceConservation(self, rel_tolerance:float=1e-3):
        """Audit and restore flamelet-wise conservation of the tabulated PD closure sources.

        The closure sources are divergences of diffusion fluxes along 1D flamelets: their
        integral over every flamelet vanishes in the flamelet data, but the table
        reconstruction (nodal sampling, linear interpolation within and between levels,
        smoothing, lookup misses) does not preserve that. The reconstructed field is integrated
        along every flamelet trajectory, replicating the SU2 lookup; the nodal values are then
        corrected by a magnitude-weighted minimum-norm projection so that each trajectory
        integral matches the flamelet's own. Raises if the audit still fails afterwards.

        :param rel_tolerance: allowed |defect| / (rho_u*u_u*range(phi)) after projection.
        :type rel_tolerance: float, optional
        """
        from Manifold_Generation.LUT.ClosureSourceConservation import (
            ClosureSourceConservation, LoadTrajectories)

        source_vars = self.__closureSourceVars()
        if len(source_vars) == 0:
            return
        # The trajectory operator replicates SU2's lookup for levels in mixture fraction and a
        # (progress variable, total enthalpy) plane; the table axes follow the configured order of
        # the controlling variables, so any other order would be projected against the wrong axes.
        expected_axes = [DefaultSettings_FGM.name_pv, FGMVars.EnthalpyTot.name, FGMVars.MixtureFraction.name]
        if list(self._Config.GetControllingVariables()) != expected_axes:
            print("WARNING: closure source conservation NOT enforced: it requires the controlling variables "
                  "in the order %s, but the table is spanned by %s."
                  % (expected_axes, list(self._Config.GetControllingVariables())))
            return
        flamelet_dir = self._Config.GetOutputDir()
        if not (os.path.isdir(os.sep.join((flamelet_dir, "freeflame_data"))) or
                os.path.isdir(os.sep.join((flamelet_dir, "burnerflame_data")))):
            print("No flamelet data found for the closure source conservation audit -- skipped.")
            return

        print("Auditing flamelet-wise conservation of the tabulated closure sources...")
        mix_header = "mixfrac_" if self._Config.GetMixtureStatus() else "phi_"
        trajectories = LoadTrajectories(flamelet_dir,
                                        self._Config.GetProgressVariableSpecies(),
                                        self._Config.GetProgressVariableWeights(),
                                        mix_status_header=mix_header)
        level_nodes = [np.asarray(nodes)[:, :2] for nodes in self._table_nodes]
        csc = ClosureSourceConservation(np.asarray(self._table_levels), level_nodes,
                                        self._table_connectivity, self._table_hullnodes)
        trajectories, _ = csc.SelectOnTableTrajectories(trajectories)
        A = csc.BuildOperator(trajectories)
        A, trajectories, _ = csc.DropEmptyRows(A, trajectories)
        print("Constraining on %i flamelets overlapping the tabulated manifold." % len(trajectories))

        for var in source_vars:
            theta0 = np.concatenate([lvl[var].values for lvl in self._data_in_table])
            targets = csc.TrajectoryTargets(trajectories, var)
            audit0 = csc.Audit(A, theta0, trajectories, var, targets)
            print("%s before projection: median |defect|/scale = %.3e, max = %.3e"
                  % (var, np.median(np.abs(audit0)), np.max(np.abs(audit0))))
            if np.max(np.abs(theta0)) == 0.0:
                continue
            # The dead-band nodes stay at zero: without this the correction is spread onto them
            # too and re-introduces a source in the fresh reactants the bands were meant to clear.
            frozen = np.concatenate(self.__deadband_nodes) if self.__deadband_nodes is not None else None
            theta = csc.Project(A, theta0, trajectories, var, targets, frozen=frozen)
            audit1 = np.abs(csc.Audit(A, theta, trajectories, var, targets))
            median_defect, max_defect = np.median(audit1), np.max(audit1)
            print("%s after projection : median |defect|/scale = %.3e, p95 = %.3e, max = %.3e"
                  % (var, median_defect, np.percentile(audit1, 95), max_defect))

            # A defect surviving on individual flamelets is a resolution property: the table
            # simply cannot represent the source profile there, and the remedy is more nodes,
            # not a different correction. A defect surviving in the MEDIAN means the whole
            # table is non-conservative, which is the failure this audit exists to catch.
            if median_defect > rel_tolerance:
                raise Exception("Closure source conservation audit FAILED for %s: median "
                                "|defect|/scale = %.3e after projection (tolerance %.1e). The "
                                "table is systematically non-conservative in this variable."
                                % (var, median_defect, rel_tolerance))
            if max_defect > rel_tolerance:
                n_over = int(np.sum(audit1 > rel_tolerance))
                print("WARNING: %s exceeds the %.1e tolerance on %i of %i flamelets (max %.3e). "
                      "These are states the table resolution cannot represent conservatively; "
                      "increase the node count per level if the CFD shows drift."
                      % (var, rel_tolerance, n_over, len(audit1), max_defect))

            for iLevel, vals in enumerate(csc.SplitToLevels(theta)):
                self._data_in_table[iLevel][var] = vals
        return

    def _postProcessTableData(self):
        """Apply the preferential diffusion closure corrections to the finished table data.

        Runs after smoothing deliberately: the smoothing pass rewrites every non-controlling
        variable, so applying the conservation projection before it would have the correction
        smoothed away again.
        """
        self.ApplyClosureSourceDeadBands()
        self.EnforceClosureSourceConservation()
        return

    def _writeAdditionalInfoToTable(self, fid):
        fid.write("Fuel:\n")
        fid.write(",".join(["%.2f:%s" % (w, sp) for w, sp in zip(self._Config.GetFuelWeights(), self._Config.GetFuelDefinition())]))
        fid.write("\nOxidizer:\n")
        fid.write(",".join(["%.2f:%s" % (w, sp) for w, sp in zip(self._Config.GetOxidizerWeights(), self._Config.GetOxidizerDefinition())]))
        fid.write("\nProgress variable:\n")
        fid.write("+".join(["(%+.3e*%s)" % (w, sp) for w, sp in zip(self._Config.GetProgressVariableWeights(), self._Config.GetProgressVariableSpecies())]))
        fid.write("\n\n")
        return