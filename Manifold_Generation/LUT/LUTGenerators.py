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
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
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

    __flamelet_data:pd.DataFrame = None
    __envelope_cache:dict = None

    __pointcloud_resolution:int = DefaultSettings_FGM.table_pointcloud_resolution
    __N_boundary_bins:int = DefaultSettings_FGM.table_boundary_bins
    __N_transverse_bins:int = DefaultSettings_FGM.table_boundary_transverse_bins

    # A closed perimiter needs a lower and an upper envelope station on either side of at least one interior station.
    __min_populated_bins:int = 3

    def __init__(self, config_in:Config_FGM):
        super().__init__(config_in)
        self.__envelope_cache = {}
        self._getFluidDataForInterpolator()
        return

    def setPointCloudResolution(self, resolution:int=DefaultSettings_FGM.table_pointcloud_resolution):
        """Specify the number of nodes per direction of the uniform reference lattice spanning each table level.
        The lattice provides the reference point cloud from which gradient normalization factors are derived.

        :param resolution: nodes per direction, defaults to DefaultSettings_FGM.table_pointcloud_resolution
        :type resolution: int, optional
        :raises Exception: if fewer than two nodes per direction are specified.
        """
        if resolution < 2:
            raise Exception("Point cloud resolution should be at least two nodes per direction.")
        self.__pointcloud_resolution = resolution
        return

    def setBoundaryBinCount(self, N_bins:int=DefaultSettings_FGM.table_boundary_bins):
        """Specify the number of bins along the progress variable used to extract the data envelope bounding
        each table level. Higher values follow the flamelet data more closely, at the cost of a perimiter that
        is more sensitive to the discrete spacing between flamelets.

        :param N_bins: number of bins, defaults to DefaultSettings_FGM.table_boundary_bins
        :type N_bins: int, optional
        :raises Exception: if fewer bins than the minimum required to close a perimiter are specified.
        """
        if N_bins < self.__min_populated_bins:
            raise Exception("Number of boundary bins should be at least %i." % self.__min_populated_bins)
        self.__N_boundary_bins = N_bins
        return

    def setBoundaryTransverseBinCount(self, N_bins:int=DefaultSettings_FGM.table_boundary_transverse_bins):
        """Specify the number of bins along the enthalpy used to extract the equilibrium boundary at maximum
        progress variable.

        That boundary is the locus of the equilibrium composition, which shifts with enthalpy. It runs almost
        parallel to the enthalpy axis, so a sweep along the progress variable places it entirely within its
        last few bins and cannot resolve it; sweeping along the enthalpy instead resolves it directly. Bins
        along the enthalpy run parallel to the flamelets, so this count should stay low enough for every bin to
        collect several flamelets.

        :param N_bins: number of bins, defaults to DefaultSettings_FGM.table_boundary_transverse_bins
        :type N_bins: int, optional
        :raises Exception: if fewer bins than the minimum required to close a perimiter are specified.
        """
        if N_bins < self.__min_populated_bins:
            raise Exception("Number of transverse boundary bins should be at least %i." % self.__min_populated_bins)
        self.__N_transverse_bins = N_bins
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
    
    def ClampSourceTerms(self, species_list:list[str], pv_frac:float=0.99, abs_tol:float=1e-3):
        """Clamp source terms of selected species to zero near the burnt (high-PV) boundary.

        For each table level the per-level progress variable maximum is computed. At every node
        where PV >= ``pv_frac * PV_max`` the net, positive, and negative source terms
        (``Y_dot_net-<sp>``, ``Y_dot_pos-<sp>``, ``Y_dot_neg-<sp>``) of each species in
        ``species_list`` are set to zero when their absolute value is below ``abs_tol``.
        In addition, the production terms are clipped to be non-negative and the consumption
        terms non-positive over all nodes, as SU2 evaluates the species source as
        ``source_prod + source_cons * y_aux``.

        Call this method after :meth:`generateTable` and before :meth:`writeSU2Table`.

        :param species_list: species names to clamp (e.g. ``['CO', 'H2']``).
        :type species_list: list[str]
        :param pv_frac: fraction of the PV maximum above which clamping is applied, defaults to 0.99
        :type pv_frac: float, optional
        :param abs_tol: source terms with an absolute value below this are set to zero, defaults to 1e-3
        :type abs_tol: float, optional
        :raises Exception: if the table data have not been generated yet.
        """
        if not species_list:
            return
        if len(self._data_in_table) == 0:
            raise Exception("Table data are not available, run generateTable before clamping source terms.")

        name_pv = DefaultSettings_FGM.name_pv
        if name_pv not in self._data_in_table[0]:
            print("ClampSourceTerms: %s is not part of the table data, skipping." % name_pv)
            return

        # Collect the source term names per role, skipping those absent from the table.
        vars_net, vars_pos, vars_neg = [], [], []
        for species in species_list:
            for name_format, bucket in (("Y_dot_net-%s", vars_net),
                                        ("Y_dot_pos-%s", vars_pos),
                                        ("Y_dot_neg-%s", vars_neg)):
                var_name = name_format % species
                if var_name in self._table_vars:
                    bucket.append(var_name)
                else:
                    print("ClampSourceTerms: variable %s is not a table variable, skipping." % var_name)

        vars_to_clamp = vars_net + vars_pos + vars_neg
        if len(vars_to_clamp) == 0:
            return

        # When the progress variable sweeps across table levels it is constant within a level, so the
        # burnt boundary follows from the progress variable maximum over all levels instead.
        pv_is_level_cv = self._is3D() and (self._table_cv_names[-1] == name_pv)
        if pv_is_level_cv:
            pv_max_global = max(np.max(data_level[name_pv].to_numpy(dtype=float)) for data_level in self._data_in_table)

        N_nodes_clamped = 0
        for iLevel in range(self._N_table_levels):
            data_level = self._data_in_table[iLevel]
            pv_nodes = data_level[name_pv].to_numpy(dtype=float)
            pv_max = pv_max_global if pv_is_level_cv else np.max(pv_nodes)
            nodes_near_products = pv_nodes >= pv_frac * pv_max
            N_nodes_clamped += int(np.sum(nodes_near_products))

            # Near the burnt boundary: zero the source terms of negligible magnitude.
            for var_name in vars_to_clamp:
                source_terms = data_level[var_name].to_numpy(dtype=float, copy=True)
                source_terms[nodes_near_products & (np.abs(source_terms) < abs_tol)] = 0.0
                data_level[var_name] = source_terms

            # Sign constraints over all nodes: production is never negative, consumption never positive.
            for var_name in vars_pos:
                data_level[var_name] = np.clip(data_level[var_name].to_numpy(dtype=float), 0.0, None)
            for var_name in vars_neg:
                data_level[var_name] = np.clip(data_level[var_name].to_numpy(dtype=float), None, 0.0)

        print("ClampSourceTerms: clamped %i nodes across %i table levels (pv_frac=%.3f, abs_tol=%.2e)." \
              % (N_nodes_clamped, self._N_table_levels, pv_frac, abs_tol))
        return

    def visualizeTableLevelPerimiter(self, level_index:int=0, detail_fraction:float=0.05, \
                                     save_path:str=None, show:bool=False):
        """Visualize the mesh of a table level and its perimiter against the flamelet data the perimiter
        was extracted from. Data outside the perimiter is data the table discards, while a perimiter
        enclosing no data marks a region the table extrapolates into.

        The figure pairs the full table level with a detail at mid-domain, where the individual cells,
        perimiter edges and data points are separable; at table resolution they overlap into solid colour.

        :param level_index: index of the table level to visualize, defaults to 0
        :type level_index: int, optional
        :param detail_fraction: width of the detail view as a fraction of the table level, defaults to 0.05
        :type detail_fraction: float, optional
        :param save_path: file path to save the figure to, defaults to None
        :type save_path: str, optional
        :param show: keep the figure open for display, defaults to False
        :type show: bool, optional
        :raises Exception: if the table has not been generated yet.
        :raises Exception: if the table level index is out of range.
        """
        if len(self._data_in_table) == 0:
            raise Exception("Table data are not available, run generateTable before visualizing the perimiter.")
        self._checkTableLevelIndex(level_index)

        level_value = self._table_levels[level_index]
        level_nodes = self._table_nodes[level_index]
        level_cells = self._table_connectivity[level_index]
        level_data = self.getFlameletDataForTableLevel(level_value)

        # An edge shared by a single cell lies on the perimiter of the mesh. The perimiter node indices are
        # not ordered along the perimiter, so the edges are collected from the connectivity instead.
        edges = np.sort(np.vstack((level_cells[:, [0, 1]], \
                                   level_cells[:, [1, 2]], \
                                   level_cells[:, [2, 0]])), axis=1)
        unique_edges, cells_per_edge = np.unique(edges, axis=0, return_counts=True)
        perimiter_segments = level_nodes[unique_edges[cells_per_edge == 1]][:, :, :2]

        x_span = level_nodes[:, 0].max() - level_nodes[:, 0].min()
        y_span = level_nodes[:, 1].max() - level_nodes[:, 1].min()
        x_detail = 0.5*(level_nodes[:, 0].min() + level_nodes[:, 0].max())
        near_middle = np.absolute(perimiter_segments[:, :, 0].mean(axis=1) - x_detail) < detail_fraction*x_span
        y_detail = perimiter_segments[near_middle][:, :, 1].min()

        fig, axes = plt.subplots(1, 2, figsize=(17, 7), constrained_layout=True)
        for ax, is_detail in zip(axes, (False, True)):
            ax.triplot(level_nodes[:, 0], level_nodes[:, 1], level_cells, \
                       lw=0.3 if is_detail else 0.15, color='0.55', zorder=2)
            if level_data is not None:
                ax.plot(level_data[:, 0], level_data[:, 1], '.', color='red', zorder=1, \
                        ms=4.0 if is_detail else 0.6, \
                        label="flamelet data (%i points)" % len(level_data))
            ax.add_collection(LineCollection(perimiter_segments, colors='k', \
                                             linewidths=2.0 if is_detail else 1.4, zorder=3))
            ax.plot([], [], '-', color='k', lw=1.6, label="mesh perimiter")
            ax.plot([], [], '-', color='0.55', lw=1.0, label="mesh (%i nodes)" % len(level_nodes))
            if is_detail:
                ax.set_xlim(x_detail - 0.8*detail_fraction*x_span, x_detail + 0.8*detail_fraction*x_span)
                ax.set_ylim(y_detail - 0.02*y_span, y_detail + 0.06*y_span)
                ax.set_title("detail of the lower perimiter", fontsize=12)
            else:
                ax.set_title("full table level", fontsize=12)
                ax.legend(markerscale=10, fontsize=10, loc='upper right')
            ax.set_xlabel(self._table_cv_names[0], fontsize=13)
            ax.set_ylabel(self._table_cv_names[1], fontsize=13)

        if self._is3D():
            fig.suptitle("Table level %i (%s = %.5f) and the flamelet data its perimiter was built from" \
                         % (level_index, self._table_cv_names[2], level_value), fontsize=14)
        else:
            fig.suptitle("Mesh perimiter and the flamelet data it was built from", fontsize=14)

        self._saveOrCloseFigure(fig, save_path, show)
        return

    def _getFluidDataForInterpolator(self):
        flamelet_data_filename = os.sep.join((self._Config.GetOutputDir(), self._Config.GetConcatenationFileHeader()+"_full.csv"))
        flameletDataPointCloud = pd.read_csv(flamelet_data_filename)
        self._state_quantities = list(flameletDataPointCloud.keys())
        self._table_vars = self._state_quantities.copy()
        self.__flamelet_data = flameletDataPointCloud
        return flameletDataPointCloud
    
    def _processTableLevels(self):
        super()._processTableLevels()

        if not self.__spansReactantsToProducts():
            # A plane not spanned by the progress variable has no reactant and product state to bound it, so
            # its extent follows the flamelet data. That data describes one manifold, not a sweep of levels.
            if self._is3D():
                raise Exception("Table levels require a table plane spanned by %s." % DefaultSettings_FGM.name_pv)
            return

        if self._is2D():
            mixture_status = self._Config.GetMixtureBounds()[0]
            if self._Config.DefineMixtureStatus():
                mixture_fraction = mixture_status
            else:
                self._Config.gas.set_equivalence_ratio(mixture_status, self._Config.GetFuelString(),self._Config.GetOxidizerString())
                mixture_fraction = self._Config.gas.mixture_fraction(self._Config.GetFuelString(),self._Config.GetOxidizerString())
            self._table_levels[0] = mixture_fraction
        return
    
    def __spansReactantsToProducts(self):
        """Whether the table plane runs from the reactant to the product state, which is the case when the
        progress variable spans it. Such a plane is bounded by states that are evaluated directly; a plane
        spanned by any other controlling variable is bounded by the flamelet data alone.
        """
        return self._table_cv_names[0] == DefaultSettings_FGM.name_pv

    def __reactantProductBounds(self, levelValue:float):
        """Evaluate the reactant and product states bounding a table level.

        :param levelValue: mixture fraction of the table level.
        :type levelValue: float
        :return: reactant and product progress variable, minimum and maximum enthalpy, and the reactant
            enthalpy at the minimum reactant temperature.
        :rtype: tuple[float]
        """
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

        return pv_unb, pv_b, h_min, h_max, h_min_unb

    def __scalePlanarCoords(self, cv_pointcloud:np.ndarray[float]):
        """Scale controlling variable data to the normalized space in which the table is discretized.
        """
        if self._is2D():
            return self._scaler_controlling_variables.transform(cv_pointcloud[:,:2])
        return self._scaler_controlling_variables.transform(cv_pointcloud)

    def getFlameletDataForTableLevel(self, levelValue:float):
        """Retrieve the flamelet data a table level is constructed from, in controlling variable values. For a
        three-dimensional table this is the band of flamelets around the level rather than the whole manifold,
        so this is the data that shaped the perimiter of that particular level.

        :param levelValue: value of the table level.
        :type levelValue: float
        :return: controlling variable values of the contributing flamelet data, or None if no data is
            available.
        :rtype: np.ndarray[float]
        """
        if self.__flamelet_data is None:
            return None

        flamelet_data = self.__flamelet_data
        if self._is3D() and len(self._table_levels) > 1:
            level_spacing = np.mean(np.diff(np.sort(self._table_levels)))
            if level_spacing > 0:
                bandwidth = DefaultSettings_FGM.table_boundary_level_bandwidth * level_spacing
                level_values = flamelet_data[self._table_cv_names[2]].to_numpy()
                flamelet_data = flamelet_data[np.absolute(level_values - levelValue) <= bandwidth]

        if len(flamelet_data) < self.__min_populated_bins:
            return None

        return np.column_stack(tuple(flamelet_data[cv].to_numpy() for cv in self._table_cv_names))

    def __flameletDataForTableLevel(self, levelValue:float):
        """Scaled plane coordinates of the flamelet data contributing to a table level.
        """
        cv_data = self.getFlameletDataForTableLevel(levelValue)
        if cv_data is None:
            return None
        return self.__scalePlanarCoords(cv_data)[:, :2]

    def __envelopeForTableLevel(self, levelValue:float):
        """Bound the flamelet data of a table level by three branches.

        The flamelet data of a table level is a family of curves that are discretely spaced in enthalpy. Tracing
        a hull along the data points would therefore cut into the gaps between neighbouring flamelets. Binning
        the data along the progress variable instead is transverse to those curves, so the resulting envelope
        follows the extremes of the data without resolving the spacing between individual flamelets. That sweep
        yields the lower and upper enthalpy branches, which bound the burner-stabilized limit and the maximum
        enthalpy.

        The boundary at maximum progress variable is the locus of the equilibrium composition, which shifts with
        enthalpy and runs almost parallel to the enthalpy axis; a progress variable sweep confines it to its last
        few bins, so it is swept along the enthalpy instead.

        :param levelValue: value of the table level.
        :type levelValue: float
        :return: the stations of the lower, upper and equilibrium branch, or None if insufficient flamelet data
            is available. Branches are simplified independently and therefore carry their own stations.
        :rtype: dict
        """
        if levelValue in self.__envelope_cache:
            return self.__envelope_cache[levelValue]

        # Cache the absence of an envelope as well, so a level without flamelet data is evaluated only once.
        self.__envelope_cache[levelValue] = None

        planar_data = self.__flameletDataForTableLevel(levelValue)
        if planar_data is None:
            return None


        pv_data = planar_data[:, 0]
        h_data = planar_data[:, 1]

        along_pv = self.__sweepEnvelope(pv_data, h_data, self.__N_boundary_bins)
        if along_pv is None:
            return None
        pv_stations, h_lower, h_upper = along_pv

        h_lower = self.__smoothEnvelope(h_lower, take_minimum=True)
        h_upper = self.__smoothEnvelope(h_upper, take_minimum=False)

        stations_lower = self.__simplifyEnvelope(pv_stations, h_lower, encloses_from_below=True)
        stations_upper = self.__simplifyEnvelope(pv_stations, h_upper, encloses_from_below=False)

        envelope = {"lower": np.column_stack((pv_stations[stations_lower], h_lower[stations_lower])),
                    "upper": np.column_stack((pv_stations[stations_upper], h_upper[stations_upper])),
                    "equilibrium": None}

        along_h = self.__sweepEnvelope(h_data, pv_data, self.__N_transverse_bins)
        if along_h is not None:
            h_stations, _, pv_equilibrium = along_h
            pv_equilibrium = self.__smoothEnvelope(pv_equilibrium, take_minimum=False)
            stations_equilibrium = self.__simplifyEnvelope(h_stations, pv_equilibrium, encloses_from_below=False)
            envelope["equilibrium"] = np.column_stack((pv_equilibrium[stations_equilibrium],
                                                       h_stations[stations_equilibrium]))

        self.__envelope_cache[levelValue] = envelope
        return envelope

    def __sweepEnvelope(self, coordinate_along:np.ndarray[float], coordinate_across:np.ndarray[float], N_bins:int):
        """Bin the data along one coordinate and collect the extremes of the other in each bin.

        :return: the bin stations with the lower and upper extreme in each, or None if too few bins are
            populated to describe a branch.
        :rtype: tuple[np.ndarray[float]]
        """
        station_min, station_max = np.min(coordinate_along), np.max(coordinate_along)
        if station_max <= station_min:
            return None

        bin_edges = np.linspace(station_min, station_max, N_bins + 1)
        bin_index = np.clip(np.digitize(coordinate_along, bin_edges) - 1, 0, N_bins - 1)

        branch_lower = np.full(N_bins, np.nan)
        branch_upper = np.full(N_bins, np.nan)
        for iBin in range(N_bins):
            data_in_bin = (bin_index == iBin)
            # Bins collecting only the tail of a single flamelet are not representative of the boundary.
            if np.count_nonzero(data_in_bin) >= DefaultSettings_FGM.table_boundary_min_points_per_bin:
                branch_lower[iBin] = np.min(coordinate_across[data_in_bin])
                branch_upper[iBin] = np.max(coordinate_across[data_in_bin])

        populated_bins = np.invert(np.isnan(branch_lower))
        if np.count_nonzero(populated_bins) < self.__min_populated_bins:
            return None

        stations = (0.5*(bin_edges[:-1] + bin_edges[1:]))[populated_bins]

        # Anchor the outermost stations on the extremes of the data rather than on the centre of their bin.
        stations[0] = station_min
        stations[-1] = station_max

        return stations, branch_lower[populated_bins], branch_upper[populated_bins]

    def __polylineIntersection(self, polyline_a:np.ndarray[float], polyline_b:np.ndarray[float]):
        """Locate the last point at which two polylines cross.

        :return: the crossing point with the index of the crossed segment in each polyline, or None.
        :rtype: tuple
        """
        crossing = None
        for iSegment in range(len(polyline_a) - 1):
            start_a = polyline_a[iSegment]
            direction_a = polyline_a[iSegment + 1] - start_a
            for jSegment in range(len(polyline_b) - 1):
                start_b = polyline_b[jSegment]
                direction_b = polyline_b[jSegment + 1] - start_b

                determinant = direction_a[0]*direction_b[1] - direction_a[1]*direction_b[0]
                if abs(determinant) < np.finfo(float).eps:
                    continue

                offset = start_b - start_a
                along_a = (offset[0]*direction_b[1] - offset[1]*direction_b[0])/determinant
                along_b = (offset[0]*direction_a[1] - offset[1]*direction_a[0])/determinant
                if 0 <= along_a <= 1 and 0 <= along_b <= 1:
                    crossing = (start_a + along_a*direction_a, iSegment, jSegment)
        return crossing

    def __simplifyEnvelope(self, pv_envelope:np.ndarray[float], h_envelope:np.ndarray[float], encloses_from_below:bool):
        """Reduce an envelope branch to the stations needed to describe it within the simplification tolerance.

        A segment spanning several stations is accepted while it stays on the enclosing side of each station it
        spans and deviates from none of them by more than the tolerance. The first condition preserves the
        flamelet data enclosed by the envelope, the second bounds how far the perimiter departs from it.

        :param pv_envelope: progress variable of the envelope stations.
        :type pv_envelope: np.ndarray[float]
        :param h_envelope: enthalpy of the envelope stations.
        :type h_envelope: np.ndarray[float]
        :param encloses_from_below: whether the branch bounds the data from below.
        :type encloses_from_below: bool
        :return: indices of the retained stations.
        :rtype: np.ndarray[int]
        """
        tolerance = DefaultSettings_FGM.table_boundary_simplification_tolerance * self._base_cell_size
        N_stations = len(pv_envelope)
        if tolerance <= 0 or N_stations < self.__min_populated_bins:
            return np.arange(N_stations)

        retained_stations = [0]
        anchor = 0
        while anchor < N_stations - 1:
            furthest_station = anchor + 1
            for candidate in range(anchor + 2, N_stations):
                spanned = slice(anchor + 1, candidate)
                interpolant = (pv_envelope[spanned] - pv_envelope[anchor])/(pv_envelope[candidate] - pv_envelope[anchor])
                h_segment = h_envelope[anchor] + interpolant*(h_envelope[candidate] - h_envelope[anchor])

                deviation = h_segment - h_envelope[spanned]
                encloses = np.all(deviation <= 0) if encloses_from_below else np.all(deviation >= 0)
                if encloses and np.max(np.absolute(deviation)) <= tolerance:
                    furthest_station = candidate
                else:
                    break
            retained_stations.append(furthest_station)
            anchor = furthest_station

        return np.array(retained_stations)

    def __smoothEnvelope(self, envelope:np.ndarray[float], take_minimum:bool):
        """Smooth an envelope with a running extremum, which displaces the envelope outward only and can
        therefore never exclude flamelet data from the table.
        """
        if DefaultSettings_FGM.table_boundary_smoothing_window <= 1:
            return envelope

        half_window = DefaultSettings_FGM.table_boundary_smoothing_window // 2
        padded_envelope = np.pad(envelope, half_window, mode="edge")
        smoothed_envelope = np.empty(len(envelope))
        for iStation in range(len(envelope)):
            window = padded_envelope[iStation:iStation + 2*half_window + 1]
            smoothed_envelope[iStation] = np.min(window) if take_minimum else np.max(window)
        return smoothed_envelope

    def _createBoundaryPolylineForTableLevel(self, levelValue:float):
        envelope = self.__envelopeForTableLevel(levelValue)
        if envelope is None:
            if self._verbosity > 0:
                print("Insufficient flamelet data to construct a data envelope at table level %.4e, "
                      "extracting the perimiter as a hull around the reference point cloud instead." % levelValue)
            return None

        lower_branch = envelope["lower"]
        upper_branch = envelope["upper"]
        equilibrium_branch = envelope["equilibrium"]

        perimiter = None
        if equilibrium_branch is not None:
            # The equilibrium branch replaces the segment closing the two enthalpy branches at maximum
            # progress variable. Both enthalpy branches are truncated where they cross it, so the corners are
            # the crossings themselves rather than the disagreement between the two sweeps.
            corner_lower = self.__polylineIntersection(equilibrium_branch, lower_branch)
            corner_upper = self.__polylineIntersection(equilibrium_branch, upper_branch)
            if corner_lower is not None and corner_upper is not None:
                crossing_lower, equilibrium_start, station_lower = corner_lower
                crossing_upper, equilibrium_end, station_upper = corner_upper
                perimiter = np.vstack((lower_branch[:station_lower + 1],
                                       [crossing_lower],
                                       equilibrium_branch[equilibrium_start + 1:equilibrium_end + 1],
                                       [crossing_upper],
                                       upper_branch[:station_upper + 1][::-1]))
            elif self._verbosity > 0:
                print("Equilibrium branch does not meet both enthalpy branches at table level %.4e, "
                      "closing the perimiter at maximum progress variable instead." % levelValue)

        if perimiter is None:
            perimiter = np.vstack((lower_branch, upper_branch[::-1]))

        return self.__mergeCoincidentPolylinePoints(perimiter)

    def __mergeCoincidentPolylinePoints(self, polyline:np.ndarray[float]):
        """Remove consecutive perimiter points that are too close together to span a mesh edge.
        """
        merge_tolerance = DefaultSettings_FGM.table_boundary_merge_tolerance * self._base_cell_size

        retained_points = [polyline[0]]
        for point in polyline[1:]:
            if np.linalg.norm(point - retained_points[-1]) > merge_tolerance:
                retained_points.append(point)

        # The closing segment is added by the mesher, so the perimiter should not end where it started.
        if len(retained_points) > self.__min_populated_bins and \
           np.linalg.norm(retained_points[-1] - retained_points[0]) <= merge_tolerance:
            retained_points.pop()

        return np.array(retained_points)

    def __referenceLattice(self, levelValue:float):
        """Uniform lattice spanning a table level, together with the node mask to fall back on when the
        flamelet data yields no envelope.

        A plane spanned by the progress variable runs between the reactant and product states, which are
        evaluated directly, and falls back on the straight line between the reactant and product enthalpy
        limits. Any other plane is spanned by the extent of the flamelet data itself, which leaves nothing to
        fall back on.

        :param levelValue: value of the table level.
        :type levelValue: float
        :return: scaled lattice nodes and the fallback node mask.
        :rtype: tuple[np.ndarray[float], np.ndarray[bool]]
        :raises Exception: if a plane spanned by the flamelet data holds no flamelet data.
        """
        if not self.__spansReactantsToProducts():
            planar_data = self.__flameletDataForTableLevel(levelValue)
            if planar_data is None:
                raise Exception("No flamelet data available to span the table plane.")

            plane_ranges = [np.linspace(np.min(planar_data[:, iPlane]), np.max(planar_data[:, iPlane]), \
                                        self.__pointcloud_resolution) for iPlane in range(2)]
            xgrid, ygrid = np.meshgrid(plane_ranges[0], plane_ranges[1])
            CV_grid_scaled = np.column_stack((xgrid.flatten(), ygrid.flatten()))

            return CV_grid_scaled, np.ones(len(CV_grid_scaled), dtype=bool)

        pv_unb, pv_b, h_min, h_max, h_min_unb = self.__reactantProductBounds(levelValue)

        # Define 2D grid between minimum and maximum progress variable and total enthalpy
        delta_pv_unb = DefaultSettings_FGM.table_reactant_pv_margin*(pv_b - pv_unb)
        pv_range = np.linspace(pv_unb-delta_pv_unb, pv_b, self.__pointcloud_resolution)
        h_range = np.linspace(h_min, h_max, self.__pointcloud_resolution)
        xgrid, ygrid = np.meshgrid(pv_range, h_range)
        zgrid = levelValue*np.ones(np.shape(xgrid))

        CV_grid_init = np.vstack((xgrid.flatten(), ygrid.flatten(), zgrid.flatten())).transpose()

        # Nodes above the straight line between the reactant and product enthalpy limits.
        pv_grid = CV_grid_init[:,0]
        h_grid = CV_grid_init[:,1]
        h_limit = ((h_min_unb - h_min) * pv_grid + (h_min*pv_unb - h_min_unb*pv_b))/(pv_unb - pv_b)

        return self.__scalePlanarCoords(CV_grid_init), h_grid >= h_limit

    def _createPointCloudForTableLevel(self, levelValue:float):
        CV_grid_scaled, nodes_within_fallback_bounds = self.__referenceLattice(levelValue)

        # Retain the nodes that lie within the table level.
        envelope = self.__envelopeForTableLevel(levelValue)
        if envelope is not None:
            idx_keep = self.__nodesWithinEnvelope(CV_grid_scaled[:, :2], envelope)
        else:
            idx_keep = nodes_within_fallback_bounds

        cv_pointcloud_scaled = CV_grid_scaled[idx_keep]
        if self._is2D():
            zcoords = np.zeros([len(cv_pointcloud_scaled), 1])
            pointCloudNodes = np.column_stack((cv_pointcloud_scaled, zcoords))
        else:
            pointCloudNodes = cv_pointcloud_scaled

        if self.__refineGradients:
            self.__computeGradientNormFactors(cv_pointcloud_scaled)

        return pointCloudNodes

    def __nodesWithinEnvelope(self, planar_nodes:np.ndarray[float], envelope:tuple):
        """Locate the nodes enclosed by the lower and upper envelope of a table level.
        """
        lower_branch = envelope["lower"]
        upper_branch = envelope["upper"]
        equilibrium_branch = envelope["equilibrium"]

        pv_nodes = planar_nodes[:, 0]
        h_nodes = planar_nodes[:, 1]

        h_lower_at_node = np.interp(pv_nodes, lower_branch[:, 0], lower_branch[:, 1])
        h_upper_at_node = np.interp(pv_nodes, upper_branch[:, 0], upper_branch[:, 1])

        within_progress_variable_range = (pv_nodes >= lower_branch[0, 0]) & (pv_nodes <= lower_branch[-1, 0])
        within_enthalpy_range = (h_nodes >= h_lower_at_node) & (h_nodes <= h_upper_at_node)

        nodes_within = within_progress_variable_range & within_enthalpy_range
        if equilibrium_branch is not None:
            pv_equilibrium_at_node = np.interp(h_nodes, equilibrium_branch[:, 1], equilibrium_branch[:, 0])
            nodes_within &= (pv_nodes <= pv_equilibrium_at_node)

        return nodes_within

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