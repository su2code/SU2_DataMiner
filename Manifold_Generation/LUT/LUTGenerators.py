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
import glob
from scipy.interpolate import PchipInterpolator
from Common.Interpolators import polylineIntersection
from Common.Properties import EntropicVars,DefaultSettings_FGM
from su2dataminer.generate_data import DataGenerator_CoolProp
from Common.DataDrivenConfig import Config_NICFD,Config_FGM
from Data_Generation.FlameletSolvers import FlameletSolverDict
from Manifold_Generation.LUT.LUTGenerator_Base import SU2TableGenerator_Base
from Manifold_Generation.LUT.MeshTools import MeshThermodynamicPlane, Mesh2DPlane
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
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
        self._orderTableVarsByAxes()
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

    __clampSourceTerms:bool = False
    __clamp_species:list[str] = []
    __clamp_pv_margin:float = 1e-2
    __clamp_abs_tol:float = 1e-3
    __pv_extent_over_levels:tuple = None

    __flamelet_data:pd.DataFrame = None
    __envelope_cache:dict = None
    __reference_lattice:np.ndarray[float] = None

    __pointcloud_resolution:int = DefaultSettings_FGM.table_pointcloud_resolution
    __N_boundary_bins:int = DefaultSettings_FGM.table_boundary_bins
    __N_transverse_bins:int = DefaultSettings_FGM.table_boundary_transverse_bins
    __min_points_per_bin:int = DefaultSettings_FGM.table_boundary_min_points_per_bin
    __level_bandwidth:float = DefaultSettings_FGM.table_boundary_level_bandwidth

    # Controlling variable data of the flamelets, kept apart per flamelet type.
    __flamelet_data_by_type:dict = None
    __simplification_tolerance:float = DefaultSettings_FGM.table_boundary_simplification_tolerance
    __merge_tolerance:float = DefaultSettings_FGM.table_boundary_merge_tolerance

    # A branch needs a station at either end and at least one in between: fewer than that cannot describe
    # a curve, so a sweep yielding fewer populated bins does not produce a branch, and a bin count below
    # this cannot be configured.
    __min_envelope_stations:int = 3

    # A polygon is closed by three points, which is what the perimiter has to retain.
    __min_perimiter_points:int = 3

    # Below this a table level holds no plane to bound.
    __min_flamelet_data_points:int = 3

    # Inserting a station only changes the two segments meeting at it, so the outstanding violations
    # shrink quickly; this bounds the work when a branch is pathological.
    __max_enclosure_passes:int = 8

    # Rounds of simplifying and re-enclosing a branch. The second round collapses what the first round
    # had to insert; further rounds change little.
    __branch_refinement_rounds:int = 2

    # Rounds of raising or lowering the control points of a fitted edge onto the data.
    __edge_fit_rounds:int = 3

    # Points a fitted edge is sampled onto for the mesher.
    __edge_sample_points:int = 40

    # Enthalpy range around a corner, as a fraction of the plane, that the corner is taken from.
    __corner_data_fraction:float = 0.02

    def __init__(self, config_in:Config_FGM):
        super().__init__(config_in)
        # Bound to the instance rather than shared with every other generator in the process.
        self.__envelope_cache = {}
        self.__reference_lattice = None
        self.__flamelet_data_by_type = None
        self._refine_for_gradients_of = []
        self._gradient_refinement_factors = []
        self._gradient_norm_factors = []
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

    def __invalidateEnvelopes(self):
        """Discard the cached data envelopes, which were shaped by the previous boundary settings."""
        self.__envelope_cache = {}
        return

    def setBoundaryBinCount(self, N_bins:int=DefaultSettings_FGM.table_boundary_bins):
        """Specify the number of bins along the progress variable used to extract the data envelope bounding
        each table level. Higher values follow the flamelet data more closely, at the cost of a perimiter that
        is more sensitive to the discrete spacing between flamelets.

        :param N_bins: number of bins, defaults to DefaultSettings_FGM.table_boundary_bins
        :type N_bins: int, optional
        :raises Exception: if fewer bins than the minimum required to close a perimiter are specified.
        """
        if N_bins < self.__min_envelope_stations:
            raise Exception("Number of boundary bins should be at least %i." % self.__min_envelope_stations)
        self.__N_boundary_bins = N_bins
        self.__invalidateEnvelopes()
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
        if N_bins < self.__min_envelope_stations:
            raise Exception("Number of transverse boundary bins should be at least %i." % self.__min_envelope_stations)
        self.__N_transverse_bins = N_bins
        self.__invalidateEnvelopes()
        return

    def setBoundaryMinPointsPerBin(self, N_points:int=DefaultSettings_FGM.table_boundary_min_points_per_bin):
        """Specify how many flamelet data points a bin should hold to contribute a station to the envelope.
        Bins collecting only the tail of a single flamelet are not representative of the boundary, so raising
        this discards more of them, at the cost of shortening the envelope where the data thins out.

        :param N_points: minimum number of data points per bin, defaults to
            DefaultSettings_FGM.table_boundary_min_points_per_bin
        :type N_points: int, optional
        :raises Exception: if fewer than one point per bin is specified.
        """
        if N_points < 1:
            raise Exception("A bin should hold at least one data point to contribute an envelope station.")
        self.__min_points_per_bin = N_points
        self.__invalidateEnvelopes()
        return

    def setBoundaryLevelBandwidth(self, bandwidth:float=DefaultSettings_FGM.table_boundary_level_bandwidth):
        """Specify the band of flamelet data around a table level that shapes its envelope, in table level
        spacings. Only three-dimensional tables draw on a band; a two-dimensional table is built from the
        whole manifold. Widening the band draws in flamelets further from the level, which smooths the
        envelope but blurs the variation of the boundary across levels.

        :param bandwidth: half-width of the band in table level spacings, defaults to
            DefaultSettings_FGM.table_boundary_level_bandwidth
        :type bandwidth: float, optional
        :raises Exception: if the bandwidth is not positive.
        """
        if bandwidth <= 0:
            raise Exception("The table level bandwidth should be positive.")
        self.__level_bandwidth = bandwidth
        self.__invalidateEnvelopes()
        return

    def setBoundarySimplificationTolerance(self, tolerance:float=DefaultSettings_FGM.table_boundary_simplification_tolerance):
        """Specify how far the perimiter may depart from the envelope stations it replaces, as a fraction of
        the coarse cell size. Simplification collapses the steps left by the discrete spacing between
        flamelets into straight segments; a tolerance of zero retains every station.

        :param tolerance: simplification tolerance as a fraction of the coarse cell size, defaults to
            DefaultSettings_FGM.table_boundary_simplification_tolerance
        :type tolerance: float, optional
        :raises Exception: if the tolerance is negative.
        """
        if tolerance < 0:
            raise Exception("The envelope simplification tolerance should not be negative.")
        self.__simplification_tolerance = tolerance
        self.__invalidateEnvelopes()
        return

    def setBoundaryMergeTolerance(self, tolerance:float=DefaultSettings_FGM.table_boundary_merge_tolerance):
        """Specify the distance below which consecutive perimiter points are merged, as a fraction of the
        coarse cell size. Points closer together than this cannot span a mesh edge.

        :param tolerance: merge tolerance as a fraction of the coarse cell size, defaults to
            DefaultSettings_FGM.table_boundary_merge_tolerance
        :type tolerance: float, optional
        :raises Exception: if the tolerance is negative.
        """
        if tolerance < 0:
            raise Exception("The perimiter merge tolerance should not be negative.")
        self.__merge_tolerance = tolerance
        self.__invalidateEnvelopes()
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

    def setSourceTermClamping(self, species_list:list[str], pv_margin:float=1e-2, abs_tol:float=1e-3):
        """Clamp the source terms of selected species to zero near the reactant and product boundaries.

        SU2 evaluates a species source as ``source_prod + source_cons * y_aux``, so production terms are
        clipped to be non-negative and consumption terms non-positive over the whole table. Within a margin
        of the reactant and of the product progress variable the net, positive and negative source terms
        (``Y_dot_net-<sp>``, ``Y_dot_pos-<sp>``, ``Y_dot_neg-<sp>``) are additionally set to zero where
        their absolute value is below ``abs_tol``, as neither boundary sustains reaction.

        Clamping the reactant side leaves burner-stabilized flamelets untouched, since their inflow
        progress variable exceeds the reactant progress variable.

        The clamp is applied during :meth:`generateTable` and may be configured at any point before it.

        :param species_list: species names to clamp (e.g. ``['CO', 'H2']``).
        :type species_list: list[str]
        :param pv_margin: width of each clamped band, as a fraction of the progress variable span between
            the reactant and product state, defaults to 1e-2
        :type pv_margin: float, optional
        :param abs_tol: source terms with an absolute value below this are set to zero, defaults to 1e-3
        :type abs_tol: float, optional
        :raises Exception: if the margin lies outside the interval [0, 0.5).
        :raises Exception: if the tolerance is negative.
        """
        if (pv_margin < 0) or (pv_margin >= 0.5):
            raise Exception("The progress variable margin should lie between 0 and 0.5.")
        if abs_tol < 0:
            raise Exception("The source term tolerance should not be negative.")
        self.__clampSourceTerms = bool(species_list)
        self.__clamp_species = list(species_list)
        self.__clamp_pv_margin = pv_margin
        self.__clamp_abs_tol = abs_tol
        return

    def _postProcessTableData(self):
        """Apply the configured source term clamp to the generated table data."""
        if self.__clampSourceTerms:
            self.__clampSpeciesSourceTerms()
        return

    def __sourceTermsToClamp(self):
        """Collect the table variables holding the source terms of the clamped species, per role.

        :return: the net, production and consumption source term names present in the table.
        :rtype: tuple[list[str]]
        """
        vars_net, vars_pos, vars_neg = [], [], []
        for species in self.__clamp_species:
            for name_format, bucket in (("Y_dot_net-%s", vars_net),
                                        ("Y_dot_pos-%s", vars_pos),
                                        ("Y_dot_neg-%s", vars_neg)):
                var_name = name_format % species
                if var_name in self._table_vars:
                    bucket.append(var_name)
                elif self._verbosity > 0:
                    print("Source term clamping: %s is not a table variable, skipping." % var_name)
        return vars_net, vars_pos, vars_neg

    def __progressVariableExtent(self, levelValue:float, pv_nodes:np.ndarray[float]):
        """Reactant and product progress variable bounding a table level.

        A plane spanned by the progress variable runs between the reactant and product states, which are
        evaluated directly. Any other plane is bounded by the flamelet data alone, so its extremes stand
        in for those states. When the progress variable sweeps the table levels instead of spanning the
        plane it is constant within a level, so those extremes are taken over all levels.

        :return: reactant and product progress variable.
        :rtype: tuple[float]
        """
        if self.__spansReactantsToProducts():
            pv_unb, pv_b, _, _, _ = self.__reactantProductBounds(levelValue)
            return pv_unb, pv_b
        if self.__pv_extent_over_levels is not None:
            return self.__pv_extent_over_levels
        return np.min(pv_nodes), np.max(pv_nodes)

    def __clampedProgressVariableBands(self, levelValue:float, pv_nodes:np.ndarray[float]):
        """Bound the progress variable bands at the reactant and at the product side of a table level.

        :return: upper bound of the reactant band and lower bound of the product band.
        :rtype: tuple[float]
        """
        pv_unb, pv_b = self.__progressVariableExtent(levelValue, pv_nodes)
        margin = self.__clamp_pv_margin * (pv_b - pv_unb)
        return pv_unb + margin, pv_b - margin

    def __clampSpeciesSourceTerms(self):
        """Zero negligible species source terms at the reactant and product boundaries and enforce the
        sign of the production and consumption terms over the whole table.
        """
        name_pv = DefaultSettings_FGM.name_pv
        if name_pv not in self._table_vars:
            print("Source term clamping: %s is not part of the table data, skipping." % name_pv)
            return

        vars_net, vars_pos, vars_neg = self.__sourceTermsToClamp()
        vars_to_clamp = vars_net + vars_pos + vars_neg
        if len(vars_to_clamp) == 0:
            return

        # The progress variable is constant within a level when it sweeps the levels, so its extent is
        # taken over the whole table instead.
        self.__pv_extent_over_levels = None
        if self._is3D() and (self._table_cv_names[-1] == name_pv):
            pv_over_levels = [data_level[name_pv].to_numpy(dtype=float) for data_level in self._data_in_table]
            self.__pv_extent_over_levels = (min(np.min(pv) for pv in pv_over_levels),
                                            max(np.max(pv) for pv in pv_over_levels))

        abs_tol = self.__clamp_abs_tol
        N_nodes_clamped = 0
        for iLevel in range(self._N_table_levels):
            data_level = self._data_in_table[iLevel]
            pv_nodes = data_level[name_pv].to_numpy(dtype=float)
            pv_reactants, pv_products = self.__clampedProgressVariableBands(self._table_levels[iLevel], pv_nodes)

            nodes_clamped = (pv_nodes <= pv_reactants) | (pv_nodes >= pv_products)
            N_nodes_clamped += int(np.count_nonzero(nodes_clamped))

            for var_name in vars_to_clamp:
                source_terms = data_level[var_name].to_numpy(dtype=float, copy=True)
                source_terms[nodes_clamped & (np.absolute(source_terms) < abs_tol)] = 0.0
                data_level[var_name] = source_terms

            # SU2 evaluates the species source as source_prod + source_cons * y_aux.
            for var_name in vars_pos:
                data_level[var_name] = np.clip(data_level[var_name].to_numpy(dtype=float), 0.0, None)
            for var_name in vars_neg:
                data_level[var_name] = np.clip(data_level[var_name].to_numpy(dtype=float), None, 0.0)

        if self._verbosity > 0:
            print("Source term clamping: clamped %i nodes across %i table levels (pv_margin=%.3e, abs_tol=%.2e)." \
                  % (N_nodes_clamped, self._N_table_levels, self.__clamp_pv_margin, abs_tol))
        return

    def visualizeTableLevelPerimiter(self, level_index:int=0, detail_fraction:float=0.05, \
                                     color_by_flamelet_type:bool=False, save_path:str=None, show:bool=False):
        """Visualize the mesh of a table level and its perimiter against the flamelet data the perimiter
        was extracted from. Data outside the perimiter is data the table discards, while a perimiter
        enclosing no data marks a region the table extrapolates into.

        The figure pairs the full table level with a detail at mid-domain, where the individual cells,
        perimiter edges and data points are separable; at table resolution they overlap into solid colour.

        :param level_index: index of the table level to visualize, defaults to 0
        :type level_index: int, optional
        :param detail_fraction: width of the detail view as a fraction of the table level, defaults to 0.05
        :type detail_fraction: float, optional
        :param color_by_flamelet_type: color the flamelet data by the flamelet type it was produced by
            rather than by whether the table covers it, defaults to False
        :type color_by_flamelet_type: bool, optional
        :param save_path: file path to save the figure to, defaults to None
        :type save_path: str, optional
        :param show: display the figure, defaults to False
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
        data_per_type = self.getFlameletDataByTypeForTableLevel(level_value) if color_by_flamelet_type else {}

        # Whether a flamelet data point is covered by the table follows from the cells themselves rather
        # than from the perimiter, which needs no ordering and accounts for any hole in the mesh.
        covered_by_table = None
        if level_data is not None:
            cell_finder = mtri.Triangulation(level_nodes[:, 0], level_nodes[:, 1], level_cells).get_trifinder()
            covered_by_table = cell_finder(level_data[:, 0], level_data[:, 1]) >= 0
            self._reportTableLevelCoverage(level_index, level_nodes, level_cells, covered_by_table)

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

        flamelet_type_colors = plt.get_cmap("tab10").colors
        fig, axes = plt.subplots(1, 2, figsize=(17, 7), constrained_layout=True)
        for ax, is_detail in zip(axes, (False, True)):
            ax.triplot(level_nodes[:, 0], level_nodes[:, 1], level_cells, \
                       lw=0.3 if is_detail else 0.15, color='0.55', zorder=2)
            if color_by_flamelet_type:
                for color, (flamelet_type, type_data) in zip(flamelet_type_colors, data_per_type.items()):
                    ax.plot(type_data[:, 0], type_data[:, 1], '.', color=color, zorder=1, \
                            ms=4.0 if is_detail else 0.8, \
                            label="%s (%i points)" % (flamelet_type, len(type_data)))
            elif level_data is not None:
                # Discarded data is what the table loses, so it is marked out from the rest.
                ax.plot(level_data[covered_by_table, 0], level_data[covered_by_table, 1], '.', \
                        color='0.72', zorder=1, ms=4.0 if is_detail else 0.6, \
                        label="flamelet data in table (%i points)" % np.count_nonzero(covered_by_table))
                if not np.all(covered_by_table):
                    ax.plot(level_data[~covered_by_table, 0], level_data[~covered_by_table, 1], '.', \
                            color='red', zorder=1, ms=5.0 if is_detail else 2.0, \
                            label="flamelet data outside table (%i points)" % np.count_nonzero(~covered_by_table))
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

    def _reportTableLevelCoverage(self, level_index:int, level_nodes:np.ndarray[float], \
                                  level_cells:np.ndarray[int], covered_by_table:np.ndarray[bool]):
        """Report how much of the flamelet data of a table level ends up in the table.

        :param level_index: index of the table level.
        :type level_index: int
        :param level_nodes: controlling variable coordinates of the table nodes.
        :type level_nodes: np.ndarray[float]
        :param level_cells: table level connectivity.
        :type level_cells: np.ndarray[int]
        :param covered_by_table: whether each flamelet data point falls within a table cell.
        :type covered_by_table: np.ndarray[bool]
        """
        if self._verbosity < 1:
            return

        cell_corners = level_nodes[level_cells][:, :, :2]
        edge_a = cell_corners[:, 1] - cell_corners[:, 0]
        edge_b = cell_corners[:, 2] - cell_corners[:, 0]
        covered_area = 0.5*np.sum(np.absolute(edge_a[:, 0]*edge_b[:, 1] - edge_a[:, 1]*edge_b[:, 0]))

        print("Table level %i: %i nodes, %i cells covering %.5e in (%s, %s)" \
              % (level_index, len(level_nodes), len(level_cells), covered_area, \
                 self._table_cv_names[0], self._table_cv_names[1]))
        print("               %i of %i flamelet data points inside the table, %i outside (%.3f%%)" \
              % (np.count_nonzero(covered_by_table), len(covered_by_table), \
                 np.count_nonzero(~covered_by_table), 100*np.mean(~covered_by_table)))
        return

    def _getFluidDataForInterpolator(self):
        # The manifold is read once and kept: the generator needs it on construction to know which
        # quantities exist, and again when the interpolator is built.
        if self.__flamelet_data is None:
            flamelet_data_filename = os.sep.join((self._Config.GetOutputDir(), self._Config.GetConcatenationFileHeader()+"_full.csv"))
            self.__flamelet_data = pd.read_csv(flamelet_data_filename)
            self._state_quantities = list(self.__flamelet_data.keys())
            self._table_vars = self._state_quantities.copy()
            self._orderTableVarsByAxes()

        return self.__flamelet_data

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
        return (self._table_cv_names[0] == DefaultSettings_FGM.name_pv and
                self._table_cv_names[1] == DefaultSettings_FGM.name_enth)

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

    def __scaleEdge(self, progress_variable:np.ndarray[float], enthalpy:np.ndarray[float], levelValue:float):
        """Scale progress variable and enthalpy of a table plane edge to the discretization space."""
        level = levelValue*np.ones(len(progress_variable))
        return self.__scalePlanarCoords(np.column_stack((progress_variable, enthalpy, level)))[:, :2]

    def __fittedEdge(self, planar_data:np.ndarray[float], h_lower:float, h_upper:float, bounds_from_below:bool, \
                     anchor_lower:float=None, anchor_upper:float=None):
        """Fit a smooth edge to the extreme of the flamelet data over a range of enthalpy.

        The edge is single valued in enthalpy, so it is swept along it. Each bin contributes the
        extreme data point it holds, a shape preserving interpolant is fitted through a subset of
        those, and its control points are then raised or lowered until the edge bounds the data.
        Correcting the control points individually keeps a single distant point from displacing the
        whole edge.

        :param planar_data: table plane coordinates of the flamelet data.
        :type planar_data: np.ndarray[float]
        :param h_lower: enthalpy the edge starts at.
        :type h_lower: float
        :param h_upper: enthalpy the edge ends at.
        :type h_upper: float
        :param bounds_from_below: whether the edge bounds the data from below.
        :type bounds_from_below: bool
        :return: plane coordinates of the edge, ordered by ascending enthalpy, or None.
        :rtype: np.ndarray[float]
        """
        within_range = (planar_data[:, 1] >= h_lower) & (planar_data[:, 1] <= h_upper)
        pv_data, h_data = planar_data[within_range, 0], planar_data[within_range, 1]
        if len(pv_data) < self.__min_flamelet_data_points:
            return None

        take_extreme = np.argmin if bounds_from_below else np.argmax
        N_bins = self.__N_boundary_bins
        bin_edges = np.linspace(h_lower, h_upper, N_bins + 1)
        bin_index = np.clip(np.digitize(h_data, bin_edges) - 1, 0, N_bins - 1)

        station_h, station_pv = [], []
        for iBin in range(N_bins):
            in_bin = np.flatnonzero(bin_index == iBin)
            if len(in_bin) >= self.__min_points_per_bin:
                extreme = in_bin[take_extreme(pv_data[in_bin])]
                station_h.append(h_data[extreme])
                station_pv.append(pv_data[extreme])
        if len(station_h) < self.__min_envelope_stations:
            return None

        station_h, station_pv = np.array(station_h), np.array(station_pv)
        # The ends of the edge reach the extremes of the data rather than the last populated bin.
        station_h[0], station_h[-1] = h_lower, h_upper

        control = np.unique(np.linspace(0, len(station_h) - 1, \
                                        DefaultSettings_FGM.table_boundary_edge_control_points).astype(int))
        control_h, control_pv = station_h[control], station_pv[control]

        # The corners the edge meets are fixed, so they are not moved by the fit.
        if anchor_lower is not None:
            control_pv[0] = anchor_lower
        if anchor_upper is not None:
            control_pv[-1] = anchor_upper

        quantile = 100*DefaultSettings_FGM.table_boundary_data_quantile

        # A control point governs the data nearest to it in enthalpy. The fit moves the control
        # points across the progress variable only, so that grouping holds for every round.
        nearest_control = np.argmin(np.absolute(h_data[:, np.newaxis] - control_h[np.newaxis, :]), axis=1)
        governed_by = [np.flatnonzero(nearest_control == iControl) for iControl in range(len(control_h))]

        for _ in range(self.__edge_fit_rounds):
            edge = PchipInterpolator(control_h, control_pv)
            outside = (edge(h_data) - pv_data) if bounds_from_below else (pv_data - edge(h_data))
            for iControl, governed in enumerate(governed_by):
                if len(governed) == 0:
                    continue
                excess = np.percentile(outside[governed], quantile)
                if excess > 0:
                    control_pv[iControl] += -excess if bounds_from_below else excess
            if anchor_lower is not None:
                control_pv[0] = anchor_lower
            if anchor_upper is not None:
                control_pv[-1] = anchor_upper

        stations = np.linspace(h_lower, h_upper, self.__edge_sample_points)
        return np.column_stack((PchipInterpolator(control_h, control_pv)(stations), stations))

    def __checkPerimiter(self, perimiter:np.ndarray[float]):
        """Verify that a perimiter is traversed anticlockwise and does not cross itself.

        :param perimiter: planar coordinates of the perimiter.
        :type perimiter: np.ndarray[float]
        :raises Exception: if the perimiter crosses itself.
        :return: the perimiter, oriented anticlockwise.
        :rtype: np.ndarray[float]
        """
        x, y = perimiter[:, 0], perimiter[:, 1]
        if 0.5*np.sum(x*np.roll(y, -1) - y*np.roll(x, -1)) < 0:
            perimiter = perimiter[::-1]

        start = perimiter
        direction = np.roll(perimiter, -1, axis=0) - perimiter
        determinant = direction[:, None, 0]*direction[None, :, 1] - direction[:, None, 1]*direction[None, :, 0]
        offset = start[None, :, :] - start[:, None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            along_a = (offset[:, :, 0]*direction[None, :, 1] - offset[:, :, 1]*direction[None, :, 0])/determinant
            along_b = (offset[:, :, 0]*direction[:, None, 1] - offset[:, :, 1]*direction[:, None, 0])/determinant
        crosses = (np.absolute(determinant) >= np.finfo(float).eps) \
                  & (along_a > 1e-9) & (along_a < 1 - 1e-9) & (along_b > 1e-9) & (along_b < 1 - 1e-9)
        if np.any(crosses):
            raise Exception("Table level perimiter crosses itself at %i of its %i segments" \
                            % (np.count_nonzero(crosses)//2, len(perimiter)))
        return perimiter

    def __rightConvexChain(self, planar_data:np.ndarray[float], start:np.ndarray[float], end:np.ndarray[float]):
        """Chain of straight segments from start to end that bounds the data towards higher progress
        variable and turns one way only.

        The equilibrium side of a table plane is convex, so a convex chain follows it without the
        overshoot an interpolant through scattered extremes would introduce.

        :param planar_data: table plane coordinates of the flamelet data.
        :type planar_data: np.ndarray[float]
        :param start: corner the chain starts at.
        :type start: np.ndarray[float]
        :param end: corner the chain ends at.
        :type end: np.ndarray[float]
        :return: plane coordinates of the chain, ordered by ascending enthalpy.
        :rtype: np.ndarray[float]
        """
        within = planar_data[(planar_data[:, 1] >= start[1]) & (planar_data[:, 1] <= end[1])]
        candidates = np.vstack((start, within, end))
        candidates = candidates[np.argsort(candidates[:, 1], kind="stable")]

        chain = []
        for point in candidates:
            while len(chain) >= 2:
                before, last = chain[-2], chain[-1]
                turn = (last[0] - before[0])*(point[1] - before[1]) \
                     - (last[1] - before[1])*(point[0] - before[0])
                if turn <= 0:
                    chain.pop()
                else:
                    break
            chain.append(point)
        return np.array(chain)

    def __parametricBoundary(self, levelValue:float):
        """Perimiter of a table plane that runs from the reactant to the product state.

        The plane is a rough triangle: the reactants sit at a single progress variable, the burnt
        products on the equilibrium side, and the two meet in a point at the coldest state. It is
        closed at the top by a straight line at the enthalpy of the hottest reactant state. The
        reactant line and the corner at which it meets the burner-stabilized edge carry the unburnt
        state, so both are exact.

        :param levelValue: value of the table level.
        :type levelValue: float
        :return: planar coordinates of the perimiter, or None if the level holds too little data.
        :rtype: np.ndarray[float]
        """
        _, _, _, _, h_min_unb = self.__reactantProductBounds(levelValue)
        h_corner = self.__scaleEdge(np.array([0.0]), np.array([h_min_unb]), levelValue)[0, 1]

        planar_data = self.__flameletDataForTableLevel(levelValue)
        if planar_data is None:
            return None

        pv_data, h_data = planar_data[:, 0], planar_data[:, 1]

        # The reactants sit at one progress variable, which the data gives exactly.
        pv_reactants = pv_data.min()

        # The plane is closed at the hottest reactant state, so that the top runs at one enthalpy from
        # the reactant line to the products. Data above it is a sliver that the table leaves out.
        on_reactant_line = pv_data <= pv_reactants + self.__corner_data_fraction*np.ptp(pv_data)
        h_top = h_data[on_reactant_line].max()

        planar_data = planar_data[planar_data[:, 1] <= h_top]
        if len(planar_data) < self.__min_flamelet_data_points:
            return None
        pv_data, h_data = planar_data[:, 0], planar_data[:, 1]
        h_bottom = h_data.min()

        # The coldest state is a single point, where the two sides of the triangle meet.
        coldest = h_data < h_bottom + self.__corner_data_fraction*(h_top - h_bottom)
        apex = np.array([pv_data[coldest].max(), h_bottom])

        topmost = h_data > h_top - self.__corner_data_fraction*(h_top - h_bottom)
        corner_right = np.array([pv_data[topmost].max(), h_top])
        corner_left = np.array([pv_reactants, h_top])

        equilibrium_edge = self.__rightConvexChain(planar_data, apex, corner_right)
        burner_edge = self.__fittedEdge(planar_data, h_bottom, h_corner, bounds_from_below=True, \
                                        anchor_lower=apex[0], anchor_upper=pv_reactants)
        if burner_edge is None:
            return None

        # From the coldest state up the equilibrium side, across the top, down the reactant line,
        # then back down the burner-stabilized side to the same point it started from.
        perimiter = np.vstack((equilibrium_edge, [corner_left], [[pv_reactants, h_corner]], burner_edge[::-1]))
        return self.__checkPerimiter(self.__mergeCoincidentPolylinePoints(perimiter))

    def getManifoldRange(self, varname:str):
        """Extent of a quantity over the flamelet data the manifold is described by.

        Useful for placing the table on the data rather than on values chosen in advance, for instance
        to span the table levels over the mixture fractions that were actually generated.

        :param varname: name of the manifold quantity.
        :type varname: str
        :return: lowest and highest value of the quantity in the flamelet data.
        :rtype: tuple[float]
        :raises Exception: if the quantity is not part of the flamelet data.
        """
        flamelet_data = self._getFluidDataForInterpolator()
        if varname not in flamelet_data:
            raise Exception("%s is not part of the flamelet data (%s)" % (varname, ", ".join(self._state_quantities)))

        values = flamelet_data[varname].to_numpy(dtype=float)
        return float(np.min(values)), float(np.max(values))

    def __levelBandwidth(self):
        """Half width in the level controlling variable of the band of flamelet data a table level is
        built from, or None when the table holds a single level that all data belongs to.
        """
        if not (self._is3D() and len(self._table_levels) > 1):
            return None
        level_spacing = np.mean(np.diff(np.sort(self._table_levels)))
        if level_spacing <= 0:
            return None
        return self.__level_bandwidth * level_spacing

    def __readFlameletDataPerType(self):
        """Read the controlling variables of every flamelet solution, kept apart per flamelet type.

        The concatenated manifold does not record which solver produced a data point, so the flamelet
        solutions themselves are read back. Only the columns the controlling variables follow from are
        read, which keeps the sweep over all solutions to seconds rather than minutes.
        """
        pv_columns = ["Y-%s" % species for species in self._Config.GetProgressVariableSpecies()]
        stored_cvs = [cv for cv in self._table_cv_names if cv != DefaultSettings_FGM.name_pv]
        columns = stored_cvs + pv_columns

        data_per_type = {}
        for flamelet_type in self._Config.getFlameletTypes():
            solver = FlameletSolverDict[flamelet_type](self._Config)
            flamelet_files = sorted(glob.glob(os.sep.join((self._Config.GetOutputDir(), \
                                                           solver.getFlameletFolder(), "*", "*.csv"))))
            cv_data = []
            for flamelet_file in flamelet_files:
                flamelet = pd.read_csv(flamelet_file, usecols=columns)
                pv = self._Config.ComputeProgressVariable(list(flamelet.keys()), flamelet.to_numpy())
                cv_data.append(np.column_stack(tuple(pv if cv == DefaultSettings_FGM.name_pv \
                                                     else flamelet[cv].to_numpy() \
                                                     for cv in self._table_cv_names)))
            if len(cv_data) > 0:
                data_per_type[solver.getPlotLabel()] = np.vstack(cv_data)
        return data_per_type

    def getFlameletDataByTypeForTableLevel(self, levelValue:float):
        """Retrieve the flamelet data a table level is constructed from, split by the flamelet type it
        was produced by. The selection is the one getFlameletDataForTableLevel makes, so the two describe
        the same data.

        :param levelValue: value of the table level.
        :type levelValue: float
        :return: controlling variable values per flamelet type.
        :rtype: dict[str, np.ndarray[float]]
        """
        if self.__flamelet_data_by_type is None:
            self.__flamelet_data_by_type = self.__readFlameletDataPerType()

        bandwidth = self.__levelBandwidth()
        data_per_type = {}
        for flamelet_type, cv_data in self.__flamelet_data_by_type.items():
            if bandwidth is not None:
                cv_data = cv_data[np.absolute(cv_data[:, 2] - levelValue) <= bandwidth]
            if len(cv_data) > 0:
                data_per_type[flamelet_type] = cv_data
        return data_per_type

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
        bandwidth = self.__levelBandwidth()
        if bandwidth is not None:
            level_values = flamelet_data[self._table_cv_names[2]].to_numpy()
            flamelet_data = flamelet_data[np.absolute(level_values - levelValue) <= bandwidth]

        if len(flamelet_data) < self.__min_flamelet_data_points:
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

        # Sweeping the progress variable yields both enthalpy branches at once, as the lowest and the
        # highest data point of each bin.
        along_pv = self.__sweepEnvelope(pv_data, h_data, self.__N_boundary_bins)
        if along_pv is None:
            return None
        (pv_at_lowest_h, lowest_h), (pv_at_highest_h, highest_h) = along_pv

        pv_lower, h_lower = self.__envelopeBranch(pv_at_lowest_h, lowest_h, pv_data, h_data, \
                                                  encloses_from_below=True)
        pv_upper, h_upper = self.__envelopeBranch(pv_at_highest_h, highest_h, pv_data, h_data, \
                                                  encloses_from_below=False)

        envelope = {"lower": np.column_stack((pv_lower, h_lower)),
                    "upper": np.column_stack((pv_upper, h_upper)),
                    "equilibrium": None}

        # The equilibrium branch is swept along the enthalpy, so its stations are enthalpies and the
        # extreme it bounds is a progress variable. Only the highest of each bin is of interest.
        along_h = self.__sweepEnvelope(h_data, pv_data, self.__N_transverse_bins)
        if along_h is not None:
            _, (h_at_highest_pv, highest_pv) = along_h
            h_equilibrium, pv_equilibrium = self.__envelopeBranch(h_at_highest_pv, highest_pv, h_data, pv_data, \
                                                                  encloses_from_below=False)
            envelope["equilibrium"] = np.column_stack((pv_equilibrium, h_equilibrium))

        self.__envelope_cache[levelValue] = envelope
        return envelope

    def __envelopeBranch(self, stations:np.ndarray[float], extremes:np.ndarray[float], \
                         data_along:np.ndarray[float], data_across:np.ndarray[float], encloses_from_below:bool):
        """Turn the extremes collected by a sweep into the stations of one envelope branch.

        A branch is first reduced to the stations needed to describe it, and then given back whatever
        stations that reduction, or the sampling itself, left data outside of. Enclosure is therefore
        settled against the flamelet data rather than against a bin width, so the result does not shift
        with the number of bins.

        :param stations: coordinate of the extreme data point collected in each bin.
        :type stations: np.ndarray[float]
        :param extremes: transverse coordinate of that same data point.
        :type extremes: np.ndarray[float]
        :param data_along: sweep coordinate of the flamelet data.
        :type data_along: np.ndarray[float]
        :param data_across: transverse coordinate of the flamelet data.
        :type data_across: np.ndarray[float]
        :param encloses_from_below: whether the branch bounds the data from below.
        :type encloses_from_below: bool
        :return: the stations of the branch and the branch value at each of them.
        :rtype: tuple[np.ndarray[float]]
        """
        # Simplifying drops stations, which can leave data outside; enclosing inserts stations, which
        # leaves the branch more detailed than it needs to be. Alternating settles on a branch that is
        # both, and ending on the enclosing pass is what makes the result binding.
        for _ in range(self.__branch_refinement_rounds):
            retained = self.__simplifyEnvelope(stations, extremes, encloses_from_below)
            stations, extremes = self.__encloseData(stations[retained], extremes[retained], \
                                                    data_along, data_across, encloses_from_below)
        return stations, extremes

    def __encloseData(self, stations:np.ndarray[float], values:np.ndarray[float], \
                      data_along:np.ndarray[float], data_across:np.ndarray[float], encloses_from_below:bool):
        """Extend a branch onto every data point it leaves outside, so that it encloses them all.

        A branch samples the boundary at its stations and runs straight between them, so data lying in
        between can fall on the wrong side of a segment. Rather than displacing the branch outward by a
        fixed margin, which would tie it to the bin width, the points that fall outside are made stations
        of the branch themselves. The branch then passes through them exactly and encloses the rest.

        :param stations: stations of the branch.
        :type stations: np.ndarray[float]
        :param values: branch value at each station.
        :type values: np.ndarray[float]
        :param data_along: sweep coordinate of the flamelet data.
        :type data_along: np.ndarray[float]
        :param data_across: transverse coordinate of the flamelet data.
        :type data_across: np.ndarray[float]
        :param encloses_from_below: whether the branch bounds the data from below.
        :type encloses_from_below: bool
        :return: the stations of the enclosing branch and its value at each of them.
        :rtype: tuple[np.ndarray[float]]
        """
        # Data beyond the ends of the branch is not enclosed by it and is left to the branches that
        # close the perimiter there.
        spanned = (data_along >= stations[0]) & (data_along <= stations[-1])
        data_along, data_across = data_along[spanned], data_across[spanned]

        # Inserting a station changes the two segments meeting at it, which can leave further points
        # outside; a handful of passes settles it, and each pass strictly encloses more than the last.
        for _ in range(self.__max_enclosure_passes):
            branch_at_data = np.interp(data_along, stations, values)
            outside = (data_across < branch_at_data) if encloses_from_below else (data_across > branch_at_data)
            if not np.any(outside):
                break

            stations = np.append(stations, data_along[outside])
            values = np.append(values, data_across[outside])

            # Ordered along the sweep, and within a shared station by value, so that the outermost of
            # any points sitting at the same station comes last. Points sharing a station would
            # otherwise make the branch ambiguous.
            order = np.lexsort((-values if encloses_from_below else values, stations))
            stations, values = stations[order], values[order]

            keep = np.append(np.diff(stations) > 0, True)
            stations, values = stations[keep], values[keep]

        return stations, values

    def __sweepEnvelope(self, coordinate_along:np.ndarray[float], coordinate_across:np.ndarray[float], N_bins:int):
        """Bin the data along one coordinate and collect the extreme data point of the other in each bin.

        A station is placed on the extreme data point itself, at both of its coordinates. Placing it at
        the centre of its bin instead displaces the envelope outward by however far the extreme sits from
        that centre, which is a bias of up to half a bin width that only disappears as the bin count
        grows. Following the data points removes it, and makes the result far less sensitive to the bin
        count: the bins then only decide how finely the boundary is sampled, not where it lies.

        The lower and upper extreme of a bin are two different data points at two different positions
        along the sweep, so each branch carries its own stations.

        :param coordinate_along: coordinate to bin the data along.
        :type coordinate_along: np.ndarray[float]
        :param coordinate_across: coordinate to collect the extremes of.
        :type coordinate_across: np.ndarray[float]
        :param N_bins: number of bins to sweep over.
        :type N_bins: int
        :return: the lower and upper branch, each as its own stations and values, or None if too few bins
            are populated to describe a branch.
        :rtype: tuple[tuple[np.ndarray[float]]]
        """
        station_min, station_max = np.min(coordinate_along), np.max(coordinate_along)
        if station_max <= station_min:
            return None

        bin_edges = np.linspace(station_min, station_max, N_bins + 1)
        bin_index = np.clip(np.digitize(coordinate_along, bin_edges) - 1, 0, N_bins - 1)

        extreme_lower, extreme_upper = [], []
        for iBin in range(N_bins):
            data_in_bin = np.flatnonzero(bin_index == iBin)
            # Bins collecting only the tail of a single flamelet are not representative of the boundary.
            if len(data_in_bin) >= self.__min_points_per_bin:
                extreme_lower.append(data_in_bin[np.argmin(coordinate_across[data_in_bin])])
                extreme_upper.append(data_in_bin[np.argmax(coordinate_across[data_in_bin])])

        if len(extreme_lower) < self.__min_envelope_stations:
            return None

        # Bins are disjoint and swept in order, so the stations of a branch are already ascending.
        return self.__branchOf(coordinate_along, coordinate_across, extreme_lower, encloses_from_below=True), \
               self.__branchOf(coordinate_along, coordinate_across, extreme_upper, encloses_from_below=False)

    def __branchOf(self, coordinate_along:np.ndarray[float], coordinate_across:np.ndarray[float], \
                   stations:list[int], encloses_from_below:bool):
        """Assemble one branch from the data points a sweep collected, spanning the full extent of the data.

        A bin holding too little data contributes no station, and the outermost bins are the ones that
        run out of data first. The branch would then stop short of the data it is meant to bound, so the
        outermost data points along the sweep are made stations of it regardless.

        :param coordinate_along: sweep coordinate of the flamelet data.
        :type coordinate_along: np.ndarray[float]
        :param coordinate_across: transverse coordinate of the flamelet data.
        :type coordinate_across: np.ndarray[float]
        :param stations: indices of the data points collected by the sweep, ordered along it.
        :type stations: list[int]
        :param encloses_from_below: whether the branch bounds the data from below.
        :type encloses_from_below: bool
        :return: the stations of the branch and the branch value at each of them.
        :rtype: tuple[np.ndarray[float]]
        """
        outermost = self.__outermostData(coordinate_along, coordinate_across, encloses_from_below)
        stations = [outermost[0]] + list(stations) + [outermost[1]]

        stations = np.array(stations)
        along, across = coordinate_along[stations], coordinate_across[stations]

        # An outermost point that a bin already contributed is now present twice.
        keep = np.append(np.diff(along) > 0, True)
        return along[keep], across[keep]

    def __outermostData(self, coordinate_along:np.ndarray[float], coordinate_across:np.ndarray[float], \
                        encloses_from_below:bool):
        """Locate the data points at either end of a sweep, taking the outermost where several share an end.

        :return: index of the data point at the start and at the end of the sweep.
        :rtype: tuple[int]
        """
        take_outermost = np.argmin if encloses_from_below else np.argmax

        at_start = np.flatnonzero(coordinate_along == np.min(coordinate_along))
        at_end = np.flatnonzero(coordinate_along == np.max(coordinate_along))
        return at_start[take_outermost(coordinate_across[at_start])], \
               at_end[take_outermost(coordinate_across[at_end])]

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
        tolerance = self.__simplification_tolerance * self._base_cell_size
        N_stations = len(pv_envelope)
        if tolerance <= 0 or N_stations < self.__min_envelope_stations:
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

    def _createBoundaryPolylineForTableLevel(self, levelValue:float):
        # A plane running from the reactant to the product state is bounded by states that are
        # evaluated rather than sampled, so its perimiter is built from those directly. Any other
        # plane is bounded by the flamelet data alone.
        if self.__spansReactantsToProducts():
            return self.__parametricBoundary(levelValue)

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
            corner_lower = polylineIntersection(equilibrium_branch, lower_branch)
            corner_upper = polylineIntersection(equilibrium_branch, upper_branch)
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
        merge_tolerance = self.__merge_tolerance * self._base_cell_size

        retained_points = [polyline[0]]
        for point in polyline[1:]:
            if np.linalg.norm(point - retained_points[-1]) > merge_tolerance:
                retained_points.append(point)

        # The closing segment is added by the mesher, so the perimiter should not end where it started.
        if len(retained_points) > self.__min_perimiter_points and \
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

        # The refinement field is sampled on the whole lattice rather than on the part of it that
        # ends up inside the table level: the mesher asks for the cell size around the perimiter too.
        self.__reference_lattice = CV_grid_scaled

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

        # The branches carry their own stations, so the plane they bound together is where they overlap.
        # Beyond that only one of them is defined and np.interp would extend it as a flat line.
        pv_start = max(lower_branch[0, 0], upper_branch[0, 0])
        pv_end = min(lower_branch[-1, 0], upper_branch[-1, 0])
        within_progress_variable_range = (pv_nodes >= pv_start) & (pv_nodes <= pv_end)
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

    def _passRefinementOptions(self, mesher:Mesh2DPlane):
        super()._passRefinementOptions(mesher)

        # Sampling the refinement on the reference lattice, once, replaces the per point evaluations
        # the mesher would otherwise make while it decides where to place nodes. Those run into the
        # hundreds of thousands for a single table level.
        if self.__reference_lattice is not None:
            mesher.setRefinementField(self.__reference_lattice, \
                                      self._refinementFactorField(self.__reference_lattice))
        return

    def _refinementFactorField(self, plane_nodes:np.ndarray[float]):
        refinement_factor = super()._refinementFactorField(plane_nodes)

        if self.__refineGradients:
            refinement_factor = np.minimum(refinement_factor, self.__gradientRefinementField(plane_nodes))

        if self.__refine_equilibrium:
            refinement_factor = np.minimum(refinement_factor, self.__equilibriumRefinementField(plane_nodes))

        return refinement_factor

    def __gradientRefinementField(self, plane_nodes:np.ndarray[float]):
        """Refinement towards the steepest gradients, evaluated for every node at once.

        :param plane_nodes: table plane coordinates to evaluate, one node per row.
        :type plane_nodes: np.ndarray[float]
        :return: refinement coefficient of each node.
        :rtype: np.ndarray[float]
        """
        jacobian = self._fluid_data_interpolator.Jacobian(plane_nodes[:, :self._nDim_table], \
                                                          varnames=self._refine_for_gradients_of)
        jacobian = np.reshape(jacobian, (self._nDim_table, len(plane_nodes), len(self._refine_for_gradients_of)))

        gradient_magnitude = np.linalg.norm(jacobian, axis=0) / self._gradient_norm_factors
        coefficients = np.asarray(self._gradient_refinement_factors)

        # Scale from no refinement where the quantity is flat, down to the coefficient where its
        # gradient is the steepest found on the level; the strongest demand over the quantities wins.
        per_quantity = np.maximum(coefficients, 1.0 - (1.0 - coefficients)*gradient_magnitude)
        return np.minimum(np.min(per_quantity, axis=1), 1.0)

    def __equilibriumRefinementField(self, plane_nodes:np.ndarray[float]):
        """Refinement in proximity of the reactants and products, evaluated for every node at once.

        :param plane_nodes: table plane coordinates to evaluate, one node per row.
        :type plane_nodes: np.ndarray[float]
        :return: refinement coefficient of each node.
        :rtype: np.ndarray[float]
        """
        refinement_factor = np.ones(len(plane_nodes))
        if self._is2D():
            near_reactants_or_products = (plane_nodes[:, 0] <= self.__margin_equilibrium) \
                                       | (plane_nodes[:, 0] >= (1 - self.__margin_equilibrium))
            refinement_factor[near_reactants_or_products] = self.__equilibrium_refinement_factor
        return refinement_factor

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