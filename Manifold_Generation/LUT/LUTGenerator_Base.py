###############################################################################################
#       #      _____ __  _____      ____        __        __  ____                   #        #
#       #     / ___// / / /__ \    / __ \____ _/ /_____ _/  |/  (_)___  ___  _____   #        #
#       #     \__ \/ / / /__/ /   / / / / __ `/ __/ __ `/ /|_/ / / __ \/ _ \/ ___/   #        #
#       #    ___/ / /_/ // __/   / /_/ / /_/ / /_/ /_/ / /  / / / / / /  __/ /       #        #
#       #   /____/\____//____/  /_____/\__,_/\__/\__,_/_/  /_/_/_/ /_/\___/_/        #        #
#       #                                                                            #        #
###############################################################################################

############################# FILE NAME: LUTGenerator_Base.py #################################
#=============================================================================================#
# author: Evert Bunschoten                                                                    |
#    :PhD Candidate ,                                                                         |
#    :Flight Power and Propulsion                                                             |
#    :TU Delft,                                                                               |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#  Base class for tabulated methods in SU2 DataMiner                                          |
#                                                                                             |
# Version: 3.2.0                                                                              |
#                                                                                             |
#=============================================================================================#
import numpy as np
import pandas as pd
import meshio
import pandas
from os import sep
from copy import copy
from sklearn.preprocessing import MinMaxScaler
from multiprocessing import Pool, get_context


from scipy.interpolate import RBFInterpolator

from Common.Interpolators import fluidDataInterpolator
from Common.DataDrivenConfig import Config
from Manifold_Generation.LUT.MeshTools import Mesh2DPlane


def createWorkerPool(n_workers:int):
    """Create a worker pool with fork semantics wherever the platform allows it.

    Linux defaults to forking, but macOS and Windows default to the "spawn" start method, under
    which every worker re-imports the ``__main__`` module. A driver script that builds and runs
    the table generator at module level (the usual style) therefore re-executes itself in each
    worker, which then fails while bootstrapping; the pool replaces the dead worker and the cycle
    repeats without meshing a single level. Avoiding that under spawn would require every driver
    script to be wrapped in an ``if __name__ == "__main__":`` guard. Falls back to the platform
    default where fork is unavailable.

    :param n_workers: number of worker processes.
    :type n_workers: int
    :return: process pool.
    """
    try:
        return get_context("fork").Pool(n_workers)
    except (ValueError, OSError):
        return Pool(n_workers)


class SU2TableGenerator_Base:
    
    _Config:Config = None
    _nDim_table:int=2
    _base_cell_size:float = 2e-2
    _target_number_of_nodes:int = None

    _table_vars:list[str] = []
    _table_nodes:list[np.ndarray[float]] = []
    _data_in_table:list[np.ndarray[float]] = []
    _table_connectivity:list[np.ndarray[int]] = []
    _table_hullnodes:list[np.ndarray[int]]  = []

    _conditional_refinement_vars:list[str] = []     # Thermophysical quantities for which refinement is applied within specified bounds
    _conditional_refinement_indices:list[int] = []  # Column indices of refinement quantities
    _conditional_lower_bound:list[float] = []       # Lower bounds of conditional refinement quantities.
    _conditional_upper_bound:list[float] = []       # Upper bounds of conditional refinement quantities.
    _conditional_refinement_factor:list[float] = [] # Refinement factors to apply within the bounds

    _normalized_refinement_dims:list[int] = []      # Controlling variable indices for which refinement is applied within a normalized range
    _normalized_lower_bound:list[float] = []        # Lower bounds of the normalized refinement ranges.
    _normalized_upper_bound:list[float] = []        # Upper bounds of the normalized refinement ranges.
    _normalized_refinement_factor:list[float] = []  # Refinement factors to apply within the normalized ranges

    _level_cv_lower:np.ndarray[float] = None        # Per-level lower bounds of the scaled controlling variables.
    _level_cv_span:np.ndarray[float] = None         # Per-level extent of the scaled controlling variables.

    _scaler_controlling_variables:MinMaxScaler = MinMaxScaler()
    _fluid_data_interpolator:fluidDataInterpolator = None

    __smoothTableData:bool = False
    __smoothingLevel:float = 0.1

    __anchorHullEdges:bool = True
    _hull_edge_anchors:list = []

    __N_nearest_neighbors:int = None
    __inverse_distance_exponent:float = None

    _N_table_levels:int = None
    _table_levels:np.ndarray[float] = []
    _tableUpperLevel:float = None
    _tableLowerLevel:float = None
    _table_level_inserts:list[float] = []

    __run_parallel:bool = False
    __N_cores:int = 1

    _state_quantities:list[str] = []

    _planarMesher:Mesh2DPlane = Mesh2DPlane()

    _verbosity:int=1

    def __init__(self, config_in:Config):
        """Initialize the table generator and retrieve settings from configuration

        :param config_in: SU2 DataMiner configuration object.
        :type config_in: Config
        """
        self._Config = copy(config_in)
        self._nDim_table = len(self._Config.GetControllingVariables())
        return
    
    def generateTable(self):
        """Initiate the table generation process.
        """
        self._processTableLevels()

        self._defineFluidDataInterpolator()

        self._generateTableLevelData()

        self._applyHullEdgeAnchors()

        if self.__smoothTableData:
            self._smoothTableLevelData()

        self._postProcessTableData()
        return

    def _varsExcludedFromSmoothing(self):
        """Table variables that must not be touched by the smoothing pass.

        The controlling variables define the node coordinates and are never smoothed.
        Subclasses extend this for quantities whose shape smoothing would destroy.

        :return: variable names to leave untouched.
        :rtype: list[str]
        """
        return list(self._Config.GetControllingVariables())

    def _unsmoothedNodes(self, level_index:int):
        """Table nodes on this level that the smoothing pass must reproduce exactly.

        Smoothing denoises nodes that have neighbours on all sides; where every neighbour lies
        on one side its stencil is one-sided and it biases instead. Subclasses return the nodes
        of such regions. None by default.

        :param level_index: table level index.
        :type level_index: int
        :return: node indices.
        :rtype: np.ndarray[int]
        """
        return np.zeros(0, dtype=int)

    def _hullEdgeAnchors(self, level_index:int):
        """Table nodes on this level whose state is known exactly, and the values they carry.

        Subclasses return a mapping from variable name to a pair of (node indices, values), so
        that different edges of the same level can anchor different quantities. Returning None
        anchors nothing, which is the default.

        :param level_index: table level index.
        :type level_index: int
        :return: {variable name: (node indices, values)}, or None.
        :rtype: dict[str, tuple[np.ndarray[int], np.ndarray[float]]] | None
        """
        return None

    def _applyHullEdgeAnchors(self):
        """Restate the exactly known values on the perimeter of every table level.

        Both stages that produce level data have one-sided support on the perimeter: the
        inverse-distance interpolation draws on the nearest point-cloud samples and the
        smoothing pass on the nearest table nodes, and on an edge every one of those neighbours
        lies inwards. A quantity that vanishes on the edge but is large just inside it is
        therefore pulled up towards the interior. The reaction source is the worst case, being
        both exactly zero on the chemically frozen edges and steep immediately inside them:
        measured on the H2 manifold the burnt edge ended up reacting at roughly 400 times the
        rate of the point cloud the table was built from, and a CFD solver reads precisely those
        nodes whenever a state is clamped back onto the manifold.

        Anchoring runs before smoothing so that the smoother sees the corrected values, and the
        anchored nodes are then fitted without regularisation so it reproduces rather than
        re-smears them.
        """
        self._hull_edge_anchors = [None] * self._N_table_levels
        if not self.__anchorHullEdges:
            return

        n_anchored = 0
        for iLevel in range(self._N_table_levels):
            anchors = self._hullEdgeAnchors(iLevel)
            if not anchors:
                continue
            anchors = {var: (np.asarray(idx, dtype=int), np.asarray(val, dtype=float))
                       for var, (idx, val) in anchors.items()
                       if var in self._data_in_table[iLevel] and len(idx) > 0}
            if not anchors:
                continue
            self._hull_edge_anchors[iLevel] = anchors
            self.__writeAnchors(iLevel, anchors)
            n_anchored += len(self.__anchoredNodeIndices(anchors))

        if n_anchored:
            self.__printmsg("Anchored %i hull edge nodes over %i table levels" % (n_anchored, self._N_table_levels))
        return

    @staticmethod
    def __anchoredNodeIndices(anchors:dict):
        """Union of the node indices touched by an anchor set."""
        return np.unique(np.concatenate([idx for idx, _ in anchors.values()]))

    def __writeAnchors(self, level_index:int, anchors:dict):
        """Write anchored values into a level's data frame by position."""
        level_data = self._data_in_table[level_index]
        for var, (idx, val) in anchors.items():
            level_data.iloc[idx, level_data.columns.get_loc(var)] = val
        return

    def _postProcessTableData(self):
        """Hook for corrections that must be applied to the finished table data.

        Called after the level data have been interpolated and smoothed, so anything applied
        here sees the values that are actually written to the table. Does nothing by default.
        """
        return
    
    def _processTableLevels(self):
        """Prepare table level values depending on the number of dimensions.

        :raises Exception: if no information is specified regarding the table level interval and number of table levels.
        """
        if self._is3D():
            if not any(self._table_levels):
                if self._tableUpperLevel and self._tableLowerLevel and self._N_table_levels:
                    self._table_levels = np.linspace(self._tableLowerLevel, self._tableUpperLevel, self._N_table_levels)
                else:
                    raise Exception("No table level information provided, aborting")

            if any(self._table_level_inserts):
                self._table_levels = np.append(self._table_levels, np.array(self._table_level_inserts))

            unique_table_levels = np.unique(self._table_levels)
            sorted_table_levels = np.sort(unique_table_levels)
            self._table_levels = sorted_table_levels
            self._N_table_levels = len(self._table_levels)
            self._tableLowerLevel = min(self._table_levels)
            self._tableUpperLevel = max(self._table_levels)
        else:
            self._table_levels = [0]
            self._N_table_levels = 1
        return
    
    
    def _defineFluidDataInterpolator(self):
        """Prepare the interpolation function for evaluating thermochemical states on the table nodes and scale the controlling variables in the state space.
        """
        stateDataFrame = self._getFluidDataForInterpolator()
        cv_data = np.column_stack(tuple(stateDataFrame[cv] for cv in self._Config.GetControllingVariables()))
        cv_data_scaled = self._scaler_controlling_variables.fit_transform(cv_data)
        self._fluid_data_interpolator = fluidDataInterpolator(cv_data_scaled, stateDataFrame, self.__N_nearest_neighbors, self.__inverse_distance_exponent)
        return
    
    def _getFluidDataForInterpolator(self):
        return
    
    def _generateTableLevelData(self):
        """Generate the connectivity for each table level and interpolate the thermochemical state data onto the table nodes.
        """
        self._table_nodes = [None]*self._N_table_levels
        self._table_connectivity = [None]*self._N_table_levels
        self._table_hullnodes = [None]*self._N_table_levels
        self._data_in_table = [None]*self._N_table_levels
        if self.__run_parallel:
            pool = createWorkerPool(self.__N_cores)
            results = pool.map(self.meshTableLevel, [i for i in range(self._N_table_levels)])
            pool.close()

            for iLevel in range(self._N_table_levels):
                self._table_nodes[iLevel] = results[iLevel][0]
                self._table_connectivity[iLevel] = results[iLevel][1]
                self._table_hullnodes[iLevel] = results[iLevel][2]
                self._data_in_table[iLevel] = results[iLevel][3]
        else:
            for iLevel in range(len(self._table_levels)):
                tableLevel = self.meshTableLevel(iLevel)
                
                self._table_nodes[iLevel] = tableLevel[0]
                self._table_connectivity[iLevel] = tableLevel[1]
                self._table_hullnodes[iLevel] = tableLevel[2]
                self._data_in_table[iLevel] = tableLevel[3]
        return
    
    def _smoothTableLevelData(self):
        self.__printmsg("Smoothing table data...")
        for iLevel in range(self._N_table_levels):
            cv_level_scaled = self._scaler_controlling_variables.transform(self._table_nodes[iLevel])
            
            data_on_table_level = self._data_in_table[iLevel].values

            # Anchored nodes are fitted exactly (zero regularisation) so that the smoother
            # reproduces the values written by _applyHullEdgeAnchors instead of pulling them
            # back towards the interior. Regularisation is per node rather than per variable,
            # so an anchored node keeps its interpolated values for every other quantity too --
            # which is what is wanted on a perimeter, where the smoother's stencil is one-sided
            # and it biases rather than denoises.
            smoothing = np.full(len(cv_level_scaled), self.__smoothingLevel, dtype=float)
            anchors = self._hull_edge_anchors[iLevel] if self._hull_edge_anchors else None
            if anchors:
                smoothing[self.__anchoredNodeIndices(anchors)] = 0.0
            smoothing[self._unsmoothedNodes(iLevel)] = 0.0

            smoothener = RBFInterpolator(cv_level_scaled[:, :2], data_on_table_level, kernel="linear",smoothing=smoothing,neighbors=100)
            smoothened_table_data = smoothener(cv_level_scaled[:, :2])
            excluded = self._varsExcludedFromSmoothing()
            for ivar, var in enumerate(list(self._data_in_table[iLevel].keys())):
                if var not in excluded:
                    self._data_in_table[iLevel][var] = smoothened_table_data[:, ivar]

            # Restated afterwards as well: the zero-regularisation fit already reproduces them,
            # but this keeps the guarantee independent of how well conditioned that fit is.
            if anchors:
                self.__writeAnchors(iLevel, anchors)
        self.__printmsg("Done")
        return
    
    def setMaximumCellSize(self, cell_size_coarse:float=1e-2):
        """Specify the coarse level cell size of the table

        :param cell_size_coarse: coarse cell size, defaults to 1e-2
        :type cell_size_coarse: float, optional
        :raises Exception: if specified cell size is negative or zero
        """
        if cell_size_coarse <= 0:
            raise Exception("Maximum cell size should be strictly positive.")
        
        self._base_cell_size = cell_size_coarse
        return
    
    def setTargetNodeCount(self, target_node_count:int=3000):
        """Specify a target number of nodes for each table level. The table generator aims to approximate the target number of nodes within 1%.

        :param target_node_count: number of nodes on each table level, defaults to 3000
        :type target_node_count: int, optional
        :raises Exception: if the specified value is not strictly positive.
        """
        if target_node_count <= 0:
            raise Exception("Target number of nodes in the table should be strictly positive.")
        self._target_number_of_nodes = target_node_count
        return
    
    def applyRefinementWithin(self, varname:str, lowerbound:float=-np.inf, upperbound:float=np.inf, coef:float=0.5):
        """Specify conditional refinement based on interpolated thermochemical state data. Cell sizes are reduced by a factor "coef" where the table data lies between the specified lower bound and upper bound (inclusive)

        :param varname: thermochemical state variable name.
        :type varname: str
        :param lowerbound: lower bound above which refinement is applied, defaults to -np.inf
        :type lowerbound: float, optional
        :param upperbound: upper bound below which refinement is applied, defaults to np.inf
        :type upperbound: float, optional
        :param coef: mesh refinement factor, defaults to 0.5
        :type coef: float, optional
        :raises Exception: if variable is not in the set of flamelet thermochemical state variables.
        :raises Exception: if lower bound value exceeds upper bound value.
        :raises Exception: if the coefficient value is negative.
        """
        
        if varname not in self._state_quantities:
            raise Exception("%s is not in the list of available thermophysical state variables" % varname)
        if lowerbound > upperbound:
            raise Exception("Upper bound value should exceed lower bound value")
        if coef <= 0:
            raise Exception("Refinement coeffcient should be positive")
        self._conditional_refinement_vars.append(varname)
        self._conditional_refinement_indices.append(self._state_quantities.index(varname))
        self._conditional_lower_bound.append(lowerbound)
        self._conditional_upper_bound.append(upperbound)
        self._conditional_refinement_factor.append(coef)
        return

    def applyRefinementWithinNormalized(self, varname:str, lowerbound:float=0.0, upperbound:float=1.0, coef:float=0.5):
        """Specify refinement within a range of a controlling variable, normalized per table level.

        Unlike :meth:`applyRefinementWithin`, which tests an absolute value, the bounds here are
        fractions of the extent that the controlling variable actually spans on the table level
        being meshed. Zero is the lower edge of that level and one is the upper edge.

        This is what makes edge-following refinement expressible on a 3D table. The controlling
        variables are scaled once, globally, over the whole point cloud, so a fixed coordinate
        does not correspond to the same feature on every level: the burnt progress variable is the
        equilibrium value at that level's mixture fraction, and on a lean or rich level it sits
        far below the global maximum that scaled coordinate 1 refers to. Testing the normalized
        coordinate instead follows the equilibrium edge across the levels. Refining the burnt 20%
        of every level is then::

            generator.applyRefinementWithinNormalized("ProgressVariable", 0.8, 1.0, 0.5)

        :param varname: controlling variable name.
        :type varname: str
        :param lowerbound: lower bound of the normalized range, defaults to 0.0
        :type lowerbound: float, optional
        :param upperbound: upper bound of the normalized range, defaults to 1.0
        :type upperbound: float, optional
        :param coef: mesh refinement factor, defaults to 0.5
        :type coef: float, optional
        :raises Exception: if the variable is not a controlling variable of the table.
        :raises Exception: if the bounds lie outside [0, 1] or the lower bound exceeds the upper bound.
        :raises Exception: if the coefficient value is negative.
        """
        controlling_variables = list(self._Config.GetControllingVariables())
        if varname not in controlling_variables:
            raise Exception("%s is not a controlling variable of this table, options are %s" % (varname, ", ".join(controlling_variables)))
        if lowerbound < 0.0 or upperbound > 1.0:
            raise Exception("Normalized refinement bounds should lie between zero and one")
        if lowerbound > upperbound:
            raise Exception("Upper bound value should exceed lower bound value")
        if coef <= 0:
            raise Exception("Refinement coeffcient should be positive")

        iDim = controlling_variables.index(varname)
        if iDim >= self._nDim_table:
            raise Exception("%s is not a controlling variable of this table" % varname)
        if self._is3D() and iDim == 2:
            raise Exception("Normalized refinement is meaningless along the table level direction %s, since a level spans a single value of it" % varname)

        self._normalized_refinement_dims.append(iDim)
        self._normalized_lower_bound.append(lowerbound)
        self._normalized_upper_bound.append(upperbound)
        self._normalized_refinement_factor.append(coef)
        return

    def _setLevelBounds(self, pointCloud_scaled:np.ndarray[float]):
        """Record the extent that each controlling variable spans on the level being meshed.

        The point cloud bounds the meshed domain -- the perimeter is a concave hull of it -- so
        its bounding box is the reference against which normalized refinement positions are
        measured.

        :param pointCloud_scaled: scaled point cloud coordinates of the current table level.
        :type pointCloud_scaled: np.ndarray[float]
        """
        finite_pts = pointCloud_scaled[np.isfinite(pointCloud_scaled).all(axis=1)]
        if len(finite_pts) == 0:
            self._level_cv_lower = None
            self._level_cv_span = None
            return
        self._level_cv_lower = finite_pts.min(axis=0)
        self._level_cv_span = finite_pts.max(axis=0) - self._level_cv_lower
        return

    def _normalizedLevelCoord(self, iDim:int, coords:np.ndarray[float]):
        """Position along a controlling variable as a fraction of the current level's extent.

        :param iDim: controlling variable index.
        :type iDim: int
        :param coords: scaled controlling variable values (x, y, z).
        :type coords: np.ndarray[float]
        :return: normalized position, or None if the level extent is unknown or degenerate.
        :rtype: float | None
        """
        if self._level_cv_span is None or iDim >= len(self._level_cv_span):
            return None
        if self._level_cv_span[iDim] <= 0.0:
            return None
        return (coords[iDim] - self._level_cv_lower[iDim]) / self._level_cv_span[iDim]

    def _applyNormalizedRefinement(self, x:float, y:float, z:float):
        """Evaluate the refinement coefficient from the normalized controlling variable ranges.

        :return: refinement coefficient
        :rtype: float
        """
        refinement_factor = 1.0
        if len(self._normalized_refinement_dims) == 0:
            return refinement_factor
        coords = np.array([x, y, z])
        for iDim, lower, upper, c in zip(self._normalized_refinement_dims, self._normalized_lower_bound, \
                                         self._normalized_upper_bound, self._normalized_refinement_factor):
            normalized_coord = self._normalizedLevelCoord(iDim, coords)
            if normalized_coord is None:
                continue
            if (normalized_coord >= lower) and (normalized_coord <= upper):
                refinement_factor = min(refinement_factor, c)
        return refinement_factor

    def meshTableLevel(self, level_index:int):
        """Discretize and interpolate the fluid data of a single table level.

        :param level_index: table level index
        :type level_index: int
        :return: list with table nodes, connectivity, perimiter indices, and interpolated fluid data.
        :rtype: list[np.ndarray[float]]
        """

        pointCloud = self._createPointCloudForTableLevel(self._table_levels[level_index])
        self._setLevelBounds(pointCloud)

        mesher = self._initiateMesher()
        if self.__run_parallel:
            mesher.setVerbosity(0)
        else:
            if self._verbosity > 1:
                mesher.setVerbosity(1)
            if self._verbosity > 2:
                mesher.setVerbosity(2)

        mesher.setInitialPointCloud(pointCloud)

        self._passRefinementOptions(mesher)

        mesher.generateMesh()

        cvTable, triangles, hullIDs, tableDataFrame = self.__extractTableLevelData(mesher)
        
        if self._is3D():
            self.__printmsg("Finished meshing table level %i at %s=%.2e with %i nodes" % (level_index, \
                                                                                    self._Config.GetControllingVariables()[2], \
                                                                                    self._table_levels[level_index], \
                                                                                    len(cvTable)))

        return [cvTable, triangles,hullIDs, tableDataFrame]
    
    def _createPointCloudForTableLevel(self, levelValue:float):
        """Generate a planar point cloud used as a reference for discretizing the table level.

        :param levelValue: value of the third table dimension corresponding to the table level.
        :type levelValue: float
        """
        return
    
    def _initiateMesher(self):
        return Mesh2DPlane()
    
    def _passRefinementOptions(self, mesher:Mesh2DPlane):
        """Transfer table refinement information to 2D meshing algorithm.

        """
        mesher.setBaseCellSize(self._base_cell_size)
        if self._target_number_of_nodes:
            mesher.setTargetNodeCount(self._target_number_of_nodes)
        mesher.setRefinementFunction(self._refinelocation)
        return
    
    def __extractTableLevelData(self, mesher:Mesh2DPlane):
        """Retrieve table level nodes and connectivity from 2D meshing algorithm.

        :param mesher: 2D meshing tool
        :type mesher: Mesh2DPlane
        :return: table nodes, connectivity, perimiter indices, and interpolated fluid data.
        :rtype: tuple
        """
        meshnodes = mesher.getMeshNodes()
        triangles = mesher.getConnectivity()
        hullIDs = mesher.getHullIDs()
        cvTable_scaled = np.zeros([np.shape(meshnodes)[0], self._nDim_table])
        for iCv in range(self._nDim_table):
            cvTable_scaled[:, iCv] = meshnodes[:, iCv]

        table_state_data = self._calculateTableStateData(cvTable_scaled)

        tableDataFrame = pd.DataFrame()
        for var in self._table_vars:
            tableDataFrame[var] = table_state_data[:, self._state_quantities.index(var)]


        cvTable = self._scaler_controlling_variables.inverse_transform(cvTable_scaled)
        for iCv, cv in enumerate(self._Config.GetControllingVariables()):
            tableDataFrame[cv] = cvTable[:, iCv]
        
        return cvTable, triangles, hullIDs, tableDataFrame
    
    def _calculateTableStateData(self, cv_table_nodes:np.ndarray[float]):
        """Interpolate fluid thermochemical data onto table nodes.

        :param cv_table_nodes: control variable values of table nodes
        :type cv_table_nodes: np.ndarray[float]
        :return: interpolated fluid data
        :rtype: pd.DataFrame
        """
        return self._fluid_data_interpolator(cv_table_nodes)
    
    def _refinelocation(self, x:float,y:float,z:float):
        """Evaluate refinement coefficient based on interpolated fluid data.

        :param x: first controlling variable.
        :type x: float
        :param y: second controlling variable.
        :type y: float
        :param z: third controlling variable.
        :type z: float
        :return: refinement coefficient
        :rtype: float
        """
        refinement_factor = self._applyNormalizedRefinement(x, y, z)
        if len(self._conditional_refinement_indices) > 0:
            cv_input = np.array([x,y,z])
            fluid_state_data = self._fluid_data_interpolator(cv_input[:self._nDim_table])
            valid_pt = len(fluid_state_data) > 0
            if valid_pt:
                for ivar, lower, upper, c in zip(self._conditional_refinement_indices, self._conditional_lower_bound, self._conditional_upper_bound, self._conditional_refinement_factor):
                    test_val = fluid_state_data[ivar]
                    if (test_val >= lower) and (test_val <= upper):
                        refinement_factor = min(refinement_factor, c)
        return refinement_factor
    
    def setTableVars(self, table_vars_in:list[str]):
        """Specify the thermocehmical state variables in the table file.

        :param table_vars_in: names of variables to include in table.
        :type table_vars_in: list[str]
        :raises Exception: if unsupported variables are included.
        """
        controlling_variables = self._Config.GetControllingVariables()
        self._table_vars = []
        for cv in controlling_variables:
            if cv not in table_vars_in:
                print("Controlling variable %s should be in table, adding" % cv)
                self._table_vars.append(cv)
        
        for user_var in table_vars_in:
            if self._checkIfVariableIsValid(user_var):
                self._table_vars.append(user_var)
            else:
                raise Exception("%s is not supported for table output variables")
        return
    
    def _checkIfVariableIsValid(self, var_to_check:str):
        return True
    
    
    def writeParaviewTable(self, filepath_out:str=None):
        """Write the table level information to vtk files.

        :param filepath_out: file name header, defaults to None
        :type filepath_out: str, optional
        """
        self.__printmsg("Writing vtk table to %s..." % filepath_out)
        for iLevel in range(self._N_table_levels):
            if self._N_table_levels > 1:
                table_level_filename = "%s_%i" % (filepath_out, iLevel)
            else:
                table_level_filename = filepath_out
            
            cv_table_nodes = self._table_nodes[iLevel]
            cv_nodes_scaled = self._scaler_controlling_variables.transform(cv_table_nodes)
            
            self.__writeParaViewVTK(cv_nodes_scaled, self._data_in_table[iLevel], self._table_connectivity[iLevel], table_level_filename)
        self.__printmsg("Done")
        return
    
    def __writeParaViewVTK(self, cv_coordinates:np.ndarray[float], table_data:pandas.DataFrame, connectivity:np.ndarray[int], file_name_out:str):
        
        pts = np.asarray(cv_coordinates, dtype=float)
        if pts.shape[1] == 2:
            pts = np.column_stack((pts, np.zeros(len(pts))))
        conn = np.asarray(connectivity, dtype=np.int64)
        if conn.min() == 1:
            conn = conn - 1
        point_data = {}
        for var in self._table_vars:
            point_data[var] = np.asarray(table_data[var])

        mesh = meshio.Mesh(
            points=pts,
            cells=[("triangle", conn)],
            point_data=point_data
        )
        mesh.write("%s.vtk" % file_name_out)

        return
    
    def _writeAdditionalInfoToTable(self, fid):
        return
    
    def writeSU2Table(self, filepath_out:str=None):
        """Write table information to SU2 drg file.

        :param filepath_out: file name header, defaults to None
        :type filepath_out: str, optional
        """
        if filepath_out:
            file_out = filepath_out + ".drg"
        else:
            file_out = sep.join((self._Config.GetOutputDir(), "LUT_%s.drg" % self._Config.GetConfigName()))

        self.__printmsg("Writing LUT file with name %s" % file_out)

        fid = open(file_out, "w+")
        fid.write("Dragon library\n\n")
        fid.write("<Header>\n\n")
        if self._is2D():
            fid.write("[Version]\n1.0.1\n\n")
        else:
            fid.write("[Version]\n1.1.0\n\n")
        
        self._writeAdditionalInfoToTable(fid)

        if self._is3D():
            fid.write("[Number of table levels]\n%i\n\n" % self._N_table_levels)
            fid.write("[Table levels]\n")
            for z in self._table_levels:
                fid.write("%+.16e\n" % z)
            fid.write("\n")

        fid.write("[Number of points]\n")
        for d in self._table_nodes:
            Np = np.shape(d)[0]
            fid.write("%i\n" % Np)
        fid.write("\n")

        fid.write("[Number of triangles]\n")
        for trias in self._table_connectivity:
            Ntria = np.shape(trias)[0]
            fid.write("%i\n" % Ntria)
        fid.write("\n")

        fid.write("[Number of hull points]\n")
        for hull in self._table_hullnodes:
            Nhull = len(hull)
            fid.write("%i\n" % Nhull)
        fid.write("\n")

        fid.write("[Number of variables]\n%i\n\n" % len(self._table_vars))
        fid.write("[Variable names]\n")
        for iVar, Var in enumerate(self._table_vars):
            fid.write("%i:%s\n" % (iVar+1, Var))
        fid.write("\n")

        fid.write("</Header>\n\n")

        fid.write("<Data>\n")
        for iLevel in range(len(self._table_nodes)):
            if self._is3D():
                fid.write("<Level>\n")
            Np = np.shape(self._table_nodes[iLevel])[0]
            for iNode in range(Np):
                line_table_data = "\t".join(["%.14e" % self._data_in_table[iLevel][var][iNode] for var in self._table_vars])
                
                fid.write(line_table_data + "\n")
            if self._is3D():
                fid.write("</Level>\n")
        fid.write("</Data>\n\n")

        fid.write("<Connectivity>\n")
        for iLevel in range(len(self._table_connectivity)):
            if self._is3D():
                fid.write("<Level>\n")
            for iCell in range(len(self._table_connectivity[iLevel])):
                fid.write("\t".join("%i" % c for c in self._table_connectivity[iLevel][iCell, :]+1) + "\n")
            if self._is3D():
                fid.write("</Level>\n")
        fid.write("</Connectivity>\n\n")

        fid.write("<Hull>\n")
        for iLevel in range(len(self._table_hullnodes)):
            if self._is3D():
                fid.write("<Level>\n")
            for iCell in range(len(self._table_hullnodes[iLevel])):
                fid.write(("%i" % (self._table_hullnodes[iLevel][iCell]+1)) + "\n")
            if self._is3D():
                fid.write("</Level>\n")
        fid.write("</Hull>\n\n")
        self.__printmsg("Done")

        fid.close()
    
    def setTableLevels(self, level_values:np.ndarray[float]):
        """Specify the level values for 3D tables.

        :param level_values: array with table level values.
        :type level_values: np.ndarray[float]
        """
        if self.__optionAppliesFor3D():
            self._table_levels = level_values
        return
    
    def setTableLimits(self, lower_limit:float, upper_limit:float):
        """Specify the upper and lower limit of the third table dimension.

        :param lower_limit: lower level value
        :type lower_limit: float
        :param upper_limit: upper level value
        :type upper_limit: float
        :raises Exception: if lower level value exceeds upper level value.
        """
        if self.__optionAppliesFor3D():
            if lower_limit >= upper_limit:
                raise Exception("Table upper level value should exceed lower level value")
            self._tableLowerLevel = lower_limit
            self._tableUpperLevel = upper_limit
        return
    
    def setNTableLevels(self, N_levels:int=10):
        """Specify the number of table levels.

        :param N_levels: number of equidistant table levels, defaults to 10
        :type N_levels: int, optional
        :raises Exception: if the specified number of levels is negative
        """
        if self.__optionAppliesFor3D():
            if N_levels <= 0:
                raise Exception("Number of table levels should be positive")
            self._N_table_levels = N_levels
        return
    
    def insertTableLevel(self, level_value:float):
        """Insert a table level for a specific value.

        :param level_value: table level value.
        :type level_value: float
        """
        if self.__optionAppliesFor3D():
            self._table_level_inserts.append(level_value)
        return
    
    def setNProcessors(self, N_cores:int=2):
        """Specify the number of parallel workers used to generate the table levels in 3D tables.

        :param N_cores: Number of processors, defaults to 2
        :type N_cores: int, optional
        :raises Exception: If fewer than 2 processors are selected.
        """
        if self.__optionAppliesFor3D():
            if N_cores <= 1:
                raise Exception("At least two cores should be used")
            self.__run_parallel = True
            self.__N_cores = N_cores
        return
    
    def setNNearestNeighbors(self, N_input:int=6):
        """Specify the number of nearest neighbors used for fluid data interpolation.

        :param N_input: number of nearest neighbors, defaults to 6
        :type N_input: int, optional
        :raises Exception: if fewer than one neighbors are selected.
        """
        if N_input < 1:
            raise Exception("Number of nearest neighbors should be strictly positive")
        self.__N_nearest_neighbors = N_input
        return
    
    def setInverseDistanceExponent(self, p_factor:float=2):
        """Specify the inverse distance exponent used for fluid data interpolation.

        :param p_factor: inverse distance exponent value, defaults to 2
        :type p_factor: float, optional
        :raises Exception: if the provided value is negative.
        """
        if p_factor <= 0:
            raise Exception("Inverse distance parameter value should be positive")
        self.__inverse_distance_exponent = p_factor
        return
    
    def setSmoothingParameter(self, smoothening_factor:float=0):
        """Apply smoothing to table data. High value = more smoothing, low value = no smoothing

        :param smoothing_factor: smoothing parameter value, defaults to 0
        :type smoothing_factor: float, optional
        """
        self.__smoothTableData = True
        self.__smoothingLevel = smoothening_factor
        return
    
    def setHullEdgeAnchoring(self, anchor_hull_edges:bool=True):
        """Restate the exactly known values on the perimeter of each table level.

        Enabled by default. Disable it only to reproduce a table generated before anchoring
        existed; see _applyHullEdgeAnchors for what it corrects.

        :param anchor_hull_edges: whether to anchor the hull edges, defaults to True
        :type anchor_hull_edges: bool, optional
        """
        self.__anchorHullEdges = anchor_hull_edges
        return

    def setVerbosity(self, verbosity_level:int=1):
        """Specify verbosity level used to print information in the terminal.

        :param verbosity_level: Verbosity level between 0 and 4, defaults to 1
        :type verbosity_level: int, optional
        :raises Exception: if specified value lies outside 0-4
        """
        if verbosity_level < 0 or verbosity_level > 4:
            raise Exception("Verbosity level should be between 0 and 4")
        self._verbosity = verbosity_level
        return
    
    def __optionAppliesFor3D(self):
        if self._is2D():
            print("Option does not apply for two-dimensional table, ignoring")
            return False
        else:
            return True

    def _is2D(self):
        return self._nDim_table==2

    def _is3D(self):
        return self._nDim_table==3
    
    
    def __printmsg(self, msg:str):
        if self._verbosity > 0:
            print(msg)
        return
    
    
    