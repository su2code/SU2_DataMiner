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
from multiprocessing import Pool
from scipy.interpolate import RBFInterpolator
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from Common.Interpolators import fluidDataInterpolator
from Common.DataDrivenConfig import Config
from Manifold_Generation.LUT.MeshTools import Mesh2DPlane


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

    _scaler_controlling_variables:MinMaxScaler = MinMaxScaler()
    _fluid_data_interpolator:fluidDataInterpolator = None

    __smoothTableData:bool = False
    __smoothingLevel:float = 0.1

    __N_nearest_neighbors:int = None
    __inverse_distance_exponent:float = None

    _table_cv_names:list[str] = None   # Controlling variables spanning the table, ordered as [plane, plane, level].
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
        self.setTableAxes(self._Config.GetControllingVariables()[:2], \
                          self._Config.GetControllingVariables()[2] if len(self._Config.GetControllingVariables()) > 2 else None)
        return
    
    def generateTable(self):
        """Initiate the table generation process.
        """
        self._processTableLevels()

        self._defineFluidDataInterpolator()

        self._generateTableLevelData()

        if self.__smoothTableData:
            self._smoothTableLevelData()
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
        cv_data = np.column_stack(tuple(stateDataFrame[cv] for cv in self._table_cv_names))
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
            pool = Pool(self.__N_cores)
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
            smoothener = RBFInterpolator(cv_level_scaled[:, :2], data_on_table_level, kernel="linear",smoothing=self.__smoothingLevel,neighbors=100)
            smoothened_table_data = smoothener(cv_level_scaled[:, :2])
            for ivar, var in enumerate(list(self._data_in_table[iLevel].keys())):
                if var not in self._Config.GetControllingVariables():
                    self._data_in_table[iLevel][var] = smoothened_table_data[:, ivar]
        self.__printmsg("Done")
        return
    
    def setTableAxes(self, plane_cv_names:list[str], level_cv_name:str=None):
        """Select which controlling variables span the table. The two plane controlling variables are the
        coordinates of the mesh generated for each table level, the level controlling variable sweeps across
        table levels. Omitting the level controlling variable yields a single, two-dimensional table, even when
        the manifold is described by more controlling variables; the remaining controlling variables are then
        interpolated onto the table as ordinary state variables.

        By default the table is spanned by the controlling variables of the configuration in their own order.

        :param plane_cv_names: the two controlling variables spanning each table level.
        :type plane_cv_names: list[str]
        :param level_cv_name: controlling variable sweeping across table levels, defaults to None
        :type level_cv_name: str, optional
        :raises Exception: if other than two plane controlling variables are specified.
        :raises Exception: if a specified name is not a controlling variable of the configuration.
        :raises Exception: if a controlling variable is used for more than one table axis.
        """
        controlling_variables = self._Config.GetControllingVariables()

        if len(plane_cv_names) != 2:
            raise Exception("A table level is spanned by exactly two controlling variables.")

        table_cv_names = list(plane_cv_names)
        if level_cv_name is not None:
            table_cv_names.append(level_cv_name)

        for cv in table_cv_names:
            if cv not in controlling_variables:
                raise Exception("%s is not a controlling variable of the manifold (%s)" % (cv, ", ".join(controlling_variables)))
        if len(set(table_cv_names)) != len(table_cv_names):
            raise Exception("Each table axis should be spanned by a different controlling variable.")

        self._table_cv_names = table_cv_names
        self._nDim_table = len(table_cv_names)
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
    
    
    def meshTableLevel(self, level_index:int):
        """Discretize and interpolate the fluid data of a single table level.

        :param level_index: table level index
        :type level_index: int
        :return: list with table nodes, connectivity, perimiter indices, and interpolated fluid data.
        :rtype: list[np.ndarray[float]]
        """

        pointCloud = self._createPointCloudForTableLevel(self._table_levels[level_index])

        mesher = self._initiateMesher()
        if self.__run_parallel:
            mesher.setVerbosity(0)
        else:
            if self._verbosity > 1:
                mesher.setVerbosity(1)
            if self._verbosity > 2:
                mesher.setVerbosity(2)

        mesher.setInitialPointCloud(pointCloud)

        boundaryPolyline = self._createBoundaryPolylineForTableLevel(self._table_levels[level_index])
        if boundaryPolyline is not None:
            mesher.setBoundaryPolyline(boundaryPolyline)

        self._passRefinementOptions(mesher)

        mesher.generateMesh()

        cvTable, triangles, hullIDs, tableDataFrame = self.__extractTableLevelData(mesher)
        
        if self._is3D():
            self.__printmsg("Finished meshing table level %i at %s=%.2e with %i nodes" % (level_index, \
                                                                                    self._table_cv_names[2], \
                                                                                    self._table_levels[level_index], \
                                                                                    len(cvTable)))

        return [cvTable, triangles,hullIDs, tableDataFrame]
    
    def _createPointCloudForTableLevel(self, levelValue:float):
        """Generate a planar point cloud used as a reference for discretizing the table level.

        :param levelValue: value of the third table dimension corresponding to the table level.
        :type levelValue: float
        """
        return

    def _createBoundaryPolylineForTableLevel(self, levelValue:float):
        """Generate an ordered, closed polyline describing the perimiter of the table level. Returning None
        leaves the perimiter to be extracted as a hull around the reference point cloud.

        :param levelValue: value of the third table dimension corresponding to the table level.
        :type levelValue: float
        :return: planar coordinates of the perimiter, or None.
        :rtype: np.ndarray[float]
        """
        return None

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
        for iCv, cv in enumerate(self._table_cv_names):
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
        refinement_factor = 1.0
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
    
    def visualizeTableLevel(self, level_index:int=0, var_to_plot:str=None, plot_3d:bool=False, \
                            show_grid:bool=True, save_path:str=None, show:bool=False):
        """Visualize the mesh of a single table level and optionally colour it with a table variable.

        When the table data have not been generated yet, the requested level is meshed on the fly,
        which allows the mesh quality and refinement settings to be inspected before committing to
        the full table.

        :param level_index: index of the table level to visualize, defaults to 0
        :type level_index: int, optional
        :param var_to_plot: name of the table variable to colour the plot with. Showing the mesh alone
            is the default.
        :type var_to_plot: str, optional
        :param plot_3d: render the variable as a 3D surface rather than a 2D colour map, defaults to False
        :type plot_3d: bool, optional
        :param show_grid: overlay the triangulation on the plot, defaults to True
        :type show_grid: bool, optional
        :param save_path: file path to save the figure to, defaults to None
        :type save_path: str, optional
        :param show: keep the figure open for display, defaults to False
        :type show: bool, optional
        :raises Exception: if the table level index is out of range.
        :raises Exception: if the variable to plot is not part of the table data.
        """
        # The index is checked against the levels that are known, so that meshing is not
        # started for a level that does not exist.
        if self._N_table_levels is not None:
            self._checkTableLevelIndex(level_index)

        table_is_generated = (level_index < len(self._data_in_table)) and (self._data_in_table[level_index] is not None)
        if not table_is_generated:
            self._processTableLevels()
            if self._fluid_data_interpolator is None:
                self._defineFluidDataInterpolator()
            self._checkTableLevelIndex(level_index)

        if table_is_generated:
            table_nodes = self._table_nodes[level_index]
            connectivity = self._table_connectivity[level_index]
            hull_nodes = self._table_hullnodes[level_index]
            data_level = self._data_in_table[level_index]
        else:
            table_nodes, connectivity, hull_nodes, data_level = self.meshTableLevel(level_index)

        if var_to_plot is not None and var_to_plot not in data_level:
            raise Exception("%s is not part of the table data (%s)" % (var_to_plot, ", ".join(list(data_level.keys()))))

        x_nodes = table_nodes[:, 0]
        y_nodes = table_nodes[:, 1]
        if self._is3D():
            title = "%s = %.4e" % (self._table_cv_names[2], self._table_levels[level_index])
        else:
            title = "Table level %i" % level_index

        if var_to_plot is None:
            fig, ax = plt.subplots(figsize=(10, 10))
            ax.triplot(x_nodes, y_nodes, connectivity, linewidth=0.5)
            ax.plot(x_nodes[hull_nodes], y_nodes[hull_nodes], 'ko', ms=3, label="Perimiter nodes")
            ax.set_xlabel(self._table_cv_names[0], fontsize=14)
            ax.set_ylabel(self._table_cv_names[1], fontsize=14)
            ax.set_title(title, fontsize=14)
            ax.legend(fontsize=12)
        elif plot_3d:
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_trisurf(x_nodes, y_nodes, data_level[var_to_plot].to_numpy(dtype=float), \
                            triangles=connectivity, cmap='viridis', alpha=0.9, \
                            edgecolor='k' if show_grid else 'none', \
                            linewidth=0.2 if show_grid else 0)
            ax.set_xlabel(self._table_cv_names[0], fontsize=14)
            ax.set_ylabel(self._table_cv_names[1], fontsize=14)
            ax.set_zlabel(var_to_plot, fontsize=14)
            ax.set_title(title, fontsize=14)
        else:
            triangulation = mtri.Triangulation(x_nodes, y_nodes, connectivity)
            fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
            colormesh = ax.tripcolor(triangulation, data_level[var_to_plot].to_numpy(dtype=float), \
                                     shading='gouraud', cmap='inferno')
            if show_grid:
                ax.triplot(triangulation, color='k', linewidth=0.2, alpha=0.4)
            colorbar = fig.colorbar(colormesh, ax=ax, pad=0.02)
            colorbar.set_label(var_to_plot, fontsize=12)
            ax.set_xlabel(self._table_cv_names[0], fontsize=14)
            ax.set_ylabel(self._table_cv_names[1], fontsize=14)
            ax.set_title("%s at %s" % (var_to_plot, title), fontsize=14)

        self._saveOrCloseFigure(fig, save_path, show)
        return

    def plotTableSlices(self, x_cv:str, y_var:str, slice_cv:str, slice_range:tuple=None, n_slices:int=10, \
                        n_x_points:int=200, level_value:float=None, save_path:str=None, show:bool=False):
        """Plot a fluid state variable against one controlling variable for a number of fixed values of
        another controlling variable, evaluated with the fluid data interpolator. Such plots reveal
        whether the interpolated manifold is smooth and free of jumps before the table is written.

        :param x_cv: controlling variable along the horizontal axis.
        :type x_cv: str
        :param y_var: state variable along the vertical axis.
        :type y_var: str
        :param slice_cv: controlling variable kept constant along each line.
        :type slice_cv: str
        :param slice_range: lower and upper value of the slice controlling variable, defaults to its data extent
        :type slice_range: tuple, optional
        :param n_slices: number of slices distributed over the slice range, defaults to 10
        :type n_slices: int, optional
        :param n_x_points: number of sample points along each line, defaults to 200
        :type n_x_points: int, optional
        :param level_value: value of the remaining controlling variable of a three-dimensional table,
            defaults to the centre of its data range
        :type level_value: float, optional
        :param save_path: file path to save the figure to, defaults to None
        :type save_path: str, optional
        :param show: keep the figure open for display, defaults to False
        :type show: bool, optional
        :raises Exception: if the horizontal or slice variable does not span the table.
        :raises Exception: if the state variable is not available in the fluid data.
        """
        for cv in (x_cv, slice_cv):
            if cv not in self._table_cv_names:
                raise Exception("%s does not span the table (%s)" % (cv, ", ".join(self._table_cv_names)))
        if x_cv == slice_cv:
            raise Exception("The horizontal and slice controlling variable should be different.")

        if self._fluid_data_interpolator is None:
            self._defineFluidDataInterpolator()
        if y_var not in self._state_quantities:
            raise Exception("%s is not available in the fluid data" % y_var)

        x_index = self._table_cv_names.index(x_cv)
        slice_index = self._table_cv_names.index(slice_cv)
        y_index = self._state_quantities.index(y_var)

        # The scaler was fitted on the controlling variable data, so it carries their extent.
        cv_min = self._scaler_controlling_variables.data_min_
        cv_max = self._scaler_controlling_variables.data_max_

        if slice_range is None:
            slice_min, slice_max = cv_min[slice_index], cv_max[slice_index]
        else:
            slice_min, slice_max = float(slice_range[0]), float(slice_range[1])

        x_values = np.linspace(cv_min[x_index], cv_max[x_index], n_x_points)
        slice_values = np.linspace(slice_min, slice_max, n_slices)

        # Query points of a three-dimensional table are taken at a single value of the third
        # controlling variable, which is the one that is neither plotted nor sliced.
        cv_query = np.zeros([n_x_points, self._nDim_table])
        if self._is3D():
            third_index = ({0, 1, 2} - {x_index, slice_index}).pop()
            if level_value is None:
                level_value = 0.5 * (cv_min[third_index] + cv_max[third_index])
            cv_query[:, third_index] = level_value
        cv_query[:, x_index] = x_values

        # Query points outside the data support are blanked rather than extrapolated.
        cv_pointcloud_scaled = self._fluid_data_interpolator.getControllingVariableNodes()
        data_support = Delaunay(cv_pointcloud_scaled[:, [x_index, slice_index]])

        colormap = plt.cm.coolwarm
        color_norm = plt.Normalize(vmin=slice_min, vmax=slice_max)
        fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)

        for slice_value in slice_values:
            cv_query[:, slice_index] = slice_value
            cv_query_scaled = self._scaler_controlling_variables.transform(cv_query)

            within_support = data_support.find_simplex(cv_query_scaled[:, [x_index, slice_index]]) >= 0
            if not np.any(within_support):
                continue

            y_values = self._fluid_data_interpolator(cv_query_scaled)[:, y_index]
            y_values[~within_support] = np.nan

            ax.plot(x_values, y_values, color=colormap(color_norm(slice_value)), linewidth=1.2)

        scalar_map = plt.cm.ScalarMappable(cmap=colormap, norm=color_norm)
        scalar_map.set_array([])
        colorbar = fig.colorbar(scalar_map, ax=ax, pad=0.02)
        colorbar.set_label(slice_cv, fontsize=11)

        ax.set_xlabel(x_cv, fontsize=12)
        ax.set_ylabel(y_var, fontsize=12)
        ax.set_title("%s(%s) for %i slices of %s" % (y_var, x_cv, n_slices, slice_cv), fontsize=13)
        ax.grid(True, linestyle='--', alpha=0.4)

        self._saveOrCloseFigure(fig, save_path, show)
        return

    def _checkTableLevelIndex(self, level_index:int):
        """Verify that a table level index addresses an existing table level.

        :param level_index: table level index.
        :type level_index: int
        :raises Exception: if the index is out of range.
        """
        if (level_index < 0) or (level_index >= self._N_table_levels):
            raise Exception("Table level index should be between 0 and %i." % (self._N_table_levels - 1))
        return

    def _saveOrCloseFigure(self, fig, save_path:str, show:bool):
        """Save a figure when a file path is provided and close it unless it is to be displayed."""
        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight')
            self.__printmsg("Saved figure %s" % save_path)
        if not show:
            plt.close(fig)
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
    
    def getTableLevels(self):
        """Values of the level controlling variable of each table level. The values are prepared when the
        table is generated, so an empty array is returned before then.

        :return: table level values.
        :rtype: np.ndarray[float]
        """
        return np.asarray(self._table_levels, dtype=float)

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
    
    
    