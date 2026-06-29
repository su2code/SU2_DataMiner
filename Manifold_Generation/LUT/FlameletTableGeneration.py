###############################################################################################
#       #      _____ __  _____      ____        __        __  ____                   #        #
#       #     / ___// / / /__ \    / __ \____ _/ /_____ _/  |/  (_)___  ___  _____   #        #
#       #     \__ \/ / / /__/ /   / / / / __ `/ __/ __ `/ /|_/ / / __ \/ _ \/ ___/   #        #
#       #    ___/ / /_/ // __/   / /_/ / /_/ / /_/ /_/ / /  / / / / / /  __/ /       #        #
#       #   /____/\____//____/  /_____/\__,_/\__/\__,_/_/  /_/_/_/ /_/\___/_/        #        #
#       #                                                                            #        #
###############################################################################################

######################### FILE NAME: FlameletTableGenerator.py ################################
#=============================================================================================#
# author: Evert Bunschoten                                                                    |
#    :PhD Candidate ,                                                                         |
#    :Flight Power and Propulsion                                                             |
#    :TU Delft,                                                                               |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#   Table generator class for generating SU2-supported tables of flamelet data.               |
# Version: 3.1.0                                                                              |
#                                                                                             |
#=============================================================================================#

import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys,os
from Common.DataDrivenConfig import Config_FGM, Config
from Common.CommonMethods import GetReferenceData
from Common.Properties import DefaultSettings_FGM
import gmsh
import pickle
from multiprocessing import Pool
from Common.Interpolators import Invdisttree

class SU2TableGenerator_Base:
    _Config = None
    _savedir:str
    _table_variables:list[str] = None
    _manifold_variables:list[str]
    _controlling_variables:list[str] = None
    _manifold_data:np.ndarray[float] = None
    _manifold_data_interpolator:Invdisttree = None
    _base_cell_size:float = 1e-2#3.7e-3      # Table level base cell size.

    _refined_cell_size:float = 1e-3#2.5e-3#1.5e-3   # Table level refined cell size.
    _refinement_radius:float = 5e-3#5e-2     # Table level radius within which refinement is applied.
    _curvature_threshold:float = 0.3    # Curvature threshold above which refinement is applied.
    _n_near:int = 4     # Number of nearest neighbors from which to evaluate flamelet data.
    _p_fac:int = 5      # Power by which to weigh distances from query point.
    _control_var_scaler:MinMaxScaler =None
    _table_nodes = []       # Progress variable, total enthalpy, and mixture fraction node values for each table level.
    _table_nodes_norm = []  # Normalized table nodes for each level.
    _table_connectivity = []    # Table node connectivity per table level.
    _table_hullnodes = []   # Hull node indices per table level.

    def __init__(self, Config_in):
        self._Config = Config_in
        self._savedir = self._Config.GetOutputDir()
        return

    def SetSaveDir(self, save_dir:str):
        if not os.path.isdir(save_dir):
            raise Exception("Output directory %s not present on current hardware." % save_dir)
        self._savedir = save_dir
        return

    def SetBaseCellSize(self, cell_size:float):
        """
        Define the base cell size for the table levels.

        :param cell_size: Normalized coarse cell size for each 2D table mesh.
        :type cell_size: float
        :raise: Exception: if cell size is lower or equal to zero.
        """
        if cell_size > 0:
            self._base_cell_size = cell_size
        else:
            raise Exception("Proviced cell size should be higher than zero.")
        return

    def SetRefinedCellSize(self, cell_size:float):
        """
        Define the refinement cell size for the table levels.

        :param cell_size: Normalized fine cell size for each 2D table mesh.
        :type cell_size: float
        :raise: Exception: if cell size is lower or equal to zero.
        """
        if cell_size > 0:
            self._refined_cell_size = cell_size
        else:
            raise Exception("Proviced cell size should be higher than zero.")
        return

    def SetRefinementThreshold(self, val_threshold:float):
        """
        Define normalized curvature threshold beyond which refinement should be applied to each table level.

        :param val_threshold: Normalized curvature threshold value. All locations in the mesh with a higher curvature receive refinement.
        :type val_threshold: float
        :raises: Exception: If the threshold value is lower than zero.
        """

        if val_threshold > 0:
            self._curvature_threshold = val_threshold
        else:
            raise Exception("Curvature threshold value should be higher than zero.")
        return


    def DefineFlameletDataInterpolator(self):

        print("Configuring KD-tree for most accurate lookups")

        print("Loading flamelet data...")
        # Define scaler for FGM controlling variables.
        full_data_file = self._Config.GetOutputDir()+"/LUT_data_full.csv"
        with open(full_data_file,'r') as fid:
            self._manifold_variables = fid.readline().strip().split(',')
        #D_full = np.loadtxt(full_data_file,delimiter=',',skiprows=1)
        self._control_var_scaler = MinMaxScaler()
        CV_full, D_full = GetReferenceData(full_data_file, self._controlling_variables, self._manifold_variables)
        data_scaler = MinMaxScaler()
        data_scaler.fit_transform(D_full)

        CV_full_scaled = self._control_var_scaler.fit_transform(CV_full)

        # Exctract train and test data
        train_data_file = self._Config.GetOutputDir()+"/"+self._Config.GetConcatenationFileHeader()+"_train.csv"
        test_data_file = self._Config.GetOutputDir()+"/"+self._Config.GetConcatenationFileHeader()+"_test.csv"

        CV_train, D_train = GetReferenceData(train_data_file, self._controlling_variables, self._manifold_variables)
        CV_test, D_test = GetReferenceData(test_data_file, self._controlling_variables, self._manifold_variables)

        CV_train_scaled = self._control_var_scaler.transform(CV_train)
        CV_test_scaled = self._control_var_scaler.transform(CV_test)
        D_train_scaled = data_scaler.transform(D_train)
        D_test_scaled = data_scaler.transform(D_test)

        print("Done!")
        print("Setting up KD-tree...")
        self._lookup_tree = Invdisttree(X=CV_train_scaled,z=D_train_scaled)
        print("Done!")

        print("Search for best tree parameters...")
        # Do brute-force search to get the optimum number of nearest neighbors and distance power.
        n_near_range = range(1, 20)
        p_range = range(1, 6)
        RMS_ppv = np.zeros([len(n_near_range), len(p_range)])
        for i in tqdm(range(len(n_near_range))):
            for j in range(len(p_range)):
                PPV_predicted = self._lookup_tree(q=CV_test_scaled, nnear=n_near_range[i], p=p_range[j])
                rms_local = np.average(np.power(PPV_predicted - D_test_scaled, 2))
                RMS_ppv[i,j] = rms_local
        [imin,jmin] = divmod(RMS_ppv.argmin(), RMS_ppv.shape[1])
        self._n_near = n_near_range[imin]
        self._p_fac = p_range[jmin]
        print("Done!")
        print("Best found number of nearest neighbors: "+str(self._n_near))
        print("Best found distance power: "+str(self._p_fac))
        print("Setting up KD-tree...")
        self._lookup_tree = Invdisttree(X=CV_full_scaled,z=D_full)
        print("Done!")

        return

    def EvaluateManifoldInterpolator(self, CV_unscaled:np.ndarray):
        CV_scaled = self._control_var_scaler.transform(CV_unscaled)
        data_interp = self._lookup_tree(q=CV_scaled,nnear=self._n_near,p=self._p_fac)
        return data_interp

    def Compute2DTable(self, CV_1:str, CV_2:str):

        return

class SU2TableGenerator:

    _Config:Config_FGM = None # Config_FGM class from which to read settings.
    _savedir:str

    _mixfrac_min:float = None     # Minimum mixture fraction value of the flamelet data.
    _mixfrac_max:float = None     # Maximum mixture fraction value of the flamelet data.

    _pv_full_norm:np.ndarray[float] = None        # Normalized progress variable values of the flamelet data.
    _enth_full_norm:np.ndarray[float] = None      # Normalized total enthalpy values of the flamelet data.
    _mixfrac_full_norm:np.ndarray[float] = None   # Normalized mixture fraction values of the flamelet data.

    _Flamelet_Variables:list[str] = None  # Variable names in the concatenated flamelet data file.
    _Flamelet_Data:np.ndarray[float] = None     # Concatenated flamelet data.

    _custom_table_limits_set:bool = False
    _level_min_table:float = None     # Lower limit of the level (sweep) controlling variable.
    _level_max_table:float = None     # Upper limit of the level (sweep) controlling variable.
    _is_2D_table:bool = False           # True when a single-level (2D plane) table is requested.

    __run_parallel:bool = False
    __Np_cores:int = 1

    _N_table_levels:int = 100   # Number of table levels.
    _level_range_table:np.ndarray[float] = None   # Level-variable values at each table level.
    _base_cell_size:float = 2e-2#3.7e-3      # Table level base cell size.

    _refined_cell_size:float = 1e-2#2.5e-3#1.5e-3   # Table level refined cell size.
    _convex_hull_cell_size:float = 1e-2   # Cell size along the convex hull boundary.
    _refinement_radius:float = 2e-2#5e-2     # Table level radius within which refinement is applied.
    _curvature_threshold:float = 0.05    # Normalized gradient/curvature threshold above which refinement is applied (0.5 = top 50% of max gradient).
    _refinement_method:str = "gradient"  # Refinement method: "gradient" or "curvature".
    _refinement_fields:list[str] = ["ProdRateTot_PV"]  # Flamelet fields used as refinement criteria; max indicator across all fields drives seed selection.
    _max_refinement_seeds:int = 500  # Maximum number of seed points passed to gmsh after subsampling.
    _curvature_grid_resolution:int = 800  # Grid resolution (NxN) for computing curvature/gradient indicators.

    _delta_z_bin:float = 0.01  # Z-bin half-width for filtering nearest neighbors to same-Z region.
    _nn_filtered_points:int = 5  # Minimum number of nearest neighbors to keep after Z-direction filtering.
    _nn_edgesize_scaling:float = 5.0  # Scaling factor: spatial filter threshold = edge_size * nn_edgesize_scaling.
    _z_filter_factor:float = 0.5  # Factor for Z-direction filter: Z_threshold = z_filter_factor * delta_z_bin.

    _table_nodes = []       # Progress variable, total enthalpy, and mixture fraction node values for each table level.
    _table_nodes_norm = []  # Normalized table nodes for each level.
    _table_connectivity = []    # Table node connectivity per table level.
    _table_hullnodes = []   # Hull node indices per table level.
    __table_insert_levels:list[float] = []

    _controlling_variables:list[str]=[DefaultSettings_FGM.name_pv,\
                                      DefaultSettings_FGM.name_enth,\
                                      DefaultSettings_FGM.name_mixfrac]  # FGM controlling variables
    _lookup_tree:Invdisttree = None     # KD tree with inverse distance weighted interpolation for flamelet data interpolation.
    _flamelet_data_scaler:MinMaxScaler = None   # Scaler for flamelet data controlling variables.
    _scaler_2d:MinMaxScaler = None      # 2D (PV, H) scaler used when _is_2D_table is True.
    _D_full:np.ndarray = None           # Full flamelet data matrix, kept for potential 2D tree rebuild.
    _level_cv_name:str = DefaultSettings_FGM.name_mixfrac  # Controlling variable swept across table levels.
    _plane_cv_names:list[str] = [DefaultSettings_FGM.name_pv, DefaultSettings_FGM.name_enth]  # In-plane CVs defining the 2D mesh axes at each level.
    _level_cv_idx:int = 2     # Column index of the level CV; derived from _level_cv_name.
    _plane_cv_idxs:list[int] = [0, 1]  # Column indices of the two in-plane CVs; derived from _plane_cv_names.
    _n_near:int = 14     # Number of nearest neighbors from which to evaluate flamelet data.
    _p_fac:int = 3      # Power by which to weigh distances from query point.
    _custom_KDtreeparams:bool = False

    _preprocessed:bool = False

    def __init__(self, Config:Config_FGM, load_file:str=None, n_near:int=None, p_fac:int=None):
        """
        Initiate table generator class.

        :param Config: Config_FGM object.
        :type Config: Config_FGM
        """

        if n_near and p_fac:
            self._custom_KDtreeparams = True
            self._n_near = n_near
            self._p_fac = p_fac

        if load_file:
            # Load an existing TableGenerator object.
            with open(load_file, "rb") as fid:
                loaded_table_generator = pickle.load(fid)
            self.__dict__ = loaded_table_generator.__dict__.copy()
        else:
            # Create new TableGenerator object.
            self._Config = Config

            self.__DefineFlameletDataInterpolator()

        self._savedir = self._Config.GetOutputDir()
        return

    def SetSaveDir(self, save_dir:str):
        if not os.path.isdir(save_dir):
            raise Exception("Output directory %s not present on current hardware." % save_dir)
        self._savedir = save_dir

    def SetNTableLevels(self, N_levels:int):
        """
        Define the number of table levels in the mixture fraction direction.

        :param N_levels: number of table levels.
        :type N_levels: int
        :raise: Exception: if number of levels is lower than 1
        """
        if N_levels >= 1:
            self._N_table_levels = N_levels
        else:
            raise Exception("Number of table levels should be at least 1.")
        return

    def SetBaseCellSize(self, cell_size:float):
        """
        Define the base cell size for the table levels.

        :param cell_size: Normalized coarse cell size for each 2D table mesh.
        :type cell_size: float
        :raise: Exception: if cell size is lower or equal to zero.
        """
        if cell_size > 0:
            self._base_cell_size = cell_size
        else:
            raise Exception("Proviced cell size should be higher than zero.")
        return

    def SetRefinedCellSize(self, cell_size:float):
        """
        Define the refinement cell size for the table levels.

        :param cell_size: Normalized fine cell size for each 2D table mesh.
        :type cell_size: float
        :raise: Exception: if cell size is lower or equal to zero.
        """
        if cell_size > 0:
            self._refined_cell_size = cell_size
        else:
            raise Exception("Proviced cell size should be higher than zero.")
        return

    def SetRefinementThreshold(self, val_threshold:float):
        """
        Define normalized curvature threshold beyond which refinement should be applied to each table level.

        :param val_threshold: Normalized curvature threshold value. All locations in the mesh with a higher curvature receive refinement.
        :type val_threshold: float
        :raises: Exception: If the threshold value is lower than zero.
        """

        if val_threshold > 0:
            self._curvature_threshold = val_threshold
        else:
            raise Exception("Curvature threshold value should be higher than zero.")
        return

    def SetRefinementFields(self, field_names:list[str]):
        """Define the flamelet data fields used as refinement criteria.
        The maximum normalized gradient/curvature across all fields drives seed point selection.
        Field names are validated against the loaded flamelet data when the table is generated.

        :param field_names: list of variable names present in the flamelet data file.
        :type field_names: list[str]
        :raises Exception: if the list is empty or a field is not a string.
        """
        if len(field_names) == 0:
            raise Exception("At least one refinement field must be specified.")
        non_str = [f for f in field_names if not isinstance(f, str)]
        if non_str:
            raise Exception("All refinement field names must be strings; got: %s" % non_str)
        if self._Flamelet_Variables is not None:
            missing = [f for f in field_names if f not in self._Flamelet_Variables]
            if missing:
                raise Exception("Refinement field(s) not found in flamelet data: %s. "
                                "Available fields: %s" % (missing, self._Flamelet_Variables))
        self._refinement_fields = list(field_names)
        return

    def SetRefinementRadius(self, radius:float):
        """Define the radius around each refinement seed point within which the fine cell size is applied.

        :param radius: refinement radius in normalized coordinates.
        :type radius: float
        :raises Exception: if radius is not positive.
        """
        if radius <= 0:
            raise Exception("Refinement radius should be positive.")
        self._refinement_radius = radius
        return

    def SetRefinementMethod(self, method:str):
        """Select the indicator used to locate refinement seed points.

        :param method: ``'gradient'`` (default) or ``'curvature'``.
        :type method: str
        :raises Exception: if an unknown method is supplied.
        """
        if method not in ("gradient", "curvature"):
            raise Exception("Refinement method must be 'gradient' or 'curvature'.")
        self._refinement_method = method
        return

    def SetMaxRefinementSeeds(self, n_max:int):
        """Set the maximum number of seed points subsampled before passing to gmsh.
        Higher values give more accurate seed coverage at the cost of meshing time.

        :param n_max: maximum number of seed points.
        :type n_max: int
        :raises Exception: if n_max is lower than one.
        """
        if n_max < 1:
            raise Exception("Maximum number of refinement seeds should be at least one.")
        self._max_refinement_seeds = n_max
        return

    def SetHullCellSize(self, cell_size:float):
        """Define the cell size applied along the convex hull boundary (PV-min and PV-max lines).

        :param cell_size: normalized hull cell size.
        :type cell_size: float
        :raises Exception: if cell size is not positive.
        """
        if cell_size <= 0:
            raise Exception("Hull cell size should be positive.")
        self._convex_hull_cell_size = cell_size
        return

    def SetCurvatureGridResolution(self, n_points:int):
        """Set the resolution of the uniform grid used for computing curvature/gradient indicators.
        The grid is NxN points. Lower values speed up table generation significantly.
        Recommended: 200-400 for coarse, 400-600 for moderate, 800+ for fine.

        :param n_points: number of grid points per dimension (grid will be n_points x n_points).
        :type n_points: int
        :raises Exception: if n_points is less than 50.
        """
        if n_points < 50:
            raise Exception("Grid resolution should be at least 50 points per dimension.")
        self._curvature_grid_resolution = n_points
        return

    def SetZBinWidth(self, delta_z:float):
        """Set the Z-bin half-width for Z-aware nearest neighbor filtering.

        Reduces extrapolation by restricting neighbors to points within Z ± delta_z.
        Helps improve table quality in regions with sparse same-Z data.

        :param delta_z: Z-bin half-width (normalized [0, 1] or physical units).
        :type delta_z: float
        :raises Exception: if delta_z is not positive.
        """
        if delta_z <= 0:
            raise Exception("Z-bin width should be positive.")
        self._delta_z_bin = delta_z
        return

    def SetNNFilteredPoints(self, n_keep:int):
        """Set the minimum number of nearest neighbors to keep after Z-direction filtering.

        :param n_keep: minimum neighbors to keep after Z-filter (default 5).
        :type n_keep: int
        :raises Exception: if n_keep is less than 1.
        """
        if n_keep < 1:
            raise Exception("Minimum filtered neighbors should be at least 1.")
        self._nn_filtered_points = n_keep
        return

    def SetNNEdgeSizeScaling(self, scaling:float):
        """Set the scaling factor for spatial filtering thresholds.

        Spatial threshold = edge_size × scaling. Higher values relax spatial filter
        (use more distant neighbors), lower values tighten it (require closer neighbors).

        :param scaling: scaling factor (default 5.0).
        :type scaling: float
        :raises Exception: if scaling is not positive.
        """
        if scaling <= 0:
            raise Exception("NN edge size scaling should be positive.")
        self._nn_edgesize_scaling = scaling
        return

    def SetZFilterFactor(self, factor:float):
        """Set the reduction factor for Z-direction filtering threshold.

        Z_threshold = factor × delta_z_bin. Controls how aggressively to filter
        neighbors in Z-direction (0.5 = allow neighbors within half the Z-bin).

        :param factor: reduction factor (default 0.5, range typically 0.1 to 1.0).
        :type factor: float
        :raises Exception: if factor is not in (0, 1].
        """
        if factor <= 0 or factor > 1.0:
            raise Exception("Z filter factor should be in (0, 1].")
        self._z_filter_factor = factor
        return

    def SetTableAxes(self, level_cv_name:str, plane_cv_names:list[str]):
        """Define which controlling variable sweeps across table levels (the level axis)
        and which two span the 2D mesh at each level (the plane axes).

        Default: level = MixtureFraction, plane = [ProgressVariable, EnthalpyTot].

        :param level_cv_name: name of the controlling variable used as the sweep (level) axis.
        :type level_cv_name: str
        :param plane_cv_names: names of the two in-plane controlling variables.
        :type plane_cv_names: list[str]
        :raises Exception: if names are not among the known controlling variables, or if the
            wrong number of plane names is given.
        """
        if len(plane_cv_names) != 2:
            raise Exception("Exactly two plane CV names must be specified.")
        self._level_cv_name = level_cv_name
        self._plane_cv_names = list(plane_cv_names)
        self._resolve_cv_indices()
        if self._is_2D_table:
            self.__Rebuild2DInterpolator()
        return

    def _resolve_cv_indices(self) -> None:
        """Derive _level_cv_idx and _plane_cv_idxs from _level_cv_name and _plane_cv_names.
        Call this once after any change to the CV name attributes."""
        if self._level_cv_name not in self._controlling_variables:
            raise Exception("Level CV '%s' not found in controlling variables: %s"
                            % (self._level_cv_name, self._controlling_variables))
        for name in self._plane_cv_names:
            if name not in self._controlling_variables:
                raise Exception("Plane CV '%s' not found in controlling variables: %s"
                                % (name, self._controlling_variables))
        self._level_cv_idx  = self._controlling_variables.index(self._level_cv_name)
        self._plane_cv_idxs = [self._controlling_variables.index(n) for n in self._plane_cv_names]
        return

    def SetMixtureFractionLimits(self, mix_frac_min:float, mix_frac_max:float):
        """
        Define the mixture fraction limits of the table.

        :param mix_frac_min: Lower mixture fraction limit.
        :type mix_frac_min: float
        :param mix_frac_max: Upper mixture fraction limit.
        :type mix_frac_max: float
        :raise: Exception: If the upper mixture fraction limit is below the lower mixture fraction limit.
        """

        self._level_min_table = mix_frac_min
        self._level_max_table = mix_frac_max
        self._is_2D_table = (mix_frac_min == mix_frac_max)
        if self._is_2D_table:
            self._N_table_levels = 1
            self.__Rebuild2DInterpolator()
        self.__PrepareTableLevels()
        return

    def SetEquivalenceRatioLimits(self, phi_min:float, phi_max:float):
        """
        Define the table extent using equivalence ratio limits. The equivalence ratios
        are converted to mixture fractions using the fuel/oxidizer definition from the
        configuration. When phi_min == phi_max a single-Z 2D (PV, H) table is generated.

        :param phi_min: Lower equivalence ratio limit.
        :type phi_min: float
        :param phi_max: Upper equivalence ratio limit.
        :type phi_max: float
        """
        fuel_str = self._Config.GetFuelString()
        ox_str = self._Config.GetOxidizerString()
        self._Config.gas.set_equivalence_ratio(phi_min, fuel_str, ox_str)
        z_min = self._Config.gas.mixture_fraction(fuel_str, ox_str)
        self._Config.gas.set_equivalence_ratio(phi_max, fuel_str, ox_str)
        z_max = self._Config.gas.mixture_fraction(fuel_str, ox_str)
        self.SetMixtureFractionLimits(mix_frac_min=z_min, mix_frac_max=z_max)
        return

    def InsertMixtureFractionLevel(self, val_mixfrac_level:float):
        self.__table_insert_levels.append(val_mixfrac_level)
        self.__PrepareTableLevels()

    def __PrepareTableLevels(self):
        self._level_range_table = np.linspace(self._level_min_table, self._level_max_table, self._N_table_levels-len(self.__table_insert_levels))
        for z in self.__table_insert_levels:
            self._level_range_table = np.append(self._level_range_table, z)
        self._level_range_table = np.unique(np.sort(self._level_range_table))
        self._N_table_levels = len(self._level_range_table)
        return

    def SetNCores(self, n_cores:int):
        """Set the number of cores and enable parallel computing of the table level connectivity generation.

        :param n_cores: number of cores to distribute tasks over.
        :type n_cores: int
        :raises Exception: if the number of cores is lower than one.
        """
        if n_cores < 1:
            raise Exception("Number of cores should be at least one.")
        self.__Np_cores = n_cores
        self.__run_parallel = True
        return

    def __DefineFlameletDataInterpolator(self):

        print("Configuring KD-tree for most accurate lookups")

        print("Loading flamelet data...")
        # Define scaler for FGM controlling variables.
        full_data_file = self._Config.GetOutputDir()+"/"+self._Config.GetConcatenationFileHeader()+"_full.csv"
        with open(full_data_file,'r') as fid:
            self._Flamelet_Variables = fid.readline().strip().split(',')
        D_full = np.loadtxt(full_data_file,delimiter=',',skiprows=1)
        self._scaler = MinMaxScaler()
        CV_full = D_full[:,:3]
        self.__min_CV, self.__max_CV = np.min(CV_full,axis=0), np.max(CV_full,axis=0)

        self._resolve_cv_indices()
        min_level_dataset = self.__min_CV[self._level_cv_idx]
        max_level_dataset = self.__max_CV[self._level_cv_idx]

        self._level_min_table = min_level_dataset + 0.1*(max_level_dataset - min_level_dataset)
        self._level_max_table = max_level_dataset - 0.1*(max_level_dataset - min_level_dataset)

        CV_full_scaled = self._scaler.fit_transform(CV_full)

        # Exctract train and test data
        train_data_file = self._Config.GetOutputDir()+"/"+self._Config.GetConcatenationFileHeader()+"_train.csv"
        test_data_file = self._Config.GetOutputDir()+"/"+self._Config.GetConcatenationFileHeader()+"_test.csv"

        var_to_test_for = "ProdRateTot_PV"

        D_train = np.loadtxt(train_data_file,delimiter=',',skiprows=1)
        D_test = np.loadtxt(test_data_file,delimiter=',',skiprows=1)

        CV_train = np.vstack(tuple(D_train[:, self._Flamelet_Variables.index(c)] for c in self._controlling_variables)).T
        CV_test = np.vstack(tuple(D_test[:, self._Flamelet_Variables.index(c)] for c in self._controlling_variables)).T

        CV_train_scaled = self._scaler.transform(CV_train)
        CV_test_scaled = self._scaler.transform(CV_test)

        PPV_test = D_test[:, self._Flamelet_Variables.index(var_to_test_for)]
        print("Done!")
        print("Setting up KD-tree...")
        self._lookup_tree = Invdisttree(X=CV_train_scaled,z=D_train)
        print("Done!")

        if not self._custom_KDtreeparams:
            print("Search for best tree parameters...")
            # Do brute-force search to get the optimum number of nearest neighbors and distance power.
            n_near_range = range(1, 25)
            p_range = range(1, 6)
            RMS_ppv = np.zeros([len(n_near_range), len(p_range)])
            for i in tqdm(range(len(n_near_range))):
                for j in range(len(p_range)):
                    PPV_predicted = self._lookup_tree(q=CV_test_scaled, nnear=n_near_range[i], p=p_range[j])[:, self._Flamelet_Variables.index(var_to_test_for)]
                    rms_local = np.average(np.power(PPV_predicted - PPV_test, 2))
                    RMS_ppv[i,j] = rms_local
            [imin,jmin] = divmod(RMS_ppv.argmin(), RMS_ppv.shape[1])
            self._n_near = n_near_range[imin]
            self._p_fac = p_range[jmin]
            print("Done!")
        print("Best found number of nearest neighbors: "+str(self._n_near))
        print("Best found distance power: "+str(self._p_fac))
        print("Setting up KD-tree...")
        self._D_full = D_full
        self._lookup_tree = Invdisttree(X=CV_full_scaled,z=D_full)
        print("Done!")
        return

    def __Rebuild2DInterpolator(self):
        """Rebuild the KD-tree using only (PV, H) coordinates for a single-Z table.

        When a dataset has only one equivalence ratio/mixture fraction, the Z column
        in the full 3D MinMax-scaled space has a range that reflects differential-diffusion
        variation along the flame (not a true Z sweep).  That tiny real-world Δ is
        stretched to [0, 1] by the scaler, making Z-distances dominate the neighbour
        search and causing unphysical jumps.  Dropping Z from the search space fixes this.

        However, if the dataset contains MULTIPLE equivalence ratios (e.g., phi=0.6, 0.7, 0.8),
        we should keep the 3D interpolator even when generating a single-level table, because
        the 3D distance properly weights neighbors by their Z-distance.
        """
        if self._D_full is None:
            return  # interpolator not yet built

        # Check if the dataset has a true Z-range (multiple equivalence ratios)
        # vs differential-diffusion variation only (single equivalence ratio)
        z_values = self._D_full[:, self._level_cv_idx]
        z_range = np.max(z_values) - np.min(z_values)
        z_mean = np.mean(z_values)
        relative_z_range = z_range / z_mean if z_mean > 1e-10 else z_range

        # Threshold: if Z varies by more than 5%, assume multiple equivalence ratios
        if relative_z_range > 0.05:
            print(f"Dataset contains multiple equivalence ratios (Z-range: {z_range:.4f}, {relative_z_range*100:.1f}% of mean).")
            print("Keeping 3D interpolator for better Z-structure preservation.")
            self._scaler_2d = None  # Clear 2D scaler flag
            return  # Keep 3D interpolator

        print(f"Rebuilding 2D ({', '.join(self._plane_cv_names)}) KD-tree for single-equivalence-ratio dataset...")
        CV_full_2d = self._D_full[:, self._plane_cv_idxs]
        self._scaler_2d = MinMaxScaler()
        CV_full_2d_scaled = self._scaler_2d.fit_transform(CV_full_2d)
        self._lookup_tree = Invdisttree(X=CV_full_2d_scaled, z=self._D_full)
        print("Done!")
        return

    def __EvaluateFlameletInterpolator(self, CV_unscaled:np.ndarray,
                                       apply_z_filtering:bool=False,
                                       mesh_nodes_norm:np.ndarray=None,
                                       mesh_simplices:np.ndarray=None):
        """
        Evaluate flamelet interpolator with optional Z-bin filtering for mesh nodes.

        :param CV_unscaled: Controlling variable values (PV, H, Z).
        :param apply_z_filtering: If True, apply 2-stage Z-bin filtering (Stage 3).
        :param mesh_nodes_norm: Normalized mesh nodes (for computing connected_edge_distance).
        :param mesh_simplices: Delaunay triangulation simplices (for computing connected_edge_distance).
        :return: Interpolated flamelet data.
        """
        if apply_z_filtering:
            print(f"[DEBUG] __EvaluateFlameletInterpolator() called with apply_z_filtering=True")
        from scipy.spatial import cKDTree

        if not apply_z_filtering:
            # Original behavior: no Z-bin filtering
            if self._is_2D_table and self._scaler_2d is not None:
                CV_scaled = self._scaler_2d.transform(CV_unscaled[:, self._plane_cv_idxs])
            else:
                CV_scaled = self._scaler.transform(CV_unscaled)
            data_interp = self._lookup_tree(q=CV_scaled, nnear=self._n_near, p=self._p_fac)
            return data_interp

        # ===== STAGE 3: Mesh Node Evaluation with Z-bin Filtering =====

        if mesh_nodes_norm is None or mesh_simplices is None:
            raise ValueError("mesh_nodes_norm and mesh_simplices required for Z-filtered evaluation")

        level_cv_idx = self._level_cv_idx
        p0_idx = self._plane_cv_idxs[0]
        p1_idx = self._plane_cv_idxs[1]

        # Extract Z values for Z-bin filtering
        z_values = self._D_full[:, level_cv_idx]
        z_target = CV_unscaled[0, level_cv_idx]  # Z value of first mesh node (should be constant for all)
        mask_z_bin = np.abs(z_values - z_target) <= self._delta_z_bin
        D_z_bin = self._D_full[mask_z_bin]

        if len(D_z_bin) == 0:
            raise ValueError(f"Empty Z-bin at Z={z_target} ± {self._delta_z_bin}. Increase delta_z_bin parameter.")

        print(f"  Z-bin filtering for mesh: {len(D_z_bin)} of {len(self._D_full)} points in Z ± {self._delta_z_bin}")

        # Visualization: Stage 3 - Z-bin used for mesh node evaluation
        try:
            vis_path_3 = os.path.join(self._savedir, "zbin_stage3_mesh_eval.png") if hasattr(self, '_savedir') else None
            print(f"  [DEBUG] Stage 3 visualization: savedir={self._savedir if hasattr(self, '_savedir') else 'N/A'}, path={vis_path_3}")
            self.__VisualizeFlameletData3D(D_z_bin, f"Stage 3: Z-bin for Mesh Evaluation (Z={z_target:.6f} ± {self._delta_z_bin:.6f})", save_path=vis_path_3)
        except Exception as e:
            print(f"  [ERROR] Stage 3 visualization failed: {e}")

        # Check if 2D scaler is available; if not (multi-Z), build temporary one from Z-bin data
        if self._scaler_2d is None:
            print(f"  2D scaler not available (multi-Z dataset). Building temporary 2D scaler from Z-bin data for filtering...")
            CV_2d_bin = D_z_bin[:, [p0_idx, p1_idx]]
            scaler_2d_temp = MinMaxScaler()
            CV_2d_bin_norm = scaler_2d_temp.fit_transform(CV_2d_bin)
        else:
            CV_2d_bin = D_z_bin[:, [p0_idx, p1_idx]]
            CV_2d_bin_norm = self._scaler_2d.transform(CV_2d_bin)
            scaler_2d_temp = self._scaler_2d

        # Build temporary 2D (PV, H) k-d tree from Z-bin data
        kdtree_2d = cKDTree(CV_2d_bin_norm)

        # Compute connected_edge_distance for each mesh node
        # connected_edge_distance = max euclidean distance to Delaunay neighbors
        n_mesh_nodes = len(mesh_nodes_norm)
        connected_edge_distances = np.zeros(n_mesh_nodes)

        for i_node in range(n_mesh_nodes):
            # Find which simplices contain this node
            simplex_idxs = np.any(mesh_simplices == i_node, axis=1)
            neighbors = np.unique(mesh_simplices[simplex_idxs])
            neighbors = neighbors[neighbors != i_node]  # Remove self

            if len(neighbors) > 0:
                # Compute distances to all connected neighbors
                dist_to_neighbors = np.linalg.norm(
                    mesh_nodes_norm[neighbors] - mesh_nodes_norm[i_node], axis=1
                )
                connected_edge_distances[i_node] = dist_to_neighbors.max()
            else:
                connected_edge_distances[i_node] = 1e-3  # Fallback for isolated nodes

        # For each mesh node, apply 2-stage filtering and interpolate
        Q_interp_list = []
        insufficient_count = 0
        for i, mesh_node in enumerate(mesh_nodes_norm):
            if (i + 1) % 1000 == 0:
                print(f"    Processed {i+1}/{len(mesh_nodes_norm)} mesh nodes")

            # For multi-Z datasets with temporary scaler: transform mesh_node coords using temp scaler
            if self._scaler_2d is None:
                mesh_node_2d_phys = mesh_node[[p0_idx, p1_idx]]
                mesh_node_2d_norm = scaler_2d_temp.transform([mesh_node_2d_phys])[0]
            else:
                # For single-Z, mesh_node is already in consistent scaler space
                mesh_node_2d_norm = mesh_node[[p0_idx, p1_idx]]

            # Query nearest neighbors in (PV, H) space (candidate pool)
            distances_2d, indices_20 = kdtree_2d.query(mesh_node_2d_norm, k=min(self._n_near, len(D_z_bin)))

            # Get Z distances for these neighbors
            z_neighbors = D_z_bin[indices_20, level_cv_idx]
            z_distances = np.abs(z_neighbors - z_target)

            # Apply 2-stage filtering with scaled connected_edge_distance as spatial threshold
            spatial_threshold = connected_edge_distances[i] * self._nn_edgesize_scaling
            filtered_idx, had_insufficient = self.__FilterNearestNeighbors_2Stage(
                indices_20, distances_2d, z_distances, spatial_threshold
            )
            if had_insufficient:
                insufficient_count += 1

            # IDW interpolate using filtered neighbors
            filtered_data = D_z_bin[filtered_idx]
            filtered_distances = distances_2d[filtered_idx]

            # Inverse distance weighting
            with np.errstate(divide='ignore', invalid='ignore'):
                weights = 1.0 / np.power(filtered_distances + 1e-16, self._p_fac)
                weights /= weights.sum()
            q_val = weights @ filtered_data

            Q_interp_list.append(q_val)

        print(f"  Mesh nodes without sufficient neighbors within scaled edge size: {insufficient_count}/{len(mesh_nodes_norm)} ({100*insufficient_count/len(mesh_nodes_norm):.1f}%)")

        Q_interp = np.array(Q_interp_list)
        return Q_interp

    def PlotTableSlices(self,
                        x_cv: str,
                        y_var: str,
                        slice_cv: str,
                        slice_range: tuple = None,
                        n_slices: int = 10,
                        n_x_points: int = 200,
                        save_path: str = None,
                        show: bool = False):
        """Plot a dependent variable against one controlling variable for a set of
        fixed slices of another controlling variable, querying the IDW interpolator.

        Example: T(Z) for 10 enthalpy levels between h_min and h_max.

        :param x_cv: Name of the controlling variable to use as the x-axis (e.g. 'MixtureFraction').
        :type x_cv: str
        :param y_var: Name of the variable to plot on the y-axis (e.g. 'Temperature').
        :type y_var: str
        :param slice_cv: Name of the controlling variable held fixed per line (e.g. 'EnthalpyTot').
        :type slice_cv: str
        :param slice_range: (min, max) values of slice_cv. If None, the data extent is used.
        :type slice_range: tuple, optional
        :param n_slices: Number of slice values uniformly distributed across slice_range.
        :type n_slices: int
        :param n_x_points: Number of x_cv sample points per slice.
        :type n_x_points: int
        :param save_path: If provided, save the figure to this path instead of displaying it.
        :type save_path: str, optional
        :raises Exception: if x_cv, slice_cv, or y_var are not found in the flamelet data.
        """
        for name in [x_cv, slice_cv]:
            if name not in self._controlling_variables:
                raise Exception("'%s' is not a controlling variable. Available: %s"
                                % (name, self._controlling_variables))
        if y_var not in self._Flamelet_Variables:
            raise Exception("'%s' not found in flamelet data. Available: %s"
                            % (y_var, self._Flamelet_Variables))

        x_idx     = self._controlling_variables.index(x_cv)
        slice_idx = self._controlling_variables.index(slice_cv)
        y_col     = self._Flamelet_Variables.index(y_var)

        # Determine x and slice extents from the loaded data.
        x_min     = float(np.min(self._D_full[:, x_idx]))
        x_max     = float(np.max(self._D_full[:, x_idx]))
        if slice_range is None:
            s_min = float(np.min(self._D_full[:, slice_idx]))
            s_max = float(np.max(self._D_full[:, slice_idx]))
        else:
            s_min, s_max = float(slice_range[0]), float(slice_range[1])

        slice_values = np.linspace(s_min, s_max, n_slices)
        x_values     = np.linspace(x_min, x_max, n_x_points)

        # Level CV value: 0.0 for 2D mode (PV = 0), or the table centre otherwise.
        n_cv = len(self._controlling_variables)
        level_val = 0.0 if self._is_2D_table else 0.5 * (
            float(np.min(self._D_full[:, self._level_cv_idx])) +
            float(np.max(self._D_full[:, self._level_cv_idx])))

        # Build a Delaunay hull of the actual 2D data so that query points outside
        # the data support (e.g. high h at Z=1 for counterflow flames) are masked.
        from scipy.spatial import Delaunay as _Delaunay
        if self._is_2D_table and self._scaler_2d is not None:
            ZH_dim  = self._D_full[:, [x_idx, slice_idx]]
            ZH_norm = self._scaler_2d.transform(ZH_dim)
            _data_hull = _Delaunay(ZH_norm)
            _has_hull  = True
            _x_scaler_range  = [float(self._D_full[:, x_idx].min()),
                                 float(self._D_full[:, x_idx].max())]
            _s_scaler_range  = [float(self._D_full[:, slice_idx].min()),
                                 float(self._D_full[:, slice_idx].max())]
        else:
            _has_hull = False

        cmap   = plt.cm.coolwarm
        norm   = plt.Normalize(vmin=s_min, vmax=s_max)

        fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)

        for s_val in slice_values:
            # Build query array: all CVs at their level value, then overwrite x and slice.
            CV_query = np.zeros([n_x_points, n_cv])
            CV_query[:, self._level_cv_idx] = level_val
            CV_query[:, x_idx]              = x_values
            CV_query[:, slice_idx]          = s_val

            # Mask x values that fall outside the convex hull of the data.
            if _has_hull:
                x_norm = (x_values - _x_scaler_range[0]) / (_x_scaler_range[1] - _x_scaler_range[0] + 1e-32)
                s_norm = (s_val    - _s_scaler_range[0]) / (_s_scaler_range[1] - _s_scaler_range[0] + 1e-32)
                query_norm = np.column_stack([x_norm, np.full(n_x_points, s_norm)])
                inside = _data_hull.find_simplex(query_norm) >= 0
            else:
                inside = np.ones(n_x_points, dtype=bool)

            if not np.any(inside):
                continue   # entire slice is outside data support — skip

            result = self.__EvaluateFlameletInterpolator(CV_query)
            y_vals = result[:, y_col]
            y_vals[~inside] = np.nan   # blank extrapolation regions

            ax.plot(x_values, y_vals, color=cmap(norm(s_val)), linewidth=1.2)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, pad=0.02)
        cb.set_label(slice_cv, fontsize=11)

        ax.set_xlabel(x_cv, fontsize=12)
        ax.set_ylabel(y_var, fontsize=12)
        ax.set_title("%s(%s)  |  %d slices of %s" % (y_var, x_cv, n_slices, slice_cv), fontsize=13)
        ax.grid(True, linestyle='--', alpha=0.4)

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight')
            print("  Saved: %s" % save_path)
        if not show:
            plt.close(fig)
        return

    def VisualizeTableLevel(self, val_mix_frac:float, var_to_plot:str=None,
                            plot_3d:bool=False, show_grid:bool=True,
                            save_path:str=None, show:bool=False):
        """Compute and visualize the table mesh and optionally a data field at a given level value.

        :param val_mix_frac: value of the level controlling variable for which to generate the table level.
        :type val_mix_frac: float
        :param var_to_plot: name of the variable to colour the plot with. If None, only the mesh is shown.
        :type var_to_plot: str, optional
        :param plot_3d: when True and var_to_plot is set, render a 3D surface; when False render a 2D colour map.
        :type plot_3d: bool
        :param show_grid: overlay the triangulation wireframe on the plot.
        :type show_grid: bool
        :param save_path: if provided, save the figure to this path instead of displaying it.
        :type save_path: str, optional
        """

        Tria, Nodes, HullIdx, level_data, XY_ref_dim = self.ComputeTableLevelMesh(val_mix_frac)
        print("Total mesh nodes: %i, refinement seed points: %i" % (len(Nodes), len(XY_ref_dim)))

        p0, p1 = self._plane_cv_idxs
        title = "%s = %s" % (self._level_cv_name, str(val_mix_frac))

        if var_to_plot is None:
            # Mesh-only plot (always 2D).
            fig, ax = plt.subplots(figsize=(10, 10))
            ax.triplot(Nodes[:, p0], Nodes[:, p1], Tria, linewidth=0.5)
            ax.plot(Nodes[HullIdx, p0], Nodes[HullIdx, p1], 'ko', ms=3, label="Hull nodes")
            ax.set_xlabel(self._plane_cv_names[0], fontsize=14)
            ax.set_ylabel(self._plane_cv_names[1], fontsize=14)
            ax.set_title(title, fontsize=14)
            ax.legend(fontsize=12)

        elif plot_3d:
            # 3D surface plot.
            var_idx = self._Flamelet_Variables.index(var_to_plot)
            fig = plt.figure(figsize=(10, 10))
            ax  = fig.add_subplot(111, projection='3d')
            ax.plot_trisurf(Nodes[:, p0], Nodes[:, p1], level_data[:, var_idx],
                            triangles=Tria, cmap='viridis', alpha=0.9,
                            edgecolor='k' if show_grid else 'none',
                            linewidth=0.2 if show_grid else 0)
            ax.set_xlabel(self._plane_cv_names[0], fontsize=14)
            ax.set_ylabel(self._plane_cv_names[1], fontsize=14)
            ax.set_zlabel(var_to_plot, fontsize=14)
            ax.set_title(title, fontsize=14)

        else:
            # 2D colour-map plot.
            import matplotlib.tri as mtri
            var_idx  = self._Flamelet_Variables.index(var_to_plot)
            values   = level_data[:, var_idx]
            triang   = mtri.Triangulation(Nodes[:, p0], Nodes[:, p1], Tria)
            fig, ax  = plt.subplots(figsize=(10, 7), constrained_layout=True)
            tc = ax.tripcolor(triang, values, shading='gouraud', cmap='inferno')
            if show_grid:
                ax.triplot(triang, color='k', linewidth=0.2, alpha=0.4)
            cb = fig.colorbar(tc, ax=ax, pad=0.02)
            cb.set_label(var_to_plot, fontsize=12)
            ax.set_xlabel(self._plane_cv_names[0], fontsize=14)
            ax.set_ylabel(self._plane_cv_names[1], fontsize=14)
            ax.set_title("%s — %s" % (var_to_plot, title), fontsize=14)

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight')
            print("  Saved: %s" % save_path)
        if not show:
            plt.close(fig)
        return

    def GenerateTableNodes(self):
        """
        Generate the table nodes and connectivity.
        """

        self.__PrepareTableLevels()

        self._table_nodes = [None] * self._N_table_levels
        self._table_nodes_norm = [None] * self._N_table_levels
        self._table_connectivity = [None] * self._N_table_levels
        self._table_hullnodes = [None] * self._N_table_levels
        self.table_data = [None] * self._N_table_levels

        flamelet_vars = []
        for var in self._Flamelet_Variables:
            flamelet_vars.append(var)
        for cv in self._controlling_variables:
            if cv in flamelet_vars:
                flamelet_vars.remove(cv)
        if "FlameletID" in flamelet_vars:
            flamelet_vars.remove("FlameletID")

        self.table_vars = flamelet_vars
        nVars = len(self.table_vars)

        # Generate the table cells for each table level.
        NTria = 0
        NHull = 0
        NNodes = 0
        if self.__run_parallel:
            pool = Pool(self.__Np_cores)
            results = pool.map(self.ComputeTableNodes, [i for i in range(self._N_table_levels)])
            pool.close()
            for iLevel in range(self._N_table_levels):
                self._table_nodes[iLevel] = results[iLevel][0]
                self._table_connectivity[iLevel] = results[iLevel][1]
                self._table_hullnodes[iLevel] = results[iLevel][2]
                data_interp = results[iLevel][3]
                table_data_level = [None] * nVars
                for iVar in range(nVars):
                    var = self.table_vars[iVar]
                    table_data_level[iVar] = data_interp[:, self._Flamelet_Variables.index(var)]
                self.table_data[iLevel] = table_data_level
                NTria += np.shape(self._table_connectivity[iLevel])[0]
                NHull += np.shape(self._table_hullnodes[iLevel])[0]
                NNodes += np.shape(self._table_nodes[iLevel])[0]
        else:
            for iLevel in range(self._N_table_levels):
                table_data_level = [None] * nVars
                result = self.ComputeTableNodes(iLevel)
                self._table_nodes[iLevel] = result[0]
                self._table_connectivity[iLevel] = result[1]
                self._table_hullnodes[iLevel] = result[2]
                data_interp = result[3]
                for iVar in range(nVars):
                    var = self.table_vars[iVar]
                    table_data_level[iVar] = data_interp[:, self._Flamelet_Variables.index(var)]
                self.table_data[iLevel] = table_data_level

                NTria += np.shape(self._table_connectivity[iLevel])[0]
                NHull += np.shape(self._table_hullnodes[iLevel])[0]
                NNodes += np.shape(self._table_nodes[iLevel])[0]

        NTria_average = int(NTria / self._N_table_levels)
        NHull_average = int(NHull / self._N_table_levels)
        NNodes_average = int(NNodes / self._N_table_levels)
        print("Average number of nodes: %i" % NNodes_average)
        print("Average number of elements: %i" % NTria_average)
        print("Average number of hull nodes: %i" % NHull_average)
        return

    def ComputeTableNodes(self, iLevel:int):
        """Compute the table connectivity for a specific table level.

        :param iLevel: table level index.
        :type iLevel: int
        :raises Exception: if the table level index is not between 0 and the number of table levels.
        :return: table nodes(dimensional), table nodes(normalized), connectivity, hull node indices
        :rtype: list[np.ndarray]
        """
        if iLevel < 0 or iLevel > self._N_table_levels:
            raise Exception("Specified table level out of bounds.")

        # Compute the connectivity, normalized node values, and hull indices for the table level at the respective
        #   mixture fraction value.
        level_val = self._level_range_table[iLevel]
        Tria, Nodes_dim, HullIdx, TableDataLevel, _ = self.ComputeTableLevelMesh(level_val)

        print("Computed triangulation on level %i out of %i with %i nodes." % (iLevel+1, self._N_table_levels, len(Nodes_dim)))

        return [Nodes_dim, Tria, HullIdx, TableDataLevel]


    def WriteTableFile(self, output_filepath:str=None):
        """
        Save the table data and connectivity as a Dragon library file. If no file name is provided,
        the table file will be named according to the Config_FGM class name.

        When _is_2D_table is True (mix_frac_min == mix_frac_max and N_table_levels == 1), a
        2D Dragon library (version 1.0.1) is written with ProgressVariable and EnthalpyTot as
        the only controlling variables and no <Level> section wrappers.

        :param output_filepath: optional output filepath for table file.
        :type output_filepath: str
        """

        if output_filepath:
            file_out = output_filepath
        else:
            file_out = self._savedir + "/LUT_"+self._Config.GetConfigName()+".drg"

        print("Writing LUT file with name " + file_out)
        fid = open(file_out, "w+")

        if self._is_2D_table:
            self.__WriteTableFile2D(fid)
        else:
            self.__WriteTableFile3D(fid)

        fid.close()
        return

    def __WriteTableFile2D(self, fid):
        """Write a 2D Dragon library file (version 1.0.1) for a single table level.
        The controlling variable names are taken from _plane_cv_names."""

        Nodes = self._table_nodes[0]          # shape (Np, n_cv)
        Connectivity = self._table_connectivity[0]
        HullNodes = self._table_hullnodes[0]
        Np = np.shape(Nodes)[0]

        cv0_name = self._plane_cv_names[0]
        cv1_name = self._plane_cv_names[1]

        fid.write("Dragon library\n\n")
        fid.write("<Header>\n\n")
        fid.write("[Version]\n1.0.1\n\n")
        fid.write("[Number of points]\n%i\n\n" % Np)
        fid.write("[Number of triangles]\n%i\n\n" % np.shape(Connectivity)[0])
        fid.write("[Number of hull points]\n%i\n\n" % np.shape(HullNodes)[0])

        if cv0_name == DefaultSettings_FGM.name_pv:
            fid.write("[Progress variable definition]\n")
            fid.write("+".join(("%+.4e * %s" % (w, s)) for w, s in zip(
                self._Config.GetProgressVariableWeights(),
                self._Config.GetProgressVariableSpecies())) + "\n\n")

        cv0_vals = Nodes[:, self._plane_cv_idxs[0]]
        cv1_vals = Nodes[:, self._plane_cv_idxs[1]]
        fid.write("[%s min]\n%e\n\n" % (cv0_name, np.min(cv0_vals)))
        fid.write("[%s max]\n%e\n\n" % (cv0_name, np.max(cv0_vals)))
        fid.write("[%s min]\n%e\n\n" % (cv1_name, np.min(cv1_vals)))
        fid.write("[%s max]\n%e\n\n" % (cv1_name, np.max(cv1_vals)))

        all_vars_2d = [cv0_name, cv1_name] + self.table_vars
        fid.write("[Number of variables]\n%i\n\n" % len(all_vars_2d))
        fid.write("[Variable names]\n")
        for var in all_vars_2d:
            fid.write(var + "\n")
        fid.write("\n")

        fid.write("</Header>\n\n")

        print("Writing table data...")
        fid.write("<Data>\n")
        for iNode in tqdm(range(Np)):
            fid.write("%+.14e %+.14e" % (cv0_vals[iNode], cv1_vals[iNode]))
            for iVar in range(len(self.table_vars)):
                fid.write(" %+.14e" % self.table_data[0][iVar][iNode])
            fid.write("\n")
        fid.write("</Data>\n\n")
        print("Done!")

        print("Writing table connectivity...")
        fid.write("<Connectivity>\n")
        for iCell in tqdm(range(len(Connectivity))):
            fid.write(" ".join("%i" % c for c in Connectivity[iCell, :] + 1) + "\n")
        fid.write("</Connectivity>\n\n")
        print("Done!")

        print("Writing hull nodes...")
        fid.write("<Hull>\n")
        for iCell in range(len(HullNodes)):
            fid.write("%i\n" % (HullNodes[iCell] + 1))
        fid.write("</Hull>\n")
        print("Done!")
        return

    def __WriteTableFile3D(self, fid):
        """Write a 3D multi-level Dragon library file (version 1.1.0)."""

        fid.write("Dragon library\n\n")
        fid.write("<Header>\n\n")
        fid.write("[Version]\n1.1.0\n\n")
        fid.write("[Progress variable definition]\n")
        fid.write("+".join(("%+.4e * %s" % (w, s)) for w, s in zip(self._Config.GetProgressVariableWeights(), self._Config.GetProgressVariableSpecies())) + "\n\n")

        fid.write("[Number of table levels]\n%i\n\n" % self._N_table_levels)
        fid.write("[Table levels]\n")
        for z in self._level_range_table:
            fid.write("%+.16e\n" % z)
        fid.write("\n")

        fid.write("[Number of points]\n")
        for Nodes in self._table_nodes:
            fid.write("%i\n" % np.shape(Nodes)[0])
        fid.write("\n")

        fid.write("[Number of triangles]\n")
        for Elements in self._table_connectivity:
            fid.write("%i\n" % np.shape(Elements)[0])
        fid.write("\n")

        fid.write("[Number of hull points]\n")
        for HullNodes in self._table_hullnodes:
            fid.write("%i\n" % np.shape(HullNodes)[0])
        fid.write("\n")

        fid.write("[Number of variables]\n%i\n\n" % len(self._Flamelet_Variables))
        fid.write("[Variable names]\n")
        for iVar, Var in enumerate(self._Flamelet_Variables):
            fid.write(str(iVar + 1)+":"+Var+"\n")
        fid.write("\n")

        fid.write("</Header>\n\n")

        print("Writing table data...")
        fid.write("<Data>\n")
        for iLevel in tqdm(range(len(self._table_nodes))):
            fid.write("<Level>\n")
            Np = np.shape(self._table_nodes[iLevel])[0]
            for iNode in range(Np):
                fid.write("\t".join("%+.14e" % cv for cv in self._table_nodes[iLevel][iNode, :]))
                for iVar in range(len(self.table_vars)):
                    fid.write("\t%+.14e" % self.table_data[iLevel][iVar][iNode])
                fid.write("\n")
            fid.write("</Level>\n")
        fid.write("</Data>\n\n")
        print("Done!")

        print("Writing table connectivity...")
        fid.write("<Connectivity>\n")
        for iLevel in tqdm(range(len(self._table_connectivity))):
            fid.write("<Level>\n")
            for iCell in range(len(self._table_connectivity[iLevel])):
                fid.write("\t".join("%i" % c for c in self._table_connectivity[iLevel][iCell, :]+1) + "\n")
            fid.write("</Level>\n")
        fid.write("</Connectivity>\n\n")
        print("Done!")

        print("Writing hull nodes...")
        fid.write("<Hull>\n")
        for iLevel in tqdm(range(len(self._table_hullnodes))):
            fid.write("<Level>\n")
            for iCell in range(len(self._table_hullnodes[iLevel])):
                fid.write(("%i" % (self._table_hullnodes[iLevel][iCell]+1)) + "\n")
            fid.write("</Level>\n")
        fid.write("</Hull>\n\n")
        print("Done!")
        return

    def ComputeTableLevelMesh(self, val_mix_frac:float):
        """
        Compute the table nodes, connectivity, and convex hull node indices of a 2D table level for a given level-variable value.

        :param val_mix_frac: Value of the level controlling variable for which to generate a 2D table level.
        :type val_mix_frac: float
        :return Connectivity: Delaunay triangulation connectivity
        :rtype Connecivity: NDarray
        :return MeshNodes:
        """
        val_level = val_mix_frac

        # For the (Z, h) plane case use the actual data hull to avoid extrapolation
        # outside the skewed parallelogram-shaped data support.
        # Handle both configurations: (Z,H) plane with PV level OR (PV,H) plane with Z level
        # BUT only if _scaler_2d is available (single-Z tables; multi-Z tables use regular path)
        is_zh_table = (self._level_cv_name == DefaultSettings_FGM.name_pv and
                       self._plane_cv_names[0] == DefaultSettings_FGM.name_mixfrac and
                       self._plane_cv_names[1] == DefaultSettings_FGM.name_enth)
        is_pvh_with_z_level = (self._level_cv_name == DefaultSettings_FGM.name_mixfrac and
                               self._plane_cv_names[0] == DefaultSettings_FGM.name_pv and
                               self._plane_cv_names[1] == DefaultSettings_FGM.name_enth)

        if (is_zh_table or is_pvh_with_z_level) and self._scaler_2d is not None:
            Coord_refinement, Coord_hull, hull_area, level_norm = \
                self.__ComputeCurvatureZH()
        else:
            Coord_refinement, Coord_hull, hull_area, level_norm, CV_mesh, table_level_data = \
                self.__ComputeCurvature(val_level)

        MeshNodes_Norm, table_level_data = self.__Compute2DMesh(XY_hull=Coord_hull, XY_refinement=Coord_refinement, val_level_norm=level_norm, level_area=hull_area)

        Tria = Delaunay(MeshNodes_Norm[:, self._plane_cv_idxs])
        HullNodes = Tria.convex_hull[:, 0]
        MeshNodes_dim = self._scaler.inverse_transform(MeshNodes_Norm)

        # Inverse-transform refinement seed coordinates for visualization.
        if len(Coord_refinement) > 0:
            n_cv = len(self._controlling_variables)
            XY_ref_3d = np.zeros([len(Coord_refinement), n_cv])
            XY_ref_3d[:, self._plane_cv_idxs[0]] = Coord_refinement[:, 0]
            XY_ref_3d[:, self._plane_cv_idxs[1]] = Coord_refinement[:, 1]
            XY_ref_3d[:, self._level_cv_idx]      = level_norm
            XY_refinement_dim = self._scaler.inverse_transform(XY_ref_3d)
        else:
            XY_refinement_dim = np.empty((0, 3))

        return Tria.simplices, MeshNodes_dim, HullNodes, table_level_data, XY_refinement_dim

    def __GetPlaneCVBounds(self, val_level:float):
        """Compute physical ranges for the two plane CVs at the given level value.

        Returns ``(plane0_unb, plane0_b, plane1_min, plane1_max, plane1_at_unb)`` where:

        * ``plane0_unb`` / ``plane0_b`` — plane-CV-0 at the unburnt and fully-burnt states.
        * ``plane1_min`` — minimum plane-CV-1 (cooled burnt state at T_lower).
        * ``plane1_max`` — maximum plane-CV-1 (reactant at T_upper).
        * ``plane1_at_unb`` — plane-CV-1 at the unburnt reactant state (T_lower),
          used to define the lower manifold boundary line.

        Currently implemented for ``level = MixtureFraction``,
        ``plane = [ProgressVariable, EnthalpyTot]``.
        Add additional ``elif`` branches to support other axis combinations.

        :param val_level: value of the level controlling variable.
        :type val_level: float
        :raises NotImplementedError: if the current axis combination is not yet supported.
        """
        if (self._level_cv_name == DefaultSettings_FGM.name_mixfrac and
                self._plane_cv_names[0] == DefaultSettings_FGM.name_pv and
                self._plane_cv_names[1] == DefaultSettings_FGM.name_enth):
            fuel = self._Config.GetFuelString()
            ox   = self._Config.GetOxidizerString()
            T_lo, T_hi = self._Config.GetUnbTempBounds()
            p = DefaultSettings_FGM.pressure

            self._Config.gas.set_mixture_fraction(val_level, fuel, ox)
            self._Config.gas.TP = T_lo, p
            plane1_at_unb = self._Config.gas.enthalpy_mass
            plane0_unb = self._Config.ComputeProgressVariable(
                variables=None, flamelet_data=None,
                Y_flamelet=self._Config.gas.Y[:, np.newaxis])[0]

            self._Config.gas.TP = T_hi, p
            plane1_max = self._Config.gas.enthalpy_mass

            self._Config.gas.set_mixture_fraction(val_level, fuel, ox)
            self._Config.gas.TP = T_lo, p
            self._Config.gas.equilibrate("HP")
            plane0_b = self._Config.ComputeProgressVariable(
                variables=None, flamelet_data=None,
                Y_flamelet=self._Config.gas.Y[:, np.newaxis])[0]
            self._Config.gas.TP = T_lo, p
            plane1_min = self._Config.gas.enthalpy_mass

            return plane0_unb, plane0_b, plane1_min, plane1_max, plane1_at_unb
        elif (self._level_cv_name == DefaultSettings_FGM.name_pv and
              self._plane_cv_names[0] == DefaultSettings_FGM.name_mixfrac and
              self._plane_cv_names[1] == DefaultSettings_FGM.name_enth):
            # (Z, h) table: derive bounds directly from the loaded data extent.
            p0_idx = self._plane_cv_idxs[0]  # MixtureFraction column index in _D_full
            p1_idx = self._plane_cv_idxs[1]  # EnthalpyTot column index in _D_full
            plane0_unb = float(np.min(self._D_full[:, p0_idx]))   # Z_min (≈ 0, pure oxidizer)
            plane0_b   = float(np.max(self._D_full[:, p0_idx]))   # Z_max (≈ 1, pure fuel)
            plane1_min = float(np.min(self._D_full[:, p1_idx]))   # h_min (coolest flame)
            plane1_max = float(np.max(self._D_full[:, p1_idx]))   # h_max (hottest flame)
            # Set plane1_at_unb == plane1_min so that h_limit = plane1_min everywhere
            # and the boundary-line filter in __ComputeCurvature keeps all grid points.
            plane1_at_unb = plane1_min
            return plane0_unb, plane0_b, plane1_min, plane1_max, plane1_at_unb
        else:
            raise NotImplementedError(
                "Plane CV bounds not implemented for level='%s', plane=%s. "
                "Currently only supported: level='%s', plane=['%s', '%s'] or "
                "level='%s', plane=['%s', '%s']." % (
                    self._level_cv_name, self._plane_cv_names,
                    DefaultSettings_FGM.name_mixfrac,
                    DefaultSettings_FGM.name_pv,
                    DefaultSettings_FGM.name_enth,
                    DefaultSettings_FGM.name_pv,
                    DefaultSettings_FGM.name_mixfrac,
                    DefaultSettings_FGM.name_enth))

    def __FilterNearestNeighbors_2Stage(self, neighbors_20_idx:np.ndarray, neighbor_distances_2d:np.ndarray,
                                         neighbor_z_distances:np.ndarray, spatial_threshold:float):
        """
        Apply 2-stage filtering to nearest neighbors: spatial then Z-direction.

        :param neighbors_20_idx: Indices of 20 nearest neighbors in the flamelet data.
        :param neighbor_distances_2d: Euclidean distances in (PV, H) space for each neighbor.
        :param neighbor_z_distances: Absolute Z-distance for each neighbor.
        :param spatial_threshold: Threshold distance in (PV, H) space (scaled edge size).
        :return: Tuple (filtered_neighbor_indices, had_insufficient_neighbors_flag)
        """
        # Sort by (PV, H) distance
        sorted_idx_spatial = np.argsort(neighbor_distances_2d)

        # FILTER STAGE 1 - Spatial: keep until (n_near - nn_filtered_points remain) OR (distance ≤ threshold)
        stage1_target = self._n_near - self._nn_filtered_points
        keep_spatial = []
        for idx in sorted_idx_spatial:
            keep_spatial.append(idx)
            if len(keep_spatial) >= stage1_target or neighbor_distances_2d[idx] <= spatial_threshold:
                break

        had_insufficient = False
        if len(keep_spatial) < stage1_target:
            # Collect all within threshold
            keep_spatial = [idx for idx in sorted_idx_spatial if neighbor_distances_2d[idx] <= spatial_threshold]
            if not keep_spatial:
                keep_spatial = list(sorted_idx_spatial[:stage1_target])
                had_insufficient = True

        # FILTER STAGE 2 - Z-direction: keep until (n_near - 2*nn_filtered_points remain) OR (|ΔZ| ≤ z_filter_factor×delta_z_bin)
        keep_spatial_arr = np.array(keep_spatial)
        z_threshold = self._z_filter_factor * self._delta_z_bin
        stage2_target = self._n_near - 2 * self._nn_filtered_points

        sorted_idx_z = keep_spatial_arr[np.argsort(neighbor_z_distances[keep_spatial_arr])]
        keep_z = []
        for idx in sorted_idx_z:
            keep_z.append(idx)
            if len(keep_z) >= stage2_target or neighbor_z_distances[idx] <= z_threshold:
                break

        if len(keep_z) < stage2_target:
            # Collect all within Z threshold
            keep_z = [idx for idx in sorted_idx_z if neighbor_z_distances[idx] <= z_threshold]
            if not keep_z:
                keep_z = list(sorted_idx_z[:stage2_target])

        return np.array(keep_z, dtype=int), had_insufficient


    def __ComputeCurvatureZH(self):
        """Hull and refinement seeds for a (Z, h) table using the actual data cloud.

        Unlike the rectangular-grid approach in __ComputeCurvature, this method
        derives the convex hull directly from the scattered (Z, h) data points
        in _D_full.  This correctly follows the skewed parallelogram shape of the
        counterflow flame data (h shifts by ~-4.7 MJ/kg from Z=0 to Z=1) and
        avoids generating mesh nodes outside the data support.

        Applies 2-stage Z-bin filtering (spatial + Z-direction) to nearest-neighbor
        queries to improve quality in sparse regions.

        Returns the same tuple as __ComputeCurvature but without CV_mesh and
        table_level_data (those are computed after meshing).
        """
        print(f"[DEBUG] __ComputeCurvatureZH() called for (Z, h) table")
        from Interpolators import Invdisttree
        from scipy.spatial import cKDTree

        p0_idx = self._plane_cv_idxs[0]   # MixtureFraction column
        p1_idx = self._plane_cv_idxs[1]   # EnthalpyTot column
        level_cv_idx = self._level_cv_idx

        # Normalized 2D data cloud (use the 2D scaler built by __Rebuild2DInterpolator).
        ZH_dim  = self._D_full[:, [p0_idx, p1_idx]]
        ZH_norm = self._scaler_2d.transform(ZH_dim)   # (N, 2) in [0,1]^2

        # Convex hull of the actual data → correct parallelogram boundary.
        data_hull  = ConvexHull(ZH_norm)
        hull_verts = ZH_norm[data_hull.vertices, :]    # (Nh, 2) normalized
        XY_hull    = hull_verts                        # passed to __Compute2DMesh

        # ===== STAGE 1: Refinement Field Evaluation on Analytical Grid =====

        # Extract Z values for Z-bin filtering
        z_values = self._D_full[:, level_cv_idx]
        z_target = 0.0  # For (Z, h) tables at Z=0 normalized
        mask_z_bin = np.abs(z_values - z_target) <= self._delta_z_bin
        D_z_bin = self._D_full[mask_z_bin]

        if len(D_z_bin) == 0:
            raise ValueError(f"Empty Z-bin at Z={z_target} ± {self._delta_z_bin}. Increase delta_z_bin parameter.")

        print(f"  Z-bin filtering: {len(D_z_bin)} of {len(self._D_full)} points in Z ± {self._delta_z_bin}")

        # Visualization: Stage 1 - Initial Z-bin data
        try:
            vis_path_1 = os.path.join(self._savedir, "zbin_stage1_initial.png") if hasattr(self, '_savedir') else None
            print(f"  [DEBUG] Stage 1 visualization: savedir={self._savedir if hasattr(self, '_savedir') else 'N/A'}, path={vis_path_1}")
            self.__VisualizeFlameletData3D(D_z_bin, f"Stage 1: Initial Z-bin Data (Z={z_target:.6f} ± {self._delta_z_bin:.6f})", save_path=vis_path_1)
        except Exception as e:
            print(f"  [ERROR] Stage 1 visualization failed: {e}")

        # Build temporary 2D (PV, H) k-d tree from Z-bin data
        CV_2d_bin = D_z_bin[:, [p0_idx, p1_idx]]
        CV_2d_bin_norm = self._scaler_2d.transform(CV_2d_bin)
        temp_tree_2d = Invdisttree(X=CV_2d_bin_norm, z=D_z_bin)

        # Refine seeds: high gradient of the refinement fields on a fine scatter
        # Sample on a regular grid clipped to the data hull.
        from scipy.spatial import Delaunay as _Delaunay
        hull_delaunay = _Delaunay(hull_verts)

        n_grid = 400
        zz = np.linspace(0.0, 1.0, n_grid)
        hh = np.linspace(0.0, 1.0, n_grid)
        ZZ, HH = np.meshgrid(zz, hh)
        grid_pts = np.column_stack([ZZ.ravel(), HH.ravel()])
        inside   = hull_delaunay.find_simplex(grid_pts) >= 0
        grid_in  = grid_pts[inside]                    # (M, 2) normalized

        grid_spacing = 1.0 / n_grid  # Normalized grid spacing
        scaled_spatial_threshold = grid_spacing * self._nn_edgesize_scaling

        # Build full 3-column CV array for evaluation
        n_cv = len(self._controlling_variables)

        # For each grid point, apply 2-stage filtering and interpolate
        print(f"  Evaluating refinement field on {len(grid_in)} grid points with Z-bin filtering...")
        print(f"    Spatial threshold (scaled): {scaled_spatial_threshold:.6f} = {grid_spacing:.6f} × {self._nn_edgesize_scaling}")
        Q_interp_list = []

        # Build k-d tree for distance queries on 2D data
        kdtree_2d = cKDTree(CV_2d_bin_norm)

        insufficient_count = 0
        for i, grid_pt in enumerate(grid_in):
            if (i + 1) % 50000 == 0:
                print(f"    Processed {i+1}/{len(grid_in)} grid points")

            # Query nearest neighbors in (PV, H) space (candidate pool)
            distances_2d, indices_20 = kdtree_2d.query(grid_pt, k=min(self._n_near, len(D_z_bin)))

            # Get Z distances for these neighbors
            z_neighbors = D_z_bin[indices_20, level_cv_idx]
            z_distances = np.abs(z_neighbors - z_target)

            # Apply 2-stage filtering
            filtered_idx, had_insufficient = self.__FilterNearestNeighbors_2Stage(
                indices_20, distances_2d, z_distances, scaled_spatial_threshold
            )
            if had_insufficient:
                insufficient_count += 1

            # IDW interpolate using filtered neighbors
            filtered_data = D_z_bin[filtered_idx]
            filtered_distances = distances_2d[filtered_idx]

            # Inverse distance weighting
            with np.errstate(divide='ignore', invalid='ignore'):
                weights = 1.0 / np.power(filtered_distances + 1e-16, self._p_fac)
                weights /= weights.sum()
            q_val = weights @ filtered_data

            Q_interp_list.append(q_val)

        Q_interp = np.array(Q_interp_list)
        print(f"  Grid points without sufficient neighbors within scaled edge size: {insufficient_count}/{len(grid_in)} ({100*insufficient_count/len(grid_in):.1f}%)")

        missing = [f for f in self._refinement_fields if f not in self._Flamelet_Variables]
        if missing:
            raise Exception("Refinement field(s) not found in flamelet data: %s" % missing)

        combined_indicator = np.zeros(len(grid_in))
        for field in self._refinement_fields:
            q_vals = Q_interp[:, self._Flamelet_Variables.index(field)]
            q_norm = (q_vals - q_vals.min()) / (q_vals.max() - q_vals.min() + 1e-32)
            # Approximate gradient magnitude via neighbour differences in the grid.
            q_grid = np.full(n_grid * n_grid, np.nan)
            q_grid[inside] = q_norm
            q_grid = q_grid.reshape(n_grid, n_grid)
            dqdy, dqdx = np.gradient(np.nan_to_num(q_grid))
            dq_mag = np.sqrt(dqdx**2 + dqdy**2)
            dq_norm = (dq_mag / (dq_mag.max() + 1e-32)).ravel()[inside]
            combined_indicator = np.maximum(combined_indicator, dq_norm)

        idx_ref = combined_indicator > self._curvature_threshold
        XY_refinement = grid_in[idx_ref]

        # level_norm in the 3D scaler space: PV=0 normalised.
        level_norm = self._scaler.transform(np.zeros([1, n_cv]))[0, self._level_cv_idx]

        # Visualization: Stage 1b - Data after 2-stage filtering (before mesh generation)
        try:
            vis_path_1b = os.path.join(self._savedir, "zbin_stage1_after_filtering.png") if hasattr(self, '_savedir') else None
            print(f"  [DEBUG] Stage 1b visualization: savedir={self._savedir if hasattr(self, '_savedir') else 'N/A'}, path={vis_path_1b}")
            self.__VisualizeFlameletData3D(D_z_bin, f"Stage 1b: Z-bin After 2-Stage Filtering (spatial+Z)", save_path=vis_path_1b)
        except Exception as e:
            print(f"  [ERROR] Stage 1b visualization failed: {e}")

        return XY_refinement, XY_hull, data_hull.area, level_norm


    def __ComputeCurvature(self, val_level:float):
        """
        Compute the curvature of the reaction rate surface at a constant level-variable slice.
        Identify the locations of high curvature where table refinement is required.

        :param val_level: value of the level controlling variable for the current table level.
        :type val_level: float
        :return XY_refinement: normalized in-plane coordinates where refinement should be applied.
        :rtype XY_refinement: array
        :return XY_hull: normalized in-plane coordinates of the convex hull of the current table level.
        :rtype XY_hull: array
        """

        # 1: Get physical ranges of the two plane CVs at this level value.
        plane0_unb, plane0_b, plane1_min, plane1_max, plane1_at_unb = self.__GetPlaneCVBounds(val_level)

        p_idxs = self._plane_cv_idxs
        n_cv   = len(self._controlling_variables)
        level_cv_idx = self._level_cv_idx

        # For multi-Z datasets: Extract and visualize Z-bin at this level
        z_target = val_level
        mask_z_bin = np.abs(self._D_full[:, level_cv_idx] - z_target) <= self._delta_z_bin
        D_z_bin = self._D_full[mask_z_bin]

        if len(D_z_bin) > 0:
            print(f"  Z-bin filtering: {len(D_z_bin)} of {len(self._D_full)} points in Z ± {self._delta_z_bin}")
            # Visualization: Stage 1 - Initial Z-bin data (multi-Z case)
            try:
                vis_path_1 = os.path.join(self._savedir, "zbin_stage1_initial_multiz.png") if hasattr(self, '_savedir') else None
                print(f"  [DEBUG] Stage 1 (multi-Z) visualization: savedir={self._savedir if hasattr(self, '_savedir') else 'N/A'}, path={vis_path_1}")
                self.__VisualizeFlameletData3D(D_z_bin, f"Stage 1: Initial Z-bin Data - Multi-Z (Z={z_target:.6f} ± {self._delta_z_bin:.6f})", save_path=vis_path_1)
            except Exception as e:
                print(f"  [ERROR] Stage 1 (multi-Z) visualization failed: {e}")
        else:
            print(f"  Warning: Empty Z-bin at Z={z_target} ± {self._delta_z_bin}")

        # Define 2D grid between minimum and maximum plane controlling variables.
        n_grid = self._curvature_grid_resolution
        plane0_range = np.linspace(plane0_unb, plane0_b, n_grid)
        plane1_range = np.linspace(plane1_min, plane1_max, n_grid)
        xgrid, ygrid = np.meshgrid(plane0_range, plane1_range)

        n_pts = xgrid.size
        CV_grid_init = np.zeros([n_pts, n_cv])
        CV_grid_init[:, p_idxs[0]]          = xgrid.flatten()
        CV_grid_init[:, p_idxs[1]]          = ygrid.flatten()
        CV_grid_init[:, self._level_cv_idx] = val_level

        # 2: Locate nodes that are above the burner-stabilized boundary line.
        plane0_grid = CV_grid_init[:, p_idxs[0]]
        plane1_grid = CV_grid_init[:, p_idxs[1]]
        h_limit = ((plane1_at_unb - plane1_min) * plane0_grid +
                   (plane1_min * plane0_unb - plane1_at_unb * plane0_b)) / (plane0_unb - plane0_b)
        idx_keep = plane1_grid >= h_limit

        CV_grid = CV_grid_init[idx_keep, :]

        CV_grid_norm_init = self._scaler.transform(CV_grid_init)
        CV_grid_norm      = self._scaler.transform(CV_grid)

        print(f"  Grid stats: {len(CV_grid_init)} initial points, {len(CV_grid)} kept after filtering ({100*len(CV_grid)/len(CV_grid_init):.1f}%)")

        # 3: Generate convex hull on the plane coordinates.
        hull   = ConvexHull(CV_grid_norm[:, p_idxs])
        x_hull = CV_grid_norm[hull.vertices, p_idxs[0]]
        y_hull = CV_grid_norm[hull.vertices, p_idxs[1]]
        print(f"  Convex hull has {len(hull.vertices)} vertices from {len(CV_grid_norm)} data points")

        # 4: Locate refinement locations based on the combined indicator across all refinement fields.
        missing = [f for f in self._refinement_fields if f not in self._Flamelet_Variables]
        if missing:
            raise Exception("Refinement field(s) not found in flamelet data: %s. "
                            "Available fields: %s" % (missing, self._Flamelet_Variables))
        print(f"  Evaluating flamelet interpolator for {len(CV_grid_init)} points with Z-bin filtering...")

        # For multi-Z datasets, compute grid Delaunay for connectivity-based thresholds
        if self._scaler_2d is None:
            # Build Delaunay on grid for connected edge distances
            grid_delaunay = Delaunay(CV_grid_norm_init[:, p_idxs])
            Q_interp = self.__EvaluateFlameletInterpolator(
                CV_unscaled=CV_grid_init,
                apply_z_filtering=True,
                mesh_nodes_norm=CV_grid_norm_init,
                mesh_simplices=grid_delaunay.simplices
            )
        else:
            # Single-Z: use standard evaluation
            Q_interp = self.__EvaluateFlameletInterpolator(CV_unscaled=CV_grid_init)

        print(f"  Done evaluating interpolator")
        combined_indicator = np.zeros(n_pts)
        for field in self._refinement_fields:
            q_grid = np.reshape(Q_interp[:, self._Flamelet_Variables.index(field)], np.shape(xgrid))
            if self._refinement_method == "curvature":
                indicator = self.__ComputeSourceTermCurvature(q_grid)
            else:
                indicator = self.__ComputeSourceTermGradient(q_grid)
            combined_indicator = np.maximum(combined_indicator, indicator)
        idx_ref = np.where(combined_indicator > self._curvature_threshold)

        x_refinement = CV_grid_norm_init[idx_ref, p_idxs[0]]
        y_refinement = CV_grid_norm_init[idx_ref, p_idxs[1]]

        # 5: Generate refinement locations at the unburnt and burnt plane0 boundaries.
        plane1_unb_range = np.linspace(plane1_at_unb, plane1_max, self._Config.GetNpTemp())
        CV_unb = np.zeros([len(plane1_unb_range), n_cv])
        CV_unb[:, p_idxs[0]]          = plane0_unb
        CV_unb[:, p_idxs[1]]          = plane1_unb_range
        CV_unb[:, self._level_cv_idx] = val_level
        CV_unb_norm = self._scaler.transform(CV_unb)

        plane1_b_range = np.linspace(plane1_min, plane1_max, self._Config.GetNpTemp())
        CV_b = np.zeros([len(plane1_b_range), n_cv])
        CV_b[:, p_idxs[0]]          = plane0_b
        CV_b[:, p_idxs[1]]          = plane1_b_range
        CV_b[:, self._level_cv_idx] = val_level
        CV_b_norm = self._scaler.transform(CV_b)

        x_refinement = np.append(x_refinement, CV_unb_norm[:, p_idxs[0]])
        x_refinement = np.append(x_refinement, CV_b_norm[:,  p_idxs[0]])
        y_refinement = np.append(y_refinement, CV_unb_norm[:, p_idxs[1]])
        y_refinement = np.append(y_refinement, CV_b_norm[:,  p_idxs[1]])

        XY_refinement = np.vstack((x_refinement, y_refinement)).T
        XY_hull       = np.vstack((x_hull, y_hull)).T

        val_level_norm = CV_grid_norm[0, self._level_cv_idx]

        # Visualization: Stage 1b - Data after grid evaluation (multi-Z case)
        if len(D_z_bin) > 0:
            try:
                vis_path_1b = os.path.join(self._savedir, "zbin_stage1b_after_filtering_multiz.png") if hasattr(self, '_savedir') else None
                print(f"  [DEBUG] Stage 1b (multi-Z) visualization: savedir={self._savedir if hasattr(self, '_savedir') else 'N/A'}, path={vis_path_1b}")
                self.__VisualizeFlameletData3D(D_z_bin, f"Stage 1b: Z-bin After Grid Evaluation - Multi-Z", save_path=vis_path_1b)
            except Exception as e:
                print(f"  [ERROR] Stage 1b (multi-Z) visualization failed: {e}")

        return XY_refinement, XY_hull, hull.area, val_level_norm, CV_grid, Q_interp

    def __ComputeSourceTermGradient(self, Q_grid:np.ndarray[float]) -> np.ndarray:
        """Return the flattened normalized gradient magnitude for a 2D field array."""
        Q_norm = (Q_grid - np.min(Q_grid)) / (np.max(Q_grid) - np.min(Q_grid) + 1e-32)
        dQdy, dQdx = np.gradient(Q_norm)
        dQ_mag = np.sqrt(np.power(dQdy, 2) + np.power(dQdx, 2))
        return (dQ_mag / (np.max(dQ_mag) + 1e-32)).flatten()

    def __VisualizeFlameletData3D(self, data:np.ndarray, title:str, save_path:str=None, skip_points:int=10):
        """Nijso: For debugging: Create a 3D scatter plot of flamelet data (PV, H, T) colored by Z.

        :param data: Flamelet data array (N, n_vars).
        :param title: Title for the plot.
        :param save_path: If provided, save to this file instead of showing.
        :param skip_points: Plot every skip_points'th point to reduce clutter.
        """
        try:
            from mpl_toolkits.mplot3d import Axes3D
        except ImportError:
            print(f"    Skipping 3D visualization (matplotlib 3D not available): {title}")
            return

        p0_idx = self._plane_cv_idxs[0]  # PV or Z
        p1_idx = self._plane_cv_idxs[1]  # H
        level_idx = self._level_cv_idx   # Z or PV

        # Assume Temperature is near the end of the data columns (check based on _Flamelet_Variables)
        try:
            temp_idx = self._Flamelet_Variables.index('Temperature')
        except ValueError:
            print(f"    Skipping 3D visualization (Temperature column not found): {title}")
            return

        subset_plot = data[::skip_points]

        if len(subset_plot) == 0:
            print(f"    Skipping 3D visualization (empty data): {title}")
            return

        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')

        scatter = ax.scatter(
            subset_plot[:, p0_idx],
            subset_plot[:, p1_idx],
            subset_plot[:, temp_idx],
            c=subset_plot[:, level_idx],
            cmap='coolwarm',
            marker='o',
            s=20,
            alpha=0.6,
            edgecolors='k',
            linewidth=0.3
        )

        ax.set_xlabel(self._plane_cv_names[0], fontsize=12)
        ax.set_ylabel(self._plane_cv_names[1], fontsize=12)
        ax.set_zlabel('Temperature (K)', fontsize=12)
        ax.set_title(f'{title}\n({len(subset_plot)}/{len(data)} points shown, SKIP={skip_points})', fontsize=14)

        cbar = fig.colorbar(scatter, ax=ax, pad=0.1, shrink=0.8)
        cbar.set_label(self._level_cv_name, fontsize=11)

        ax.view_init(elev=20, azim=45)
        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"    Saved visualization: {save_path}")
        else:
            plt.show()
        plt.close(fig)

    def __ComputeSourceTermCurvature(self, Q_grid:np.ndarray[float]) -> np.ndarray:
        """Return the flattened normalized curvature magnitude for a 2D field array."""
        Q_norm = (Q_grid - np.min(Q_grid)) / (np.max(Q_grid) - np.min(Q_grid) + 1e-32)
        dQdy, dQdx = np.gradient(Q_norm)
        dQ_mag = np.sqrt(np.power(dQdy, 2) + np.power(dQdx, 2))
        dQ_norm = dQ_mag / (np.max(dQ_mag) + 1e-32)
        d2Qdy2, d2Qdx2 = np.gradient(dQ_norm)
        d2Q_mag = np.sqrt(np.power(d2Qdy2, 2) + np.power(d2Qdx2, 2))
        return (d2Q_mag / (np.max(d2Q_mag) + 1e-32)).flatten()

    def __Compute2DMesh(self, XY_hull:np.ndarray, XY_refinement:np.ndarray, val_level_norm:float, level_area:float):
        """
        Generate a 2D mesh for the current table level.

        :param XY_hull: Array containing normalized pv and enth coordinates of the outline of the table level.
        :type XY_hull: NDArray
        :param XY_refinement: Array containing normalized pv and enth coordinates where refinement should be applied.
        :type XY_refinement: NDArray
        :return: mesh nodes of the 2D table mesh.
        :rtype: NDArray
        """
        # Initialize gmsh or clear existing model
        if not gmsh.isInitialized():
            gmsh.initialize()

        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 1)
        gmsh.model.add("table_level")
        factory = gmsh.model.geo

        base_cell_size = self._base_cell_size #* level_area
        refined_cell_size = self._refined_cell_size #* level_area
        hull_cell_size = self._convex_hull_cell_size
        refinement_radius = self._refinement_radius #* np.sqrt(level_area)
        print("Generating 2D mesh with base cell size %.4f, hull cell size %.4f and refined cell size %.4f" % (base_cell_size, hull_cell_size, refined_cell_size))

        # Create hull boundary points
        hull_pts = []
        for i in range(len(XY_hull)):
            hull_pts.append(factory.addPoint(XY_hull[i, 0], XY_hull[i, 1], 0, hull_cell_size))

        # Subsample refinement seed points to avoid excessive PointsList size.
        N_max_seeds = self._max_refinement_seeds
        if len(XY_refinement) > N_max_seeds:
            idx_sub = np.round(np.linspace(0, len(XY_refinement) - 1, N_max_seeds)).astype(int)
            XY_refinement_sub = XY_refinement[idx_sub]
        else:
            XY_refinement_sub = XY_refinement

        print(f"  Creating {len(XY_refinement_sub)} refinement seed points")
        embed_pts = []
        for i in range(len(XY_refinement_sub)):
            pt_idx = factory.addPoint(XY_refinement_sub[i, 0], XY_refinement_sub[i, 1], 0, refined_cell_size)
            embed_pts.append(pt_idx)
        print(f"  Created {len(embed_pts)} embed points")

        # Create individual line segments for the hull boundary
        hull_lines = []
        for i in range(len(hull_pts)):
            line = factory.addLine(hull_pts[i], hull_pts[(i+1) % len(hull_pts)])
            hull_lines.append(line)

        CL = factory.addCurveLoop(hull_lines)

        surf = factory.addPlaneSurface([CL])
        gmsh.model.addPhysicalGroup(1, hull_lines, name="hull_boundary")
        gmsh.model.addPhysicalGroup(2, [surf], name="table_level")
        gmsh.model.geo.synchronize()

        # Let the background mesh field fully control element sizes.
        # Disabling these prevents the coarse hull boundary size from spreading
        # inward and overriding the Threshold field in the refined region.
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)

        print(f"  Setting up mesh fields with {len(embed_pts)} refinement points...")

        # Only create refinement distance field if we have refinement points
        if len(embed_pts) > 0:
            gmsh.model.mesh.field.add("Distance", 1)
            gmsh.model.mesh.field.setNumbers(1, "PointsList", embed_pts)
            gmsh.model.mesh.field.setNumber(1, "Sampling", 100)
            gmsh.model.mesh.field.add("Threshold", 2)
            gmsh.model.mesh.field.setNumber(2, "InField", 1)
            gmsh.model.mesh.field.setNumber(2, "SizeMin", refined_cell_size)
            gmsh.model.mesh.field.setNumber(2, "SizeMax", base_cell_size)
            gmsh.model.mesh.field.setNumber(2, "DistMin", refinement_radius)
            gmsh.model.mesh.field.setNumber(2, "DistMax", 1.5*refinement_radius)

        # Refine along the convex hull boundary curves.
        gmsh.model.mesh.field.add("Distance", 3)
        gmsh.model.mesh.field.setNumbers(3, "CurvesList", hull_lines)
        gmsh.model.mesh.field.setNumber(3, "Sampling", 200)
        gmsh.model.mesh.field.add("Threshold", 4)
        gmsh.model.mesh.field.setNumber(4, "InField", 3)
        gmsh.model.mesh.field.setNumber(4, "SizeMin", hull_cell_size)
        gmsh.model.mesh.field.setNumber(4, "SizeMax", base_cell_size)
        gmsh.model.mesh.field.setNumber(4, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(4, "DistMax", refinement_radius)

        gmsh.model.mesh.field.add("Min", 7)
        if len(embed_pts) > 0:
            gmsh.model.mesh.field.setNumbers(7, "FieldsList", [2, 4])
        else:
            gmsh.model.mesh.field.setNumbers(7, "FieldsList", [4])
        gmsh.model.mesh.field.setAsBackgroundMesh(7)

        print("  Generating mesh...")
        gmsh.option.setNumber("Mesh.Algorithm", 5)
        gmsh.model.mesh.generate(2)
        print("  Mesh generated successfully!")
        print("  Extracting nodes...")
        nodes = gmsh.model.mesh.getNodes(dim=2, tag=-1, includeBoundary=True, returnParametricCoord=False)[1]
        print(f"  Extracted {len(nodes)//3} nodes")
        MeshPoints = np.array([nodes[::3], nodes[1::3]]).T
        print(f"  MeshPoints shape: {MeshPoints.shape}")

        # Don't finalize or clear - gmsh has issues with cleanup in some environments
        # This causes a small memory leak but avoids segfaults
        print("  Skipping gmsh cleanup to avoid segfault...")

        # Build the full CV array in controlling-variable column order.
        print("  Building CV array...")
        plane0_norm = MeshPoints[:, 0]
        plane1_norm = MeshPoints[:, 1]
        level_norm  = val_level_norm * np.ones(len(plane0_norm))

        n_cv = len(self._controlling_variables)
        CV_level_norm = np.zeros([len(plane0_norm), n_cv])
        CV_level_norm[:, self._plane_cv_idxs[0]] = plane0_norm
        CV_level_norm[:, self._plane_cv_idxs[1]] = plane1_norm
        CV_level_norm[:, self._level_cv_idx]      = level_norm
        CV_level_dim = self._scaler.inverse_transform(CV_level_norm)

        # Compute Delaunay triangulation for connectivity (used for Stage 3 filtering)
        print("  Computing Delaunay triangulation...")
        Tria = Delaunay(MeshPoints[:, self._plane_cv_idxs])
        mesh_simplices = Tria.simplices

        print(f"  Evaluating flamelet interpolator for {len(CV_level_dim)} mesh points...")
        table_level_data = self.__EvaluateFlameletInterpolator(
            CV_level_dim,
            apply_z_filtering=True,
            mesh_nodes_norm=MeshPoints,
            mesh_simplices=mesh_simplices
        )
        print("  Done evaluating")


        return CV_level_norm, table_level_data

    def __GetStochMixtureFraction(self):
        fuel_definition = self._Config.GetFuelDefinition()
        fuel_weights = self._Config.GetFuelWeights()
        fuel_string = ",".join(fuel_definition[i] + ":" + str(fuel_weights[i]) for i in range(len(fuel_definition)))

        ox_definition = self._Config.GetOxidizerDefinition()
        ox_weights = self._Config.GetOxidizerWeights()
        ox_string = ",".join(ox_definition[i] + ":" + str(ox_weights[i]) for i in range(len(ox_definition)))

        self._Config.gas.set_equivalence_ratio(1.0, fuel_string, ox_string)
        mixfrac_stoch = self._Config.gas.mixture_fraction(fuel_string, ox_string)
        return mixfrac_stoch

    def ClampSourceTerms(self, species_list:list[str], pv_frac:float=0.99, abs_tol:float=1e-3):
        """Clamp source terms of selected species to zero near the burnt (high-PV) boundary.

        For each table level the per-level PV maximum is computed.  At every node
        where PV >= ``pv_frac * PV_max`` the net, positive, and negative source
        terms (``Y_dot_net-<sp>``, ``Y_dot_pos-<sp>``, ``Y_dot_neg-<sp>``) for each
        species in ``species_list`` are forced to zero when their absolute value is
        below ``abs_tol``.

        This method must be called after :meth:`GenerateTableNodes` and before
        :meth:`WriteTableFile`.

        :param species_list: species names to clamp (e.g. ``['CO', 'H2']``).
        :type species_list: list[str]
        :param pv_frac: fraction of per-level PV_max above which clamping is applied; default 0.99.
        :type pv_frac: float
        :param abs_tol: source terms with absolute value below this are set to zero; default 1e-3.
        :type abs_tol: float
        """
        if not species_list:
            return
        pv_col    = self._plane_cv_idxs[0]   # column of ProgressVariable in _table_nodes

        # Build per-species index lists for each source-term role.
        idx_net, idx_pos, idx_neg = [], [], []
        for sp in species_list:
            for fmt, bucket in [("Y_dot_net-%s", idx_net),
                                 ("Y_dot_pos-%s", idx_pos),
                                 ("Y_dot_neg-%s", idx_neg)]:
                vname = fmt % sp
                if vname in self.table_vars:
                    bucket.append(self.table_vars.index(vname))
                else:
                    print("ClampSourceTerms: variable '%s' not found in table_vars, skipping." % vname)

        clamp_indices = idx_net + idx_pos + idx_neg
        if not clamp_indices:
            return

        n_clamped_total = 0
        for iLevel in range(self._N_table_levels):
            pv_nodes = self._table_nodes[iLevel][:, pv_col]
            pv_max   = np.max(pv_nodes)
            mask     = pv_nodes >= pv_frac * pv_max
            n_clamped_total += int(np.sum(mask))

            # Near the burnt boundary: zero small-magnitude net/pos/neg source terms.
            for iVar in clamp_indices:
                src = self.table_data[iLevel][iVar]
                src[mask & (np.abs(src) < abs_tol)] = 0.0
                self.table_data[iLevel][iVar] = src

            # Sign constraints (all nodes):
            #   source_prod (Y_dot_pos) >= 0  — production is never negative.
            #   source_cons (Y_dot_neg) <= 0  — consumption coefficient is never positive
            #   (SU2 uses:  source = source_prod + source_cons * y_aux).
            for iVar in idx_pos:
                src = self.table_data[iLevel][iVar]
                np.clip(src, 0.0, None, out=src)
                self.table_data[iLevel][iVar] = src
            for iVar in idx_neg:
                src = self.table_data[iLevel][iVar]
                np.clip(src, None, 0.0, out=src)
                self.table_data[iLevel][iVar] = src

        print("ClampSourceTerms: clamped %i nodes across %i levels (pv_frac=%.3f, abs_tol=%.2e)." %
              (n_clamped_total, self._N_table_levels, pv_frac, abs_tol))
        return

    def SaveTableGenerator(self, file_name:str):
        """Save the current TableGenerator object settings such that subsequent tables can be
        generated faster.

        :param file_name: file path and name to which to save the current TableGenerator.
        :type file_name: str
        """
        file = open(self._savedir + "/"+file_name +".tgen", "wb")
        pickle.dump(self, file)
        file.close()

    def Inverse_LookUp_T(self, val_pv, val_mixfrac, val_T, val_h_start=2000):
        CV_array = np.array([[val_pv, val_h_start, val_mixfrac]])
        delta = 1e32
        while np.abs(delta) > 1e-2:
            Q_interp = self.__EvaluateFlameletInterpolator(CV_array)
            val_T_interp = Q_interp[0, self._Flamelet_Variables.index("Temperature")]
            val_cp_interp  = Q_interp[0, self._Flamelet_Variables.index("Cp")]
            delta = val_T - val_T_interp
            delta_h = val_cp_interp * delta
            CV_array[0,1] += delta_h
        return CV_array[0,1]

if __name__ == "__main__":
    config_input_file = sys.argv[-2]
    N_cores = int(sys.argv[-1])
    Config = Config_FGM(config_input_file)
    T = SU2TableGenerator(Config)
    if N_cores > 1:
        T.SetNCores(N_cores)
    T.SetMixtureFractionLimits(mix_frac_min=0.009, mix_frac_max=0.022)
    T.InsertMixtureFractionLevel(0.01446751783896619)
    T.InsertMixtureFractionLevel(0.01447)
    T.InsertMixtureFractionLevel(0.01445)
    T.VisualizeTableLevel(0.01446751783896619)
    T.SetNTableLevels(200)
    T.GenerateTableNodes()
    T.WriteTableFile()
    # #T.InterpolateTableData()
    # T.WriteTableFile()
    # T.SaveTableGenerator("LUT_"+Config.GetConfigName())
