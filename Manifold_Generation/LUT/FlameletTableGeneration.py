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
# Version: 3.2.0                                                                              |
#                                                                                             |
#=============================================================================================#

import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys,os
from Common.DataDrivenConfig import Config_FGM, Config
from Common.CommonMethods import readStateDataFromFile
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

        self._control_var_scaler = MinMaxScaler()
        CV_full, D_full = readStateDataFromFile(full_data_file, self._controlling_variables, self._manifold_variables)
        data_scaler = MinMaxScaler()
        data_scaler.fit_transform(D_full)

        CV_full_scaled = self._control_var_scaler.fit_transform(CV_full)

        # Exctract train and test data
        train_data_file = self._Config.GetOutputDir()+"/"+self._Config.GetConcatenationFileHeader()+"_train.csv"
        test_data_file = self._Config.GetOutputDir()+"/"+self._Config.GetConcatenationFileHeader()+"_test.csv"

        CV_train, D_train = readStateDataFromFile(train_data_file, self._controlling_variables, self._manifold_variables)
        CV_test, D_test = readStateDataFromFile(test_data_file, self._controlling_variables, self._manifold_variables)

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
        print(min(RMS_ppv))
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

    _user_vars:list[str] = []
    _table_vars:list[str] = None
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

    _table_nodes = []       # Progress variable, total enthalpy, and mixture fraction node values for each table level.
    _table_data:list[np.ndarray[float]] = []
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
        self.__target_node_count = False
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
        full_data_file = "/".join((self._Config.GetOutputDir(), "%s_full.csv" % self._Config.GetConcatenationFileHeader()))

        with open(full_data_file,'r') as fid:
            self._Flamelet_Variables = fid.readline().strip().split(',')
        D_full = np.loadtxt(full_data_file,delimiter=',',skiprows=1)
        self._scaler = MinMaxScaler()
        CV_full = D_full[:,[self._Flamelet_Variables.index(cv) for cv in self._Config.GetControllingVariables()]]
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
        self._lookup_tree = Invdisttree(samples_state_space=CV_train_scaled,state_data=D_train)
        print("Done!")

        if not self._custom_KDtreeparams:
            # Do brute-force search to get the optimum number of nearest neighbors and distance power.
            n_near_range = range(2, 30)
            p_range = range(1, 6)
            RMS_ppv = np.zeros([len(n_near_range), len(p_range)])
            for i in tqdm(range(len(n_near_range)), desc="Searching for best tree parameters"):
                for j in range(len(p_range)):
                    PPV_predicted = self._lookup_tree.query(query_samples=CV_test_scaled, nearest_neighbors=n_near_range[i], inverse_distance_exponent=p_range[j])[:, self._Flamelet_Variables.index(var_to_test_for)]
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
        self._lookup_tree = Invdisttree(samples_state_space=CV_full_scaled,state_data=D_full)
        print("Done!")
        return

    def __Rebuild2DInterpolator(self):
        """Rebuild the KD-tree using only (PV, H) coordinates for a single-Z table.

        When a dataset has only one equivalence ratio/mixture fraction, the Z column
        in the full 3D MinMax-scaled space has a range that reflects differential-diffusion
        variation along the flame (not a true Z sweep).  That tiny real-world Δ is
        stretched to [0, 1] by the scaler, making Z-distances dominate the neighbour
        search and causing unphysical jumps.  Dropping Z from the search space fixes this.
        """
        if self._D_full is None:
            return  # interpolator not yet built
        print("Rebuilding 2D (%s) KD-tree for single-level table..." % ", ".join(self._plane_cv_names))
        CV_full_2d = self._D_full[:, self._plane_cv_idxs]
        self._scaler_2d = MinMaxScaler()
        CV_full_2d_scaled = self._scaler_2d.fit_transform(CV_full_2d)
        self._lookup_tree = Invdisttree(samples_state_space=CV_full_2d_scaled, state_data=self._D_full)
        print("Done!")
        return

    def __EvaluateFlameletInterpolator(self, CV_unscaled:np.ndarray):
        if self._is_2D_table and self._scaler_2d is not None:
            CV_scaled = self._scaler_2d.transform(CV_unscaled[:, self._plane_cv_idxs])
        else:
            CV_scaled = self._scaler.transform(CV_unscaled)
        data_interp = self._lookup_tree.query(query_samples=CV_scaled,nearest_neighbors=self._n_near,inverse_distance_exponent=self._p_fac)
        return data_interp

    def VisualizeTableLevel(self, val_mix_frac:float, var_to_plot:str=None):
        """Compute and visualize the table connectivity for a certain mixture fraction value.

        :param val_mix_frac: mixture fraction value for which to compute the table connectivity.
        :type val_mix_frac: float
        :raises Exception: if the mixture fraction value lies outside the flamelet data range.
        """

        Tria, Nodes, HullIdx, level_data, XY_ref_dim = self.ComputeTableLevelMesh(val_mix_frac)
        print("Total mesh nodes: %i, refinement seed points: %i" % (len(Nodes), len(XY_ref_dim)))

        p0, p1 = self._plane_cv_idxs
        if var_to_plot == None:
            _ = plt.figure(figsize=[10,10])
            ax = plt.axes()
            ax.triplot(Nodes[:, p0], Nodes[:, p1], Tria)
            ax.plot(Nodes[HullIdx, p0], Nodes[HullIdx, p1], 'ko', label=r"Hull nodes")
            ax.set_xlabel(self._plane_cv_names[0], fontsize=20)
            ax.set_ylabel(self._plane_cv_names[1], fontsize=20)
            ax.legend(fontsize=20)
            ax.set_title(self._level_cv_name + " = " + str(val_mix_frac))
            plt.show()
        else:
            var_idx = self._Flamelet_Variables.index(var_to_plot)
            fig = plt.figure(figsize=[10,10])
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_trisurf(Nodes[:, p0], Nodes[:, p1], level_data[:, var_idx],
                            triangles=Tria, cmap='viridis', alpha=0.9, edgecolor='k', linewidth=0.2)
            ax.set_xlabel(self._plane_cv_names[0], fontsize=20)
            ax.set_ylabel(self._plane_cv_names[1], fontsize=20)
            ax.set_title(self._level_cv_name + " = " + str(val_mix_frac))
            plt.show()
        return

    def generateTableNodes(self):
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
                self._table_data[iLevel] = results[iLevel][3]
                for icv, cv in enumerate(self._Config.GetControllingVariables()):
                    self._table_data[iLevel][:, self._Flamelet_Variables.index(cv)] = self._table_nodes[iLevel][:, icv]

                NTria += np.shape(self._table_connectivity[iLevel])[0]
                NHull += np.shape(self._table_hullnodes[iLevel])[0]
                NNodes += np.shape(self._table_nodes[iLevel])[0]
        else:
            for iLevel in range(self._N_table_levels):
                result = self.ComputeTableNodes(iLevel)
                self._table_nodes[iLevel] = result[0]
                self._table_connectivity[iLevel] = result[1]
                self._table_hullnodes[iLevel] = result[2]
                self._table_data[iLevel] = result[3]
                for icv, cv in enumerate(self._Config.GetControllingVariables()):
                    self._table_data[iLevel][:, self._Flamelet_Variables.index(cv)] = self._table_nodes[iLevel][:, icv]
                
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


    def writeSU2Table(self, output_filepath:str=None):
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

        if len(self._user_vars) > 0:
            table_vars = self._user_vars.copy()
        else:
            table_vars = self._Flamelet_Variables.copy()

        print("Writing LUT file with name " + file_out)
        fid = open(file_out, "w+")

        if self._is_2D_table:
            self.__writeSU2Table2D(fid)
        else:
            self.__writeSU2Table3D(fid)

        fid.close()
        return

    def __writeSU2Table2D(self, fid):
        """Write a 2D Dragon library file (version 1.0.1) for a single mixture-fraction level.
        The controlling variables are ProgressVariable and EnthalpyTot only."""

        Nodes = self._table_nodes[0]          # shape (Np, 3): [PV, H, Z]
        Connectivity = self._table_connectivity[0]
        HullNodes = self._table_hullnodes[0]
        Np = np.shape(Nodes)[0]

        fid.write("Dragon library\n\n")
        fid.write("<Header>\n\n")
        fid.write("[Version]\n1.0.1\n\n")
        fid.write("[Number of points]\n%i\n\n" % Np)
        fid.write("[Number of triangles]\n%i\n\n" % np.shape(Connectivity)[0])
        fid.write("[Number of hull points]\n%i\n\n" % np.shape(HullNodes)[0])

        fid.write("[Progress variable definition]\n")
        fid.write("+".join(("%+.4e * %s" % (w, s)) for w, s in zip(
            self._Config.GetProgressVariableWeights(),
            self._Config.GetProgressVariableSpecies())) + "\n\n")

        pv_vals = Nodes[:, self._plane_cv_idxs[0]]
        h_vals  = Nodes[:, self._plane_cv_idxs[1]]
        fid.write("[ProgressVariable min]\n%e\n\n" % np.min(pv_vals))
        fid.write("[ProgressVariable max]\n%e\n\n" % np.max(pv_vals))
        fid.write("[EnthalpyTot min]\n%e\n\n" % np.min(h_vals))
        fid.write("[EnthalpyTot max]\n%e\n\n" % np.max(h_vals))

        all_vars_2d = ["ProgressVariable", "EnthalpyTot"] + self.table_vars
        fid.write("[Number of variables]\n%i\n\n" % len(all_vars_2d))
        fid.write("[Variable names]\n")
        for var in all_vars_2d:
            fid.write(var + "\n")
        fid.write("\n")

        fid.write("</Header>\n\n")

        print("Writing table data...")
        fid.write("<Data>\n")
        for iNode in tqdm(range(Np)):
            fid.write("%+.14e %+.14e" % (pv_vals[iNode], h_vals[iNode]))
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

    def __writeSU2Table3D(self, fid):
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

        fid.write("[Number of variables]\n%i\n\n" % len(table_vars))
        fid.write("[Variable names]\n")
        for iVar, Var in enumerate(table_vars):
            fid.write(str(iVar + 1)+":"+Var+"\n")
        fid.write("\n")

        fid.write("</Header>\n\n")

        fid.write("<Data>\n")
        for iLevel in tqdm(range(len(self._table_nodes)), desc="Writing table data"):
            fid.write("<Level>\n")
            Np = np.shape(self._table_nodes[iLevel])[0]
            for iNode in range(Np):

                for var in table_vars:
                    fid.write("\t%+.14e" % self._table_data[iLevel][iNode, self._Flamelet_Variables.index(var)])
                
                fid.write("\n")
            fid.write("</Level>\n")
        fid.write("</Data>\n\n")
        print("Done!")

        fid.write("<Connectivity>\n")
        for iLevel in tqdm(range(len(self._table_connectivity)),desc="Writing table connectivity"):
            fid.write("<Level>\n")
            for iCell in range(len(self._table_connectivity[iLevel])):
                fid.write("\t".join("%i" % c for c in self._table_connectivity[iLevel][iCell, :]+1) + "\n")
            fid.write("</Level>\n")
        fid.write("</Connectivity>\n\n")
        print("Done!")

        fid.write("<Hull>\n")
        for iLevel in tqdm(range(len(self._table_hullnodes)),desc="Writing hull nodes"):
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
        Coord_refinement, Coord_hull, hull_area, level_norm, CV_mesh, table_level_data = self.__ComputeCurvature(val_level)
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
        else:
            raise NotImplementedError(
                "Plane CV bounds not implemented for level='%s', plane=%s. "
                "Currently only supported: level='%s', plane=['%s', '%s']." % (
                    self._level_cv_name, self._plane_cv_names,
                    DefaultSettings_FGM.name_mixfrac,
                    DefaultSettings_FGM.name_pv,
                    DefaultSettings_FGM.name_enth))

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

        # Define 2D grid between minimum and maximum plane controlling variables.
        plane0_range = np.linspace(plane0_unb, plane0_b, 800)
        plane1_range = np.linspace(plane1_min, plane1_max, 800)
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

        # 3: Generate convex hull on the plane coordinates.
        hull   = ConvexHull(CV_grid_norm[:, p_idxs])
        x_hull = CV_grid_norm[hull.vertices, p_idxs[0]]
        y_hull = CV_grid_norm[hull.vertices, p_idxs[1]]

        # 4: Locate refinement locations based on the combined indicator across all refinement fields.
        missing = [f for f in self._refinement_fields if f not in self._Flamelet_Variables]
        if missing:
            raise Exception("Refinement field(s) not found in flamelet data: %s. "
                            "Available fields: %s" % (missing, self._Flamelet_Variables))
        Q_interp = self.__EvaluateFlameletInterpolator(CV_unscaled=CV_grid_init)
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

        return XY_refinement, XY_hull, hull.area, val_level_norm, CV_grid, Q_interp

    def __ComputeSourceTermGradient(self, Q_grid:np.ndarray[float]) -> np.ndarray:
        """Return the flattened normalized gradient magnitude for a 2D field array."""
        Q_norm = (Q_grid - np.min(Q_grid)) / (np.max(Q_grid) - np.min(Q_grid) + 1e-32)
        dQdy, dQdx = np.gradient(Q_norm)
        dQ_mag = np.sqrt(np.power(dQdy, 2) + np.power(dQdx, 2))
        return (dQ_mag / (np.max(dQ_mag) + 1e-32)).flatten()

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

        :param XY_hull: Progress variable-total enthalpy locations of the table level outline.
        :type XY_hull: np.ndarray
        :param XY_refinement: Locations where refinement is to be applied.
        :type XY_refinement: np.ndarray
        :param val_mixfrac_norm: Scaled value of the mixture fraction
        :type val_mixfrac_norm: float
        :param level_area: area of the table level.
        :type level_area: float
        :return: table nodes, connectivity, boundary indices, interpolated flamelet data
        :rtype: tuple
        """
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

        hull_pts = []
        for i in range(int(len(XY_hull)/2)):
            hull_pts.append(factory.addPoint(XY_hull[i, 0], XY_hull[i, 1], 0, hull_cell_size))
        hull_pts_2 = [hull_pts[-1]]
        for i in range(int(len(XY_hull)/2), len(XY_hull)):
            hull_pts_2.append(factory.addPoint(XY_hull[i, 0], XY_hull[i, 1], 0, hull_cell_size))
        hull_pts_2.append(hull_pts[0])

        # Subsample refinement seed points to avoid excessive PointsList size.
        N_max_seeds = self._max_refinement_seeds
        if len(XY_refinement) > N_max_seeds:
            idx_sub = np.round(np.linspace(0, len(XY_refinement) - 1, N_max_seeds)).astype(int)
            XY_refinement_sub = XY_refinement[idx_sub]
        else:
            XY_refinement_sub = XY_refinement

        embed_pts = []
        for i in range(len(XY_refinement_sub)):
            pt_idx = factory.addPoint(XY_refinement_sub[i, 0], XY_refinement_sub[i, 1], 0, refined_cell_size)
            embed_pts.append(pt_idx)

        hull_curve_1 = factory.addPolyline(hull_pts)
        hull_curve_2 = factory.addPolyline(hull_pts_2)

        CL = factory.addCurveLoop([hull_curve_1, hull_curve_2])

        surf = factory.addPlaneSurface([CL])
        gmsh.model.addPhysicalGroup(1, [hull_curve_1], name="hull_curve_1")
        gmsh.model.addPhysicalGroup(1, [hull_curve_2], name="hull_curve_2")
        gmsh.model.addPhysicalGroup(2, [surf], name="table_level")
        gmsh.model.geo.synchronize()

        # Let the background mesh field fully control element sizes.
        # Disabling these prevents the coarse hull boundary size from spreading
        # inward and overriding the Threshold field in the refined region.
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)

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
        gmsh.model.mesh.field.setNumbers(3, "CurvesList", [hull_curve_1, hull_curve_2])
        gmsh.model.mesh.field.setNumber(3, "Sampling", 200)
        gmsh.model.mesh.field.add("Threshold", 4)
        gmsh.model.mesh.field.setNumber(4, "InField", 3)
        gmsh.model.mesh.field.setNumber(4, "SizeMin", hull_cell_size)
        gmsh.model.mesh.field.setNumber(4, "SizeMax", base_cell_size)
        gmsh.model.mesh.field.setNumber(4, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(4, "DistMax", refinement_radius)

        gmsh.model.mesh.field.add("Min", 7)
        gmsh.model.mesh.field.setNumbers(7, "FieldsList", [2, 4])
        gmsh.model.mesh.field.setAsBackgroundMesh(7)

        gmsh.option.setNumber("Mesh.Algorithm", 5)
        gmsh.model.mesh.generate(2)
        nodes = gmsh.model.mesh.getNodes(dim=2, tag=-1, includeBoundary=True, returnParametricCoord=False)[1]
        MeshPoints = np.array([nodes[::3], nodes[1::3]]).T

        # we need finalize
        gmsh.finalize()

        # Build the full CV array in controlling-variable column order.
        plane0_norm = MeshPoints[:, 0]
        plane1_norm = MeshPoints[:, 1]
        level_norm  = val_level_norm * np.ones(len(plane0_norm))

        n_cv = len(self._controlling_variables)
        CV_level_norm = np.zeros([len(plane0_norm), n_cv])
        CV_level_norm[:, self._plane_cv_idxs[0]] = plane0_norm
        CV_level_norm[:, self._plane_cv_idxs[1]] = plane1_norm
        CV_level_norm[:, self._level_cv_idx]      = level_norm
        CV_level_dim = self._scaler.inverse_transform(CV_level_norm)

        table_level_data = self.__EvaluateFlameletInterpolator(CV_level_dim)

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

        This method must be called after :meth:`generateTableNodes` and before
        :meth:`writeSU2Table`.

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
    T.generateTableNodes()
    T.writeSU2Table()
    # #T.InterpolateTableData()
    # T.writeSU2Table()
    # T.SaveTableGenerator("LUT_"+Config.GetConfigName())
