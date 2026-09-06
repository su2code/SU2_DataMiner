###############################################################################################
#       #      _____ __  _____      ____        __        __  ____                   #        #
#       #     / ___// / / /__ \    / __ \____ _/ /_____ _/  |/  (_)___  ___  _____   #        #
#       #     \__ \/ / / /__/ /   / / / / __ `/ __/ __ `/ /|_/ / / __ \/ _ \/ ___/   #        #
#       #    ___/ / /_/ // __/   / /_/ / /_/ / /_/ /_/ / /  / / / / / /  __/ /       #        #
#       #   /____/\____//____/  /_____/\__,_/\__/\__,_/_/  /_/_/_/ /_/\___/_/        #        #
#       #                                                                            #        #
###############################################################################################

############################### FILE NAME: Interpolators.py ###################################
#=============================================================================================#
# author: Evert Bunschoten                                                                    |
#    :PhD Candidate ,                                                                         |
#    :Flight Power and Propulsion                                                             |
#    :TU Delft,                                                                               |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#  Interpolator functions used during data processing steps.                                  |
#                                                                                             |
# Version: 3.2.0                                                                              |
#                                                                                             |
#=============================================================================================#

from scipy.spatial import cKDTree as KDTree
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.preprocessing import MinMaxScaler
class Invdisttree:
    """ inverse-distance-weighted interpolation using KDTree:
invdisttree = Invdisttree( X, z )  -- data points, values
interpol = invdisttree( q, nnear=3, eps=0, p=1, weights=None, stat=0 )
    interpolates z from the 3 points nearest each query point q;
    For example, interpol[ a query point q ]
    finds the 3 data points nearest q, at distances d1 d2 d3
    and returns the IDW average of the values z1 z2 z3
        (z1/d1 + z2/d2 + z3/d3)
        / (1/d1 + 1/d2 + 1/d3)
        = .55 z1 + .27 z2 + .18 z3  for distances 1 2 3

    q may be one point, or a batch of points.
    eps: approximate nearest, dist <= (1 + eps) * true nearest
    p: use 1 / distance**p
    weights: optional multipliers for 1 / distance**p, of the same shape as q
    stat: accumulate wsum, wn for average weights

How many nearest neighbors should one take ?
a) start with 8 11 14 .. 28 in 2d 3d 4d .. 10d; see Wendel's formula
b) make 3 runs with nnear= e.g. 6 8 10, and look at the results --
    |interpol 6 - interpol 8| etc., or |f - interpol*| if you have f(q).
    I find that runtimes don't increase much at all with nnear -- ymmv.

p=1, p=2 ?
    p=2 weights nearer points more, farther points less.
    In 2d, the circles around query points have areas ~ distance**2,
    so p=2 is inverse-area weighting. For example,
        (z1/area1 + z2/area2 + z3/area3)
        / (1/area1 + 1/area2 + 1/area3)
        = .74 z1 + .18 z2 + .08 z3  for distances 1 2 3
    Similarly, in 3d, p=3 is inverse-volume weighting.

Scaling:
    if different X coordinates measure different things, Euclidean distance
    can be way off.  For example, if X0 is in the range 0 to 1
    but X1 0 to 1000, the X1 distances will swamp X0;
    rescale the data, i.e. make X0.std() ~= X1.std() .

A nice property of IDW is that it's scale-free around query points:
if I have values z1 z2 z3 from 3 points at distances d1 d2 d3,
the IDW average
    (z1/d1 + z2/d2 + z3/d3)
    / (1/d1 + 1/d2 + 1/d3)
is the same for distances 1 2 3, or 10 20 30 -- only the ratios matter.
In contrast, the commonly-used Gaussian kernel exp( - (distance/h)**2 )
is exceedingly sensitive to distance and to h.

    """
# anykernel( dj / av dj ) is also scale-free
# error analysis, |f(x) - idw(x)| ? todo: regular grid, nnear ndim+1, 2*ndim
    __state_space_data:np.ndarray[float]=None
    __nDim_input:int=None
    __number_of_states:int=None
    __KD_tree:KDTree = None
    __single_query_sample:bool=False
    __number_of_nearest_neighbors:int=None
    __inverse_distance_exponent:float=None
    
    def __init__( self, samples_state_space:np.ndarray[float], state_data:np.ndarray[float], leafsize:int=10):
        assert len(samples_state_space) == len(state_data), "len(X) %d != len(z) %d" % (len(samples_state_space), len(state_data))
        self.__nDim_input = np.shape(samples_state_space)[1]
        self.__number_of_states = np.shape(state_data)[1]
        self.__KD_tree = KDTree( samples_state_space, leafsize=leafsize )
        self.__state_space_data = state_data

    
    def query( self, query_samples:np.ndarray[float], nearest_neighbors:int=3, inverse_distance_exponent:float=2 ):
        """Interpolate based based on the inverse-distance-weighted values of a specified number of nearest neighbors.

        :param query_samples: query samples
        :type query_samples: np.ndarray[float]
        :param nearest_neighbors: number of nearest neighbors, defaults to 3
        :type nearest_neighbors: int, optional
        :param inverse_distance_exponent: exponent of the inverse distance to the nearest neighbors, defaults to 2
        :type inverse_distance_exponent: float, optional
        :return: interpolated data from nearest neighbors
        :rtype: np.ndarray[float]
        """
        self.__number_of_nearest_neighbors = nearest_neighbors
        self.__inverse_distance_exponent = inverse_distance_exponent

        formatted_query_samples = self.__parse_query_input(query_samples)

        return self.__interpolate_from_inverse_distance(formatted_query_samples)
    
    def __parse_query_input(self, input_query_samples:np.ndarray[float]):
        formatted_query_samples = self.__format_query_input(input_query_samples)
        self.__check_correct_input_dimensions(formatted_query_samples)
        return formatted_query_samples
    
    def __format_query_input(self, unformatted_query_samples:np.ndarray[float]):
        number_of_dimensions = unformatted_query_samples.ndim
        if number_of_dimensions == 1:
            self.__single_query_sample = True
            return unformatted_query_samples.reshape(1, -1)
        else:
            self.__single_query_sample = False
            return unformatted_query_samples
    
    def __check_correct_input_dimensions(self, query_samples:np.ndarray[float]):
        if query_samples.shape[1] != self.__nDim_input:
            raise Exception("Incorrect input dimensionality")
        return
    
    def __interpolate_from_inverse_distance(self, query_input_samples:np.ndarray[float]):
        distances, nearest_neighbor_indices = self.__evaluate_KD_tree(query_input_samples, self.__number_of_nearest_neighbors)
        interpolated_data = self.__interpolate_from_nearest_neighbors(nearest_neighbor_indices, distances)
        if self.__single_query_sample:
            return interpolated_data[0]
        else:
            return interpolated_data
    
    def queryJacobian(self, query_input_samples:np.ndarray[float], nearest_neighbors:int=3, inverse_distance_exponent:float=2):
        
        dx = 1e-2
        grads = np.zeros([self.__nDim_input, query_input_samples.shape[0], self.__state_space_data.shape[1]])
        for i in range(self.__nDim_input):
            q_copy = query_input_samples.copy()
            if query_input_samples.ndim==1:
                q_copy[i] += dx
            else:
                q_copy[:, i] += dx
            interp_data_plus = self.query(q_copy, nearest_neighbors, inverse_distance_exponent)
            if query_input_samples.ndim==1:
                q_copy[i] -= 2*dx
            else:
                q_copy[:, i] -= 2*dx
            interp_data_minus = self.query(q_copy, nearest_neighbors, inverse_distance_exponent)
            
            grads[i] = (interp_data_plus - interp_data_minus) / (2*dx)
        return grads
    
    def __evaluate_KD_tree(self, query_input_samples:np.ndarray[float], number_of_nearest_neighbors:int):
        distances_to_nearest_neighbors, nearest_neighbor_indices = self.__KD_tree.query(query_input_samples, k=number_of_nearest_neighbors,eps=0)
        return distances_to_nearest_neighbors, nearest_neighbor_indices
    
    def __interpolate_from_nearest_neighbors(self, nearest_neighbor_indices:np.ndarray[int], distances_to_neighbors:np.ndarray[float]):
        number_of_query_samples = np.shape(distances_to_neighbors)[0]
        interpolated_state_data = np.zeros([number_of_query_samples, self.__number_of_states])
        for j, (dist, ix) in enumerate(zip( distances_to_neighbors, nearest_neighbor_indices )):
            if self.__single_nearest_neighbor():
                interpolated_state_data[j] = self.__state_space_data[ix]
            elif self.__nearest_sample_is_too_close(dist):
                interpolated_state_data[j] = self.__state_space_data[ix[0]]
            else:  # weight z s by 1/dist --
                nearest_neighbor_coefficients = np.power(dist, -self.__inverse_distance_exponent)
                nearest_neighbor_coefficients /= np.sum(nearest_neighbor_coefficients)
                interpolated_state_data[j] = np.dot( nearest_neighbor_coefficients, self.__state_space_data[ix])

        return interpolated_state_data
    
    def __single_nearest_neighbor(self):
        return self.__number_of_nearest_neighbors==1
    
    def __nearest_sample_is_too_close(self, ordered_distances:np.ndarray[float]):
        return ordered_distances[0] < 1e-10
    
    def query_with_local_parameters( self, q:np.ndarray[float], nnear_vals:np.ndarray[int], p_vals:np.ndarray[float]):
        """Interpolate with number of nearest neighbors and p-factor for each sample.

        :param q: query samples
        :type q: np.ndarray[float]
        :param nnear_vals: number of nearest neighbors for each query sample
        :type nnear_vals: np.array[int]
        :param p_vals: distance exponent factor for each query sample
        :type p_vals: np.array[float]
        :return: interpolated data from nearest neighbors
        :rtype: np.ndarray[float]
        """
            # nnear nearest neighbours of each query point --
        q = np.asarray(q)
        max_nnear = max(max(nnear_vals), 2)
        qdim = q.ndim
        if qdim == 1:
            q = np.expand_dims(q, -1)
            nnear_vals = np.expand_dims(nnear_vals, -1)
            p_vals = np.expand_dims(p_vals, -1)

        eps=0
        self.distances, self.ix = self.__KD_tree.query( q, k=max_nnear, eps=eps )
        interpol = np.zeros([np.shape(q)[0], np.shape(self.z)[1]])
        jinterpol = 0
        for p, nnear, dist, ix in zip(p_vals, nnear_vals, self.distances, self.ix ):
            if nnear == 1 or dist[0] < 1e-10:
                wz = self.z[ix[0]]
            else:  # weight z s by 1/dist --
                w = np.power(dist[:nnear], -p)
                w /= np.sum(w)
                wz = np.dot(w,self.z[ix[:nnear]])
            interpol[jinterpol, :] = wz
            jinterpol += 1
        return interpol if qdim > 1  else interpol[0]
    
class fluidDataInterpolator:

    __controlling_variable_nodes:np.ndarray[float] = None
    __lookup_data:np.ndarray[float] = None
    __lookup_tree:Invdisttree = None
    __n_near_single:int = 6
    __p_fac_single:int = 2
    __state_vars:list[str] = None

    def __init__(self,cv_data:np.ndarray[float], state_data:pd.DataFrame, number_of_nearest_neighbors:int=None, inverse_distance_exponent:float=None):
        """Inverse distance weighted interpolation algorithm for the interpolation of thermodynamic state data from unstructured point clouds.

        :param cv_data: scaled controlling variable data from point cloud samples
        :type cv_data: np.ndarray[float]
        :param state_data: interpolation data from point cloud samples.
        :type state_data: pd.DataFrame
        :param number_of_nearest_neighbors: number of nearest neighbors from which fluid data is interpolated, defaults to None
        :type number_of_nearest_neighbors: int, optional
        :param inverse_distance_exponent: exponent applied to the inverse distance of the nearest neighbors, defaults to None
        :type inverse_distance_exponent: float, optional
        """
        self.__controlling_variable_nodes = cv_data
        self.__lookup_data = state_data.values
        self.__state_vars = list(state_data.keys())
        self.__lookup_tree = Invdisttree(cv_data, state_data.values)
        self.__n_near_single = number_of_nearest_neighbors
        self.__p_fac_single = inverse_distance_exponent

        if not self.__n_near_single or not self.__p_fac_single:
            self.tuneTreeParameters()
        return

    def getControllingVariableNodes(self):
        """Scaled controlling variable values of the point cloud samples the interpolator was built from.

        :return: scaled controlling variable data.
        :rtype: np.ndarray[float]
        """
        return self.__controlling_variable_nodes

    def tuneTreeParameters(self):

        Np_total = np.shape(self.__controlling_variable_nodes)[0]
        Np_train = int(0.8*Np_total)

        n_control_vars = np.shape(self.__controlling_variable_nodes)[1]
        scaler = MinMaxScaler()
        lookup_data_scaled = scaler.fit_transform(self.__lookup_data)
        full_data_set = np.column_stack((self.__controlling_variable_nodes, lookup_data_scaled))
        np.random.shuffle(full_data_set)

        train_data = full_data_set[:Np_train]
        test_data = full_data_set[Np_train:]

        cv_data_train= train_data[:, :n_control_vars]
        cv_data_test= test_data[:, :n_control_vars]

        train_tree = Invdisttree(cv_data_train, cv_data_train)
        n_near_range = range(1, 20)
        p_range = np.linspace(0, 6, 10)
        RMS_ppv = np.zeros([len(n_near_range), len(p_range)])
        for i in tqdm(range(len(n_near_range)),desc="Searching for optimal tree parameters..."):
            for j in range(len(p_range)):
                interpolated_state = train_tree.query(cv_data_test, nearest_neighbors=n_near_range[i], inverse_distance_exponent=p_range[j])
                rms_local = np.sum(np.power(interpolated_state - cv_data_test, 2))
                RMS_ppv[i,j] = rms_local
        [imin,jmin] = divmod(RMS_ppv.argmin(), RMS_ppv.shape[1])
        self.__n_near_single = n_near_range[imin]
        self.__p_fac_single = p_range[jmin]

        print("Number of nearest neighbors: %i" % self.__n_near_single)
        print("Inverse distance exponent: %.2f" % self.__p_fac_single)
        return

    def __call__(self, cv_data_query:np.ndarray[float],varnames:list[str]=[]):
        ndim_input = cv_data_query.ndim
        if ndim_input < 2:
            cv_data_query = cv_data_query[np.newaxis]

        fluid_data_interp = self.__lookup_tree.query(cv_data_query, nearest_neighbors=self.__n_near_single, inverse_distance_exponent=self.__p_fac_single)
        nvars =len(varnames)
        if nvars > 0:
            var_indices = [self.__state_vars.index(var) for var in varnames]

            if fluid_data_interp.shape[0] == 1:
                return fluid_data_interp[0,var_indices]
            else:
                return fluid_data_interp[:, var_indices]
        else:
            if fluid_data_interp.shape[0] == 1:
                return fluid_data_interp[0]
            else:
                return fluid_data_interp

    def Jacobian(self, cv_data_query:np.ndarray[float],varnames:list[str]=[]):
        ndim_input = cv_data_query.ndim
        if ndim_input < 2:
            cv_data_query = cv_data_query[np.newaxis]
        
        grads = self.__lookup_tree.queryJacobian(cv_data_query, nearest_neighbors=self.__n_near_single, inverse_distance_exponent=self.__p_fac_single)
        nvars =len(varnames)
        if nvars > 0:
            var_indices = [self.__state_vars.index(var) for var in varnames]

            if grads.shape[1] == 1:
                return grads[:,0,var_indices]
            else:
                return grads[:,:,var_indices]
        else:
            if grads.shape[1] == 1:
                return grads[:,0,:]
            else:
                return grads