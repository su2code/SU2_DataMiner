import numpy as np
import gmsh
from concave_hull import concave_hull, concave_hull_indexes
from collections.abc import Callable
from Common.CommonMethods import shoelace, FiniteDifferenceDerivative
from shapely.geometry import Point, Polygon

def default_refinement_function(this, x:float, y:float, z:float):
    return 1.0

class Mesh2DPlane:

    __pointCloud:np.ndarray[float] = None
    _pointCloud_hullNodes:np.ndarray[float] = None
    __mesh_along_coords:list[int] = [0, 1]

    _gmsh_geo:gmsh.model.geo = None
    _gmsh_mesher:gmsh.model.mesh = None
    _gmsh_verbosity:int = 0
    _mesher_verbosity:int = 1

    _base_cell_size:float = 2e-2

    refinement_function = default_refinement_function

    __plane_area:float = None
    __perimiter_tag:int=None
    __use_target_node_count:bool = False
    __target_node_count:int=3000

    __meshNodes:np.ndarray[float] = None
    __connectivity:np.ndarray[int] = None
    __hullNodeIDs:np.ndarray[int] = None


    def __init__(self):
        return
    
    def generateMesh(self):
        self.__checkPlanarCoords()

        self._initializeGmesh()
        self.__createGeometry()

        if self.__use_target_node_count:
            self.__optimizeMeshResolution()
            
        self.__meshNodes, self.__connectivity, self.__hullNodeIDs = self.__mesh2D()
        self._gmsh_geo.synchronize()
        return
    
    def __optimizeMeshResolution(self):
        self.__printMessage("Refining table based on target node count")
        self.__printMessage("Target node count: %i" % self.__target_node_count)

        sufficient_refinement = False
        niter_max = 20
        relaxation = 0.35

        # Initial guess
        self._base_cell_size = 2 * np.sqrt(self.__plane_area / self.__target_node_count)
        iter = 0
        self.__printMessage("| Iteration | # Nodes | Cell size | Relative % diff |")
        while not sufficient_refinement and iter < niter_max:
            nodes, _, _ = self.__mesh2D()
            n_nodes = len(nodes)
            rel_diff = abs(float(n_nodes - self.__target_node_count)/self.__target_node_count)

            self.__printMessage("| %i | %i | %.3e | %.1f |" % (iter, n_nodes, self._base_cell_size, 100*rel_diff))
            # Terminate if relative difference is less than 1%
            if rel_diff > 0.01:
                self._base_cell_size *= (1 + relaxation * (float(n_nodes / self.__target_node_count) - 1))
            else:
                sufficient_refinement = True
            iter += 1
        
        if iter == niter_max:
            self.__printMessage("Target node count was not reached within %i iterations" % niter_max)
        return
    
    def __createGeometry(self):
        

        curvloop_table_outline, perimiter_tags = self.__createHullPerimiter()
        self._gmsh_geo.synchronize()
        self._createInternalGeometry(curvloop_table_outline, perimiter_tags)
        return
    
    def __mesh2D(self):
        
        self.__prepareMesher()

        self._gmsh_mesher.generate(2)

        if self._mesher_verbosity > 1:
            gmsh.fltk.run()

        meshNodeCoords, triaTags = self.__extractNodesTriangles()

        hullTags = self.__extractHulltags()
        return meshNodeCoords, triaTags, hullTags
    
    def __checkPlanarCoords(self):
        all_dims = [0, 1, 2]
        for iDim in self.__mesh_along_coords:
            all_dims.remove(iDim)
        const_dim = all_dims[0]
        equal_Z_coords = np.all(self.__pointCloud[:, const_dim]==self.__pointCloud[0, const_dim])
        if not equal_Z_coords:
            raise Exception("Point cloud does not contain planar coordinates")
        return
    
    def __createHullPerimiter(self):
        self._pointCloud_hullNodes = self.__findHullNodes()
        planeCurvLoop, perimiterTags = self.__createHullCurvLoop(self._pointCloud_hullNodes)
        self.__perimiter_tag = self._gmsh_geo.addPhysicalGroup(1, perimiterTags, name="perimeter")
        return planeCurvLoop, perimiterTags
    

    def __findHullNodes(self):
        
        planarCoords_PointCloud = self.__pointCloud[:, self.__mesh_along_coords]
        hullIndices = concave_hull_indexes(planarCoords_PointCloud, length_threshold=self._base_cell_size)

        planarCoords_hull = planarCoords_PointCloud[hullIndices]
        self.__plane_area = shoelace(planarCoords_hull)

        return self.__pointCloud[hullIndices]
    
    def __createHullCurvLoop(self, hull_coords:np.ndarray[float]):

        hull_pts = self.__createHullPointEntities(hull_coords)
        hull_perimiter = self.__connectHullPoints(hull_pts)

        curvloop = self._gmsh_geo.addCurveLoop(hull_perimiter)

        self._gmsh_geo.synchronize()

        return curvloop, hull_perimiter
  
    def __createHullPointEntities(self, coords_hull_pts:np.ndarray[float]):
        hull_pts = []
        for i in range(int(len(coords_hull_pts))):
            hull_pts.append(self._gmsh_geo.addPoint(coords_hull_pts[i, 0], coords_hull_pts[i, 1], coords_hull_pts[i, 2]))
        return hull_pts
    
    def __connectHullPoints(self, pointIDs:list[int]):
        hull_lines = []
        for i in range(len(pointIDs)):
            j = (i + 1) % len(pointIDs)
            hull_lines.append(self._gmsh_geo.addLine(pointIDs[i], pointIDs[j]))
        return hull_lines
    
    def setInitialPointCloud(self, point_could_3D:np.ndarray[float]):
        """Provide the coordinates of the point cloud which space is discretized.

        :param point_could_3D: normalized coordinates of the initial point cloud.
        :type point_could_3D: np.ndarray[float]
        """
        self.__pointCloud = self.__filterNansFromPointCloud(point_could_3D)

        return
    
    def __filterNansFromPointCloud(self, point_cloud_3D:np.ndarray[float]):
        unique_pts_init = np.unique(point_cloud_3D, axis=0)
        nans = np.isnan(unique_pts_init).any(1)
        pts_wo_nans = unique_pts_init[np.invert(nans)]
        return pts_wo_nans

    
    def _createInternalGeometry(self,curvloop_perimiter:int, perimiter_tags:list[int]):

        planeID = self._gmsh_geo.addPlaneSurface([curvloop_perimiter])
        self._gmsh_geo.addPhysicalGroup(2, [planeID])
        return
    
    def __extractNodesTriangles(self):
        nodeTags, coords, _ = self._gmsh_mesher.getNodes()
        nodeTags = np.asarray(nodeTags, dtype=np.int64)
        MeshPoints = np.asarray(coords, dtype=float).reshape(-1, 3)

        order = np.argsort(nodeTags)
        nodeTags_sorted = nodeTags[order]

        # 2) 2D elements
        elemTypes, _, elemNodeTags = self._gmsh_mesher.getElements(2)
       
        tris = []
        quads = []

        for et, nodes_flat in zip(elemTypes, elemNodeTags):
            if et == 2:  # triangles with 3 nodes
                tri_tags = np.asarray(nodes_flat, dtype=np.int64).reshape(-1, 3)
                tris.append(self.__map_tags(tri_tags, nodeTags_sorted,order).reshape(-1, 3))
            elif et == 3:  # quad with 4 nodes
                quad_tags = np.asarray(nodes_flat, dtype=np.int64).reshape(-1, 4)
                quads.append(self.__map_tags(quad_tags, nodeTags_sorted,order).reshape(-1, 4))

        tris = np.vstack(tris) if tris else np.zeros((0, 3), dtype=np.int64)

        if quads:
            quads = np.vstack(quads)
            tris = np.vstack([
                tris,
                quads[:, [0, 1, 2]],
                quads[:, [0, 2, 3]],
            ])
        return MeshPoints, tris
    
    def __map_tags(self,tags,nodeTags_sorted,order):
        tags = np.asarray(tags, dtype=np.int64).ravel()
        pos = np.searchsorted(nodeTags_sorted, tags)
        return order[pos]


    def __extractHulltags(self):
        perimiter_curve_tags = gmsh.model.getEntitiesForPhysicalGroup(1, self.__perimiter_tag)
        allHullTags = []
        for line in perimiter_curve_tags:
            nodeTags, nodes, _ = self._gmsh_mesher.getNodes(includeBoundary=True,tag=line,dim=1)
            nodes = np.asarray(nodes, dtype=float).reshape(-1, 3)
            nodeTags = np.asarray(nodeTags, dtype=np.int64)

            allHullTags = np.append(allHullTags, nodeTags)
        hullTags = np.unique(np.int64(allHullTags))
        return hullTags
    
    def setBaseCellSize(self, cell_size_input:float):
        """Provide the size metric for the coarse cells in the 2D mesh.

        :param cell_size_input: coarse cell size
        :type cell_size_input: float
        :raises Exception: if the input is not strictly positive.
        """
        if cell_size_input <= 0:
            raise Exception("Cell size should be strictly positive")
        self._base_cell_size = cell_size_input
        return
    
    def _initializeGmesh(self):
        gmsh.initialize()
        gmsh.model.add("2D plane")
        gmsh.option.setNumber("General.Verbosity", self._gmsh_verbosity)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints",0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary",0)

        self._gmsh_geo = gmsh.model.geo
        self._gmsh_mesher = gmsh.model.mesh
        return
    
    def __prepareMesher(self):
        self._gmsh_geo.synchronize()
        self._gmsh_mesher.clear()
        gmsh.option.setNumber("Mesh.MeshSizeMax", self._base_cell_size)

        def meshSizeCallback(dim,tag,x,y,z,lc):
            refinement_factor = self.refinement_function(x, y, z)
            return self._base_cell_size * refinement_factor
        
        self._gmsh_mesher.setSizeCallback(meshSizeCallback)
        return
    
    def setRefinementFunction(self, function_input:Callable):
        """Provide a function which dictates the local refinement in the 2D table.

        :param function_input: Function that takes the normalized 3D mesh coordinates (x, y, z) as input and outputs the mesh refinement factor (0.0-1.0).
        :type function_input: Callable
        """
        self.refinement_function = function_input
        return
    
    def setTargetNodeCount(self, Np_target:int=3000):
        """Specify the target number of nodes in the 2D table. The number of nodes will be approximated within 1%.

        :param Np_target: targeted number of nodes, defaults to 3000
        :type Np_target: int, optional
        :raises Exception: if the specified number is not strictly positive.
        """
        if Np_target <= 0:
            raise Exception("Target number of nodes should be positive")
        self.__target_node_count = Np_target
        self.__use_target_node_count = True
        return
    
    def getMeshNodes(self):
        """Returns the 3D mesh node locations.

        :return: mesh node coordinates.
        :rtype: np.ndarray[float]
        """
        return self.__meshNodes
    
    def getConnectivity(self):
        """Returns the connectivity information of the 2D table.

        :return: mesh node connectivity information.
        :rtype: np.ndarray[int]
        """
        return self.__connectivity
    
    def getHullIDs(self):
        """Returns the node indices of the perimiter of the 2D table.

        :return: node indices of the perimiter.
        :rtype: np.ndarray[int]
        """
        return self.__hullNodeIDs

    def __printMessage(self, msg:str):
        if self._mesher_verbosity > 0:
            print(msg)
        return
    
    def setVerbosity(self, verbosity_level:int=1):
        """Specify verbosity level of the table generator. No information will be displayed with a value 0. For a value 1, information on the mesh generation process are displayed in the terminal. For higher values, the 2D mesh is visualized through the Gmsh interface.

        :param verbosity_level: verbosity level, defaults to 1
        :type verbosity_level: int, optional
        """
        self._mesher_verbosity = verbosity_level
        if self._mesher_verbosity > 1:
            self._gmsh_verbosity = 4

        return


class MeshThermodynamicPlane(Mesh2DPlane):

    __includeSaturationCurve:bool=False
    __saturation_curve_points:np.ndarray[float] = None

    def __init__(self):
        super().__init__()
        return
    
    def _createInternalGeometry(self, perimiter_curvloop:int, perimiter_tags:list[int]):
        if self.__includeSaturationCurve:
            saturation_curve_crvloops = self.__createSaturationCurveGeom()
            fluid_plane_crvloops = [perimiter_curvloop] + [-c for c in saturation_curve_crvloops]
            fluid_surf = self._gmsh_geo.addPlaneSurface(fluid_plane_crvloops)
            sat_surfs = [self._gmsh_geo.addPlaneSurface([c]) for c in saturation_curve_crvloops]
            self._gmsh_geo.addPhysicalGroup(2, [fluid_surf] + sat_surfs)
        else:
            super()._createInternalGeometry(perimiter_curvloop, perimiter_tags)
        return
    
    def __createSaturationCurveGeom(self):
        
        sat_curve_pts, sat_curve_lower, sat_curve_upper = self.__saturationOffsetCurves()
        if len(sat_curve_lower) > 0:
            sat_curve_lower_pts, sat_curve_upper_pts, segment_lengths = self.__samplePointsAlongSaturationCurve(sat_curve_pts, sat_curve_lower, sat_curve_upper)

            return self.__encloseSaturationCurve(sat_curve_lower_pts, sat_curve_upper_pts, segment_lengths)
        else:
            return []
        
    def __encloseSaturationCurve(self, lower_offset_pts:np.ndarray[int], upper_offset_pts:np.ndarray[int], segment_lengths:list[float]):
        connecting_lines = []
        Npoints = len(lower_offset_pts)
        for i in range(Npoints):
            connecting_lines.append(self._gmsh_geo.addLine(lower_offset_pts[i], upper_offset_pts[i]))

        sat_curve_crvloops = []
        for i in range(Npoints-1):
            if segment_lengths[i] < 2*self._base_cell_size:
                tangent_line_upper=self._gmsh_geo.addLine(upper_offset_pts[i],upper_offset_pts[i+1])
                tangent_line_lower=self._gmsh_geo.addLine(lower_offset_pts[i],lower_offset_pts[i+1])

                sat_curve_crvloops.append(self._gmsh_geo.addCurveLoop([tangent_line_upper, connecting_lines[i+1], tangent_line_lower, connecting_lines[i]], reorient=True))
        self._gmsh_geo.synchronize()
        return sat_curve_crvloops
    
    def __samplePointsAlongSaturationCurve(self, saturation_curve_points:np.ndarray[float], offset_lower:np.ndarray[float], offset_upper:np.ndarray[float]):
        i = 0
        j = 1
        sat_curve_upper_pts = []
        sat_curve_lower_pts = []
        segment_lengths = []
        sat_curve_upper_pts.append(self._gmsh_geo.addPoint(offset_upper[0,0], offset_upper[0,1],0))
        sat_curve_lower_pts.append(self._gmsh_geo.addPoint(offset_lower[0,0], offset_lower[0,1],0))
        while j < len(saturation_curve_points):
            dist = np.sqrt(np.sum(np.power(saturation_curve_points[j] - saturation_curve_points[i],2)))
            local_refinement_factor = self.refinement_function(saturation_curve_points[j, 0], saturation_curve_points[j, 1], 0)
            if dist < 0.5*local_refinement_factor*self._base_cell_size:
                j += 1
            else:
                i = j
                j += 1
                segment_lengths.append(dist)
                sat_curve_upper_pts.append(self._gmsh_geo.addPoint(offset_upper[i, 0], offset_upper[i,1],0))
                sat_curve_lower_pts.append(self._gmsh_geo.addPoint(offset_lower[i, 0], offset_lower[i,1],0))
        return sat_curve_lower_pts, sat_curve_upper_pts, segment_lengths
    
    def __saturationOffsetCurves(self):

        norm_vector = self.__getSaturationCurveNormal()

        # Create offset curves to ensure that no nodes are generated on the saturation curve itself.
        sat_curve_rhoe_lower, sat_curve_rhoe_upper = self.__createSaturationCurveOffsets(norm_vector)

        valid_pts = self.__clipSaturationCurveToPlane(sat_curve_rhoe_lower, sat_curve_rhoe_upper)
        sat_curve_pts = self.__saturation_curve_points[valid_pts, :]
        norm_vector = norm_vector[valid_pts, :]
        sat_curve_rhoe_lower = sat_curve_rhoe_lower[valid_pts]
        sat_curve_rhoe_upper = sat_curve_rhoe_upper[valid_pts]

        return sat_curve_pts, sat_curve_rhoe_lower, sat_curve_rhoe_upper
    
    def __getSaturationCurveNormal(self):
        dedrho_sat_norm = FiniteDifferenceDerivative(self.__saturation_curve_points[:,1], self.__saturation_curve_points[:,0])
        norm_vector = np.column_stack((-dedrho_sat_norm, np.ones(len(dedrho_sat_norm))))
        norm_vector = norm_vector / np.sqrt(np.sum(np.power(norm_vector, 2), axis=1))[:,np.newaxis]
        
        return norm_vector
    
    def __createSaturationCurveOffsets(self, normal_vector:np.ndarray[float]):
        offset = 5e-4*self._base_cell_size

        sat_curve_rhoe_upper = self.__saturation_curve_points + offset * normal_vector
        sat_curve_rhoe_lower = self.__saturation_curve_points - offset * normal_vector
        
        return sat_curve_rhoe_lower, sat_curve_rhoe_upper
    

    def __clipSaturationCurveToPlane(self, pts_offset_lower:np.ndarray[float], pts_offset_upper:np.ndarray[float]):

        hull_pts_orig = concave_hull(self._pointCloud_hullNodes, length_threshold=self._base_cell_size)
        polygon_thermodynamic_plane = Polygon(hull_pts_orig)
        within_hull = np.zeros(len(self.__saturation_curve_points),dtype=bool)

        for i in range(len(self.__saturation_curve_points)):
            pt_upper = Point(pts_offset_upper[i])
            pt_lower = Point(pts_offset_lower[i])
            within_hull[i] = polygon_thermodynamic_plane.contains(pt_upper) and polygon_thermodynamic_plane.contains(pt_lower)

        return within_hull
    
    def setSaturationCurvePoints(self, saturation_curve_rhoe:np.ndarray[float]):
        self.__saturation_curve_points = saturation_curve_rhoe
        self.__includeSaturationCurve = True
        return
    
    