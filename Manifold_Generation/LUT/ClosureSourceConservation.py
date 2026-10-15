###############################################################################################
#       #      _____ __  _____      ____        __        __  ____                   #        #
#       #     / ___// / / /__ \    / __ \____ _/ /_____ _/  |/  (_)___  ___  _____   #        #
#       #     \__ \/ / / /__/ /   / / / / __ `/ __/ __ `/ /|_/ / / __ \/ _ \/ ___/   #        #
#       #    ___/ / /_/ // __/   / /_/ / /_/ / /_/ /_/ / /  / / / / / /  __/ /       #        #
#       #   /____/\____//____/  /_____/\__,_/\__/\__,_/_/  /_/_/_/ /_/\___/_/        #        #
#       #                                                                            #        #
###############################################################################################

##################### FILE NAME: ClosureSourceConservation.py #################################
#=============================================================================================#
# Conservation audit and integral-preserving correction for the preferential-diffusion       |
# closure sources S_PV, S_Z, S_h (Eq. 16 of Schepers & van Oijen, C&F 280 (2025), FGM B2)    |
# after interpolation onto the (PV, h, Z) table levels.                                      |
#                                                                                            |
# The closure sources are divergences of diffusion fluxes along 1D flamelets: their          |
# integral over every flamelet vanishes. The table reconstruction does NOT preserve this —   |
# the sources are stiff dipoles concentrated in the flame front, and sampling them at table  |
# nodes and interpolating back leaks front-magnitude values into the slow preheat and        |
# post-flame regions, leaving a spurious net source of order a percent of the inlet          |
# convective flux at every resolved flame.                                                   |
#                                                                                            |
# This module (i) AUDITS the tabulated sources by integrating the reconstructed field along  |
# every 1D flamelet trajectory, replicating CLookUpTable::LookUp_XYZ of the solver, and      |
# (ii) PROJECTS the nodal values so those integrals match the flamelet's own, via a          |
# magnitude-weighted minimum-norm least-squares correction — the smallest weighted           |
# perturbation of the tabulated field that restores conservation along all flamelets at      |
# once. Conservation is exact along the fitted flamelets; states far from any flamelet keep  |
# a reconstruction error set by how well the table resolves the source profile.              |
#=============================================================================================#

import numpy as np
from os import listdir, path

from Common.Properties import FGMVars, PDClosureSourceName

_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def _collect_flamelet_files(flamelet_data_dir: str, mix_status_header: str = "phi_"):
    """List all free-flame and burner-flame flamelet CSV files below the flamelet data directory."""
    files = []
    for sub in ("freeflame_data", "burnerflame_data"):
        base = path.join(flamelet_data_dir, sub)
        if not path.isdir(base):
            continue
        for mix_folder in sorted(listdir(base)):
            mix_path = path.join(base, mix_folder)
            if not path.isdir(mix_path) or not mix_folder.startswith(mix_status_header):
                continue
            for fname in sorted(listdir(mix_path)):
                if fname.endswith(".csv"):
                    files.append(path.join(mix_path, fname))
    return files


def _load_trajectory(csv_file: str, pv_species: list, pv_weights: list):
    """Load the (PV, h, Z) trajectory, quadrature weights, and flux scales of a 1D flamelet.

    Returns None for degenerate (single-row or non-burning) files.
    """
    with open(csv_file, "r") as fid:
        variables = fid.readline().strip().split(",")
    D = np.loadtxt(csv_file, delimiter=",", skiprows=1)
    if D.ndim != 2 or np.shape(D)[0] < 3:
        return None
    iv = {v: i for i, v in enumerate(variables)}
    x = D[:, iv[FGMVars.Distance.name]]
    Z = D[:, iv[FGMVars.MixtureFraction.name]]
    h = D[:, iv["EnthalpyTot"]]
    T = D[:, iv[FGMVars.Temperature.name]]
    rho = D[:, iv[FGMVars.Density.name]]
    cp = D[:, iv[FGMVars.Cp.name]]
    u = D[:, iv["Velocity"]]
    pv = np.zeros(len(x))
    for sp, w in zip(pv_species, pv_weights):
        pv += w * D[:, iv["Y-" + sp]]

    # Trapezoidal quadrature weights along the physical coordinate.
    dx = np.diff(x)
    w_quad = np.zeros(len(x))
    w_quad[0] = 0.5 * dx[0]
    w_quad[-1] = 0.5 * dx[-1]
    w_quad[1:-1] = 0.5 * (dx[:-1] + dx[1:])

    mass_flux = rho[0] * u[0]
    scales = {
        PDClosureSourceName(FGMVars.ProgressVariable.name): mass_flux * max(np.ptp(pv), 1e-32),
        PDClosureSourceName(FGMVars.MixtureFraction.name): mass_flux * max(abs(np.mean(Z)), np.ptp(Z), 1e-32),
        PDClosureSourceName(FGMVars.EnthalpyTot.name): mass_flux * max(np.ptp(h), cp[0] * np.ptp(T), 1e-32),
    }
    stored = {}
    for var in (PDClosureSourceName(FGMVars.ProgressVariable.name),
                PDClosureSourceName(FGMVars.MixtureFraction.name),
                PDClosureSourceName(FGMVars.EnthalpyTot.name)):
        if var in iv:
            stored[var] = D[:, iv[var]]
    return {"file": csv_file, "pv": pv, "h": h, "Z": Z, "w_quad": w_quad,
            "scales": scales, "stored": stored}


class ClosureSourceConservation:
    """Audit and enforce flamelet-wise conservation of tabulated closure sources.

    The reconstruction replicated here follows CLookUpTable::LookUp_XYZ of the SU2 flamelet
    LUT: the progress variable is re-normalized to per-level coordinates (enthalpy is passed
    raw), each bracketing level is interpolated barycentrically on its stored triangulation
    with an inverse-distance fallback on the two nearest hull nodes, the two results are
    blended linearly in mixture fraction, and the closure sources are dropped entirely at
    states the solver flags as a miss.
    """

    def __init__(self, z_levels: np.ndarray, level_nodes: list, level_connectivity: list,
                 level_hulls: list = None):
        """
        :param z_levels: mixture fraction values of the table levels.
        :param level_nodes: per level, (n_nodes, 2) array of (PV, EnthalpyTot) node coordinates.
        :param level_connectivity: per level, (n_tria, 3) zero-based triangle node indices.
        :param level_hulls: per level, indices of the hull nodes. Required to replicate the
            solver's nearest-neighbour fallback for queries outside a level's triangulation.
        """
        from matplotlib.tri import Triangulation
        self._z_levels = np.asarray(z_levels)
        self._level_nodes = [np.asarray(n, dtype=float) for n in level_nodes]
        self._node_offset = np.concatenate([[0], np.cumsum([len(n) for n in self._level_nodes])])
        self._n_nodes_total = int(self._node_offset[-1])
        self._tri, self._trifinder = [], []
        for nodes, conn in zip(self._level_nodes, level_connectivity):
            tri = Triangulation(nodes[:, 0], nodes[:, 1], np.asarray(conn, dtype=int))
            self._tri.append(tri)
            self._trifinder.append(tri.get_trifinder())
        self._x_lim = np.array([[n[:, 0].min(), n[:, 0].max()] for n in self._level_nodes])
        self._y_lim = np.array([[n[:, 1].min(), n[:, 1].max()] for n in self._level_nodes])
        if level_hulls is None:
            self._hulls = [np.arange(len(n)) for n in self._level_nodes]
        else:
            self._hulls = [np.asarray(hh, dtype=int) for hh in level_hulls]

    def _lookup_xy(self, iLevel: int, pv: np.ndarray, h: np.ndarray):
        """Replicate CLookUpTable::LookUp_XY on one level.

        Inside the triangulation the value is the barycentric interpolant of the containing
        triangle; outside, the solver falls back to inverse-distance weighting of the two
        nearest HULL nodes in per-level normalized coordinates (it does not return zero).

        :return: (inside flags, node ids (m,3), weights (m,3))
        """
        tri = self._tri[iLevel]
        tids = self._trifinder[iLevel](pv, h)
        inside = tids >= 0
        nid = np.zeros([len(pv), 3], dtype=int)
        w = np.zeros([len(pv), 3])
        if np.any(inside):
            t = tri.triangles[tids[inside]]
            x1, y1 = tri.x[t[:, 0]], tri.y[t[:, 0]]
            x2, y2 = tri.x[t[:, 1]], tri.y[t[:, 1]]
            x3, y3 = tri.x[t[:, 2]], tri.y[t[:, 2]]
            det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
            det = np.where(np.abs(det) < 1e-300, 1e-300, det)
            l1 = ((y2 - y3) * (pv[inside] - x3) + (x3 - x2) * (h[inside] - y3)) / det
            l2 = ((y3 - y1) * (pv[inside] - x3) + (x1 - x3) * (h[inside] - y3)) / det
            nid[inside] = t
            w[inside] = np.column_stack([l1, l2, 1.0 - l1 - l2])
        out = ~inside
        if np.any(out):
            hull = self._hulls[iLevel]
            nodes = self._level_nodes[iLevel]
            cx = 1.0 / (self._x_lim[iLevel, 1] - self._x_lim[iLevel, 0])
            cy = 1.0 / (self._y_lim[iLevel, 1] - self._y_lim[iLevel, 0])
            hx, hy = nodes[hull, 0] * cx, nodes[hull, 1] * cy
            qx, qy = pv[out] * cx, h[out] * cy
            d2 = (qx[:, None] - hx[None, :]) ** 2 + (qy[:, None] - hy[None, :]) ** 2
            two = np.argpartition(d2, 1, axis=1)[:, :2]
            d = np.sqrt(np.take_along_axis(d2, two, axis=1))
            d = np.maximum(d, 1e-300)
            ww = (1.0 / d) / np.sum(1.0 / d, axis=1, keepdims=True)
            nid[out, 0] = hull[two[:, 0]]
            nid[out, 1] = hull[two[:, 1]]
            nid[out, 2] = hull[two[:, 0]]
            w[out, 0] = ww[:, 0]
            w[out, 1] = ww[:, 1]
            w[out, 2] = 0.0
        return inside, nid, w

    def BuildTrajectoryRow(self, traj: dict):
        """Accumulate the linear functional A_t such that A_t . theta = int S(x) dx along the
        trajectory for nodal values theta (global node numbering across levels).

        Trajectory points whose mixture fraction lies outside the tabulated range contribute
        nothing: those states are off-manifold, the SU2 lookup reports a miss there and the
        closure sources are not applied. Constraining the fit with clamped off-table states
        would distort the nodal values to satisfy conditions the CFD never evaluates.

        :return: (node_ids, weights) of the aggregated sparse row.
        """
        nid_all, w_all, valid = self.ReconstructionWeights(traj)
        w_quad = np.where(valid, traj["w_quad"], 0.0)
        cols, wgts = [], []
        for nid, w in zip(nid_all, w_all):
            keep = w_quad > 0
            if np.any(keep):
                cols.append(nid[keep].ravel())
                wgts.append((w_quad[keep][:, None] * w[keep]).ravel())
        if len(cols) == 0:
            return np.zeros(0, dtype=int), np.zeros(0)
        cols = np.concatenate(cols)
        wgts = np.concatenate(wgts)
        # Aggregate duplicate node contributions (a node can be hit by both levels).
        uniq, inv = np.unique(cols, return_inverse=True)
        agg = np.zeros(len(uniq))
        np.add.at(agg, inv, wgts)
        return uniq, agg

    def ReconstructionWeights(self, traj: dict, z_weight_tol: float = 1e-3,
                              bbox_rtol: float = 1e-3):
        """Replicate CLookUpTable::LookUp_XYZ for every point of a flamelet trajectory.

        Mirrors the solver exactly: the progress variable is re-normalized to per-level
        coordinates (the enthalpy is passed raw to both levels, being a conserved variable),
        queries that fall numerically onto one level are resolved by a single 2D lookup there,
        and a query counts as on-manifold when EITHER bracketing cross-section contains it and
        the progress variable lies within the interpolated bounding box. Off-manifold points
        are where the solver reports a miss and skips the closure sources entirely.

        :return: (list of (n_pts,3) node id arrays, list of (n_pts,3) weight arrays, valid mask)
        """
        pv, h, z = traj["pv"], traj["h"], traj["Z"]
        n = len(pv)
        zl = self._z_levels
        within_z = (z >= zl[0]) & (z <= zl[-1])

        iup = np.clip(np.searchsorted(zl, z), 1, len(zl) - 1)
        ilo = iup - 1
        dz = zl[iup] - zl[ilo]
        z_frac = np.where(within_z, (z - zl[ilo]) / dz, 0.0)

        snap = within_z & ((z_frac <= z_weight_tol) | (z_frac >= 1.0 - z_weight_tol))
        blend = within_z & ~snap

        nid_lo = np.zeros([n, 3], dtype=int); w_lo = np.zeros([n, 3])
        nid_up = np.zeros([n, 3], dtype=int); w_up = np.zeros([n, 3])
        valid = np.zeros(n, dtype=bool)

        # --- queries resolved on a single (snap) level ---
        if np.any(snap):
            snap_level = np.where(z_frac >= 1.0 - z_weight_tol, iup, ilo)
            for iL in np.unique(snap_level[snap]):
                sel = np.where(snap & (snap_level == iL))[0]
                x_min, x_max = self._x_lim[iL]
                x_tol = bbox_rtol * (x_max - x_min)
                x_eps = 1e-9 * (x_max - x_min)
                pv_q = pv[sel].copy()
                inband = (pv_q >= x_min - x_tol) & (pv_q <= x_max + x_tol)
                pv_q[inband] = np.clip(pv_q[inband], x_min + x_eps, x_max - x_eps)
                ins, nid, w = self._lookup_xy(iL, pv_q, h[sel])
                nid_lo[sel], w_lo[sel] = nid + self._node_offset[iL], w
                valid[sel] = ins

        # --- queries blended between two levels ---
        if np.any(blend):
            x_min_q = self._x_lim[ilo, 0] + (self._x_lim[iup, 0] - self._x_lim[ilo, 0]) * z_frac
            x_max_q = self._x_lim[ilo, 1] + (self._x_lim[iup, 1] - self._x_lim[ilo, 1]) * z_frac
            bbox_tol = bbox_rtol * (x_max_q - x_min_q)
            in_bbox = (pv >= x_min_q - bbox_tol) & (pv <= x_max_q + bbox_tol)
            x_norm = np.clip((pv - x_min_q) / (x_max_q - x_min_q), 1e-9, 1.0 - 1e-9)

            ins_lo = np.zeros(n, dtype=bool)
            ins_up = np.zeros(n, dtype=bool)
            for arr, nid_t, w_t, ins_t in ((ilo, nid_lo, w_lo, ins_lo), (iup, nid_up, w_up, ins_up)):
                for iL in np.unique(arr[blend]):
                    sel = np.where(blend & (arr == iL))[0]
                    pv_lv = self._x_lim[iL, 0] + (self._x_lim[iL, 1] - self._x_lim[iL, 0]) * x_norm[sel]
                    ins, nid, w = self._lookup_xy(iL, pv_lv, h[sel])
                    nid_t[sel], w_t[sel] = nid + self._node_offset[iL], w
                    ins_t[sel] = ins
            valid[blend] = (ins_lo[blend] | ins_up[blend]) & in_bbox[blend]
            w_lo[blend] *= (1.0 - z_frac[blend])[:, None]
            w_up[blend] *= z_frac[blend][:, None]

        w_lo[~valid] = 0.0
        w_up[~valid] = 0.0
        return [nid_lo, nid_up], [w_lo, w_up], valid

    def SelectOnTableTrajectories(self, trajectories: list, min_coverage: float = 0.25):
        """Keep only flamelets with a substantial portion inside the tabulated mixture fraction
        range. Fully off-table flamelets (e.g. ultra-lean ones below the Z floor) carry no
        information about the tabulated field and must not constrain the fit.

        :param min_coverage: minimum fraction of the flamelet length inside the Z range.
        :return: (kept trajectories, kept index array)
        """
        keep = []
        for it, t in enumerate(trajectories):
            z = t["Z"]
            on = (z >= self._z_levels[0]) & (z <= self._z_levels[-1])
            if np.sum(t["w_quad"][on]) > min_coverage * np.sum(t["w_quad"]):
                keep.append(it)
        keep = np.array(keep, dtype=int)
        return [trajectories[i] for i in keep], keep

    def BuildOperator(self, trajectories: list, verbose: bool = True):
        """Assemble the sparse trajectory-integral operator A (n_traj x n_nodes_total)."""
        from scipy import sparse
        rows, cols, vals = [], [], []
        for it, traj in enumerate(trajectories):
            nid, wgt = self.BuildTrajectoryRow(traj)
            rows.append(np.full(len(nid), it, dtype=int))
            cols.append(nid)
            vals.append(wgt)
            if verbose and (it + 1) % 500 == 0:
                print("  trajectory operator: %i / %i" % (it + 1, len(trajectories)), flush=True)
        A = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                              shape=(len(trajectories), self._n_nodes_total))
        return A

    def DropEmptyRows(self, A, trajectories: list):
        """Remove flamelets whose reconstruction carries no on-manifold weight.

        Such rows constrain nothing and make the dense Gram matrix of the projection singular.
        """
        keep = np.asarray(np.abs(A).sum(axis=1)).ravel() > 0
        return A[keep], [t for t, k in zip(trajectories, keep) if k], keep

    def TrajectoryTargets(self, trajectories: list, var_name: str):
        """Target integral of each trajectory: the flamelet's OWN integral over the portion of
        the trajectory that lies inside the tabulated mixture fraction range.

        For flamelets fully inside the range this is zero (the closure source is a divergence
        of a flux vanishing at both flamelet ends). For flamelets that only partially overlap
        the table, the reconstruction should reproduce the redistribution the flamelet performs
        over the covered portion, which is the corresponding partial integral -- demanding zero
        there would inject an error equal to the flux difference across the covered interval.
        """
        b = np.zeros(len(trajectories))
        for it, t in enumerate(trajectories):
            if var_name not in t["stored"]:
                continue
            z = t["Z"]
            on = (z >= self._z_levels[0]) & (z <= self._z_levels[-1])
            b[it] = np.sum(np.where(on, t["w_quad"], 0.0) * t["stored"][var_name])
        return b

    def Audit(self, A, theta: np.ndarray, trajectories: list, var_name: str,
              targets: np.ndarray = None):
        """Normalized trajectory integral defects of the reconstructed source field."""
        if targets is None:
            targets = self.TrajectoryTargets(trajectories, var_name)
        scales = np.array([t["scales"][var_name] for t in trajectories])
        return (A.dot(theta) - targets) / scales

    def Project(self, A, theta0: np.ndarray, trajectories: list, var_name: str,
                targets: np.ndarray = None, rel_floor: float = 1e-3, tikhonov: float = 1e-10,
                verbose: bool = True, frozen: np.ndarray = None):
        """Magnitude-weighted minimum-norm correction enforcing A (theta0 + delta) = targets.

        Solves min ||W^-1 delta||_2 subject to the trajectory constraints, with
        W = diag(|theta0| + floor), so nodes carrying large source magnitudes absorb most of
        the correction and near-zero (bulk state) nodes stay near zero:
            delta = W^2 A^T (A W^2 A^T + lam I)^-1 (targets - A theta0).

        The Tikhonov weight only conditions the dense Gram solve. A 5-fold study over the H2
        manifold (lam/trace from 1e-10 to 3) showed the defect on UNSEEN flamelets increasing
        monotonically with lam, so the near-unregularized value is deliberate: the constraint
        rows carry genuine information rather than noise, and damping them only degrades the
        result. Note that the projection makes the reconstruction conservative along the
        flamelets it is fitted on; states far from any flamelet retain a reconstruction error
        set by the table's resolution of the (stiff) source profile.

        Nodes flagged in ``frozen`` (e.g. those a dead band has switched off) get zero weight, so
        the correction is carried entirely by the remaining nodes and leaves them untouched.
        """
        if targets is None:
            targets = self.TrajectoryTargets(trajectories, var_name)
        floor = rel_floor * max(np.max(np.abs(theta0)), 1e-32)
        w2 = (np.abs(theta0) + floor) ** 2
        if frozen is not None:
            w2 = np.where(frozen, 0.0, w2)
        AW = A.multiply(w2[np.newaxis, :]).tocsr()      # A W^2
        G = (AW @ A.T).toarray()                        # A W^2 A^T (n_traj x n_traj)
        lam = tikhonov * np.trace(G) / max(G.shape[0], 1)
        G[np.diag_indices_from(G)] += lam
        r0 = A.dot(theta0) - targets
        y = np.linalg.solve(G, r0)
        delta = -(AW.T).dot(y)
        theta = theta0 + delta
        if verbose:
            scales = np.array([t["scales"][var_name] for t in trajectories])
            before = np.abs(r0 / scales)
            after = np.abs((A.dot(theta) - targets) / scales)
            print("%s conservation projection: |defect|/scale median %.3e -> %.3e, max %.3e -> %.3e"
                  % (var_name, np.median(before), np.median(after), np.max(before), np.max(after)))
            print("%s nodal correction: max |delta| = %.3e (field max |S| = %.3e)"
                  % (var_name, np.max(np.abs(delta)), np.max(np.abs(theta0))))
        return theta

    def SplitToLevels(self, theta: np.ndarray):
        """Split a global nodal vector back into per-level arrays."""
        return [theta[self._node_offset[i]:self._node_offset[i + 1]]
                for i in range(len(self._level_nodes))]


def LoadTrajectories(flamelet_data_dir: str, pv_species: list, pv_weights: list,
                     mix_status_header: str = "phi_", verbose: bool = True):
    """Load all flamelet trajectories below a flamelet data directory."""
    files = _collect_flamelet_files(flamelet_data_dir, mix_status_header)
    trajectories = []
    for i, f in enumerate(files):
        t = _load_trajectory(f, pv_species, pv_weights)
        if t is not None:
            trajectories.append(t)
        if verbose and (i + 1) % 500 == 0:
            print("  loaded %i / %i flamelet trajectories" % (i + 1, len(files)), flush=True)
    if verbose:
        print("Loaded %i flamelet trajectories." % len(trajectories))
    return trajectories
