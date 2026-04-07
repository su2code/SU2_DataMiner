"""
Compute a single burner-stabilized flame with settings matching SU2_DataMiner
and save the result to a CSV file.

Usage:
    python run_single_burner_flame.py
    # or override defaults on the command line:
    python run_single_burner_flame.py --phi 0.60 --T 300 --mdot 0.060
    # with a restart file as initial guess:
    python run_single_burner_flame.py --phi 0.60 --T 300 --mdot 0.045 --restart-file burner_phi0.60_T300_mdot0.060.csv
"""

import argparse
import csv
import numpy as np
import cantera as ct

# ── Settings matching 0_generate_config.py / DataGenerator_FGM defaults ──────
MECHANISM   = "gri30.yaml"
FUEL        = "CH4"
OXIDIZER    = "O2:0.21,N2:0.79"
TRANSPORT   = "mixture-averaged"    # same as Config.SetTransportModel
PRESSURE    = ct.one_atm

# Grid settings matching DataGenerator_FGM.__initial_grid* defaults
GRID_LENGTH = 0.1      # metres (DefaultSettings_FGM.initial_grid_length)
GRID_NP     = 20        # points  (DefaultSettings_FGM.initial_grid_Np)

# Refine criteria matching compute_SingleBurnerFlame
RATIO  = 3
#SLOPE  = 0.02
#CURVE  = 0.02
#PRUNE  = 0.01
SLOPE  = 0.05
CURVE  = 0.05
PRUNE  = 0.01

# PV species and weights from 0_generate_config.py
PV_SPECIES  = ["H2O", "CO2", "CO", "H2"]
PV_WEIGHTS  = [1.0,   1.0,   1.0,  1.0]


def compute_pv(gas: ct.Solution, Y: np.ndarray) -> float:
    """Sum weighted species mass fractions to get the progress variable."""
    pv = 0.0
    for sp, w in zip(PV_SPECIES, PV_WEIGHTS):
        idx = gas.species_index(sp)
        pv += w * float(Y[idx])
    return pv


def _load_csv_restart(flame: ct.BurnerFlame, csv_path: str):
    """Set T, velocity, and species profiles on *flame* from a SU2_DataMiner CSV file.

    The CSV rows are ordered from burner (unburnt) to products (burnt), matching
    the order Cantera expects. A uniform normalised position [0, 1] is used because
    the CSV does not store grid coordinates.
    """
    with open(csv_path, "r") as f:
        header = f.readline().strip().split(",")
    D = np.loadtxt(csv_path, delimiter=",", skiprows=1, ndmin=2)
    N = len(D)
    locs = np.linspace(0.0, 1.0, N)

    flame.set_profile("T", locs, D[:, header.index("Temperature")])

    if "Velocity" in header:
        flame.set_profile("velocity", locs, D[:, header.index("Velocity")])

    for sp in flame.gas.species_names:
        col = f"Y-{sp}"
        if col in header:
            flame.set_profile(sp, locs, D[:, header.index(col)])

    print(f"  Loaded {N}-point restart: T_max={D[:, header.index('Temperature')].max():.1f} K")


def run_burner_flame(phi: float, T_burner: float, mdot: float,
                     restart_file: str = None) -> ct.BurnerFlame:
    gas = ct.Solution(MECHANISM)
    gas.TP = T_burner, PRESSURE
    gas.set_equivalence_ratio(phi, FUEL, OXIDIZER)

    grid = np.linspace(0, GRID_LENGTH, GRID_NP)
    flame = ct.BurnerFlame(gas, grid=grid)
    flame.transport_model = TRANSPORT
    flame.set_refine_criteria(ratio=RATIO, slope=SLOPE, curve=CURVE, prune=PRUNE)
    flame.burner.mdot = mdot

    if restart_file is not None:
        print(f"  Loading initial guess from: {restart_file}")
        if restart_file.endswith(".csv"):
            _load_csv_restart(flame, restart_file)
        else:
            # YAML or HDF5 saved by a previous Cantera flame.save() call
            restart_data = ct.SolutionArray(gas)
            if restart_file.endswith(".h5") or restart_file.endswith(".hdf"):
                restart_data.read_hdf(restart_file)
            else:
                restart_data.restore(restart_file)
            flame.set_initial_guess(data=restart_data)

    print(f"Solving burner flame: phi={phi:.4f}  T={T_burner:.1f} K  mdot={mdot:.4e} kg/m2/s")
    flame.solve(loglevel=1, refine_grid=True, auto=False)
    print(f"  T_max = {flame.T.max():.1f} K   grid points = {len(flame.grid)}")
    return flame


def save_flame(flame: ct.BurnerFlame, filepath: str):
    """Save flame data in the exact column order used by SU2_DataMiner.__SaveFlameletData."""
    gas = flame.gas
    N   = len(flame.grid)
    MW  = gas.molecular_weights        # kg/kmol, shape (n_species,)

    # ── per-grid-point arrays ────────────────────────────────────────────────
    net_rates  = flame.net_production_rates          # (n_species, N) kmol/m3/s
    dest_rates = flame.destruction_rates
    pos_rates  = net_rates - dest_rates

    Y_dot_net = (net_rates  * MW[:, np.newaxis]).T    # (N, n_sp) kg/m3/s
    Y_dot_pos = (pos_rates  * MW[:, np.newaxis]).T
    Y_dot_neg = (dest_rates * MW[:, np.newaxis]).T / (flame.Y.T + 1e-11)

    cp_i   = (flame.partial_molar_cp.T   / MW)       # (N, n_sp)
    enth_i = (flame.partial_molar_enthalpies.T / MW)  # (N, n_sp)

    Le_i = (flame.thermal_conductivity / flame.cp_mass / flame.density_mass
            / (flame.mix_diff_coeffs + 1e-15))        # (N, n_sp)

    # mixture fraction via Bilger
    try:
        mf = flame.mixture_fraction("Bilger")         # (N,)
    except Exception:
        mf = np.array([gas.mixture_fraction(FUEL, OXIDIZER)
                       for _ in range(N)])

    # ── build header exactly as __SaveFlameletData ───────────────────────────
    sp = gas.species_names
    header  = ["Distance", "Velocity"]
    header += [f"Y-{s}"         for s in sp]
    header += [f"Y_dot_net-{s}" for s in sp]
    header += [f"Y_dot_pos-{s}" for s in sp]
    header += [f"Y_dot_neg-{s}" for s in sp]
    header += [f"Cp-{s}"        for s in sp]
    header += [f"h-{s}"         for s in sp]
    header += [f"Le-{s}"        for s in sp]
    header += ["EnthalpyTot", "MixtureFraction", "Temperature",
               "Density", "MolarWeightMix", "Cp", "Conductivity",
               "ViscosityDyn", "Heat_Release"]

    # ── build data rows ───────────────────────────────────────────────────────
    rows = []
    for i in range(N):
        gas.TPY = flame.T[i], PRESSURE, flame.Y[:, i]
        row  = [flame.grid[i], flame.velocity[i]]
        row += list(flame.Y[:, i])
        row += list(Y_dot_net[i])
        row += list(Y_dot_pos[i])
        row += list(Y_dot_neg[i])
        row += list(cp_i[i])
        row += list(enth_i[i])
        row += list(Le_i[i])
        row += [
            gas.enthalpy_mass,
            float(mf[i]),
            flame.T[i],
            gas.density,
            gas.mean_molecular_weight,
            gas.cp_mass,
            gas.thermal_conductivity,
            gas.viscosity,
            flame.heat_release_rate[i],
        ]
        rows.append(row)

    with open(filepath, "w", newline="") as fid:
        writer = csv.writer(fid)
        writer.writerow(header)
        writer.writerows(rows)

    print(f"Saved {N} points → {filepath}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phi",  type=float, default=0.60)
    parser.add_argument("--T",    type=float, default=300.0)
    parser.add_argument("--mdot", type=float, default=0.060)
    parser.add_argument("--restart-file", type=str, default=None,
                        help="CSV file written by SU2_DataMiner (or run_single_burner_flame.py) "
                             "used as the initial guess via set_profile(). "
                             "YAML / HDF5 files saved by flame.save() are also accepted.")
    args = parser.parse_args()

    flame = run_burner_flame(phi=args.phi, T_burner=args.T, mdot=args.mdot,
                             restart_file=args.restart_file)

    fname = f"burner_phi{args.phi:.2f}_T{int(args.T)}_mdot{args.mdot:.3f}.csv"
    save_flame(flame, fname)


if __name__ == "__main__":
    main()
