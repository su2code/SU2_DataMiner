"""
Plot free flame solutions from flamelet_data/freeflame_data.

For every equivalence-ratio (or mixture-fraction) folder found, one figure per
variable in PLOT_VARIABLES is created.  Each figure overlays all flame solutions
found in that folder, plotted as a function of the progress variable (PV).

PV definition (from 0_generate_config.py):
    PV = sum(w_i * Y_i)  with  species=['H2O','CO2','CO','H2'], weights=[1,1,1,1]

Figures are saved to  flamelet_data/freeflame_plots/<folder>/<variable>.png
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Progress variable definition (must match 0_generate_config.py)
# ---------------------------------------------------------------------------
PV_SPECIES  = ['H2O', 'CO2', 'CO', 'H2']
PV_WEIGHTS  = [1.0,   1.0,   1.0,  1.0]

# ---------------------------------------------------------------------------
# Variables to plot (column names as they appear in the CSV header)
# ---------------------------------------------------------------------------
PLOT_VARIABLES = [
    'Temperature',
    'Heat_Release',
    'Density',
    'Cp',
    'Conductivity',
    'ViscosityDyn',
    'EnthalpyTot',
    'MixtureFraction',
    'Velocity',
]


def compute_pv(variables: list[str], data: np.ndarray) -> np.ndarray:
    """Compute progress variable as weighted sum of species mass fractions."""
    pv = np.zeros(data.shape[0])
    for sp, w in zip(PV_SPECIES, PV_WEIGHTS):
        col = 'Y-' + sp
        if col in variables:
            pv += w * data[:, variables.index(col)]
        else:
            print(f"  Warning: species column '{col}' not found – skipping.")
    return pv


def plot_folder(folder_path: str, folder_name: str, output_dir: str):
    """Load all CSV files in one mixture-fraction/phi folder and plot."""
    csv_files = sorted([
        f for f in os.listdir(folder_path)
        if f.endswith('.csv')
    ])
    if not csv_files:
        return

    # ---- collect data from all files in this folder -------------------------
    flames: list[dict] = []
    variables: list[str] | None = None

    for csv_file in csv_files:
        fpath = os.path.join(folder_path, csv_file)
        with open(fpath, 'r') as fid:
            header = fid.readline().strip()
        cols = header.split(',')
        data = np.loadtxt(fpath, delimiter=',', skiprows=1, ndmin=2)
        if data.shape[0] == 0:
            continue
        if variables is None:
            variables = cols
        pv = compute_pv(cols, data)
        # sort by PV for cleaner lines
        order = np.argsort(pv)
        flames.append({'name': csv_file, 'data': data[order], 'pv': pv[order], 'vars': cols})

    if not flames or variables is None:
        return

    os.makedirs(output_dir, exist_ok=True)

    # ---- one figure per variable --------------------------------------------
    for var in PLOT_VARIABLES:
        # Check at least one flame has this column
        avail = [f for f in flames if var in f['vars']]
        if not avail:
            print(f"  Variable '{var}' not found in any file of {folder_name}, skipping.")
            continue

        fig, ax = plt.subplots(figsize=(7, 4.5))
        for f in avail:
            idx = f['vars'].index(var)
            label = f['name'].replace('freeflamelet_', '').replace('.csv', '')
            ax.plot(f['pv'], f['data'][:, idx], lw=1.2, label=label)

        ax.set_xlabel('Progress Variable  [ - ]')
        ax.set_ylabel(var)
        ax.set_title(f'{folder_name}  –  {var}')
        ax.legend(fontsize=6, ncol=2, loc='best')
        ax.grid(True, ls='--', alpha=0.4)
        fig.tight_layout()

        out_file = os.path.join(output_dir, f'{var}.png')
        fig.savefig(out_file, dpi=150)
        plt.close(fig)
        print(f'  Saved: {out_file}')


def main():
    script_dir  = os.path.dirname(os.path.abspath(__file__))
    freeflame_root = os.path.join(script_dir, 'flamelet_data', 'freeflame_data')
    plots_root     = os.path.join(script_dir, 'flamelet_data', 'freeflame_plots')

    if not os.path.isdir(freeflame_root):
        print(f"Error: directory not found:\n  {freeflame_root}")
        sys.exit(1)

    folders = sorted([
        d for d in os.listdir(freeflame_root)
        if os.path.isdir(os.path.join(freeflame_root, d))
    ])

    if not folders:
        print("No equivalence-ratio folders found in freeflame_data/.")
        sys.exit(0)

    print(f"Found {len(folders)} folder(s): {folders}")
    for folder in folders:
        folder_path = os.path.join(freeflame_root, folder)
        output_dir  = os.path.join(plots_root, folder)
        print(f"\nProcessing folder: {folder}")
        plot_folder(folder_path, folder, output_dir)

    print("\nDone.")


if __name__ == '__main__':
    main()
