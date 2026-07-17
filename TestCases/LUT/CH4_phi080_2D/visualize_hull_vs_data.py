#!/usr/bin/env python3
"""
Diagnostic script to visualize actual flamelet data points vs the constructed mesh hull.
This helps identify whether the convex hull is cutting off concave data regions.

Usage:
    python visualize_hull_vs_data.py <path_to_config_file> <mixture_fraction_value>

Example:
    python visualize_hull_vs_data.py TestCases/LUT/Flamelet_Table_CH4_phi080/TableGeneration.cfg 0.055
"""

import sys
from Manifold_Generation.LUT.FlameletTableGeneration import SU2TableGenerator
from Common.DataDrivenConfig import Config_FGM

def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print("\nUsing default values:")
        config_path = "TableGeneration.cfg"
        z_value = 0.055
        print(f"  Config: {config_path}")
        print(f"  Mixture fraction: {z_value}")
    else:
        config_path = sys.argv[1]
        z_value = float(sys.argv[2])
    
    print("\n" + "="*70)
    print("Hull vs Data Visualization")
    print("="*70)
    print(f"Loading configuration: {config_path}")
    
    # Load configuration
    Config = Config_FGM(config_path)
    
    # Initialize table generator
    print("Initializing SU2TableGenerator...")
    Tgen = SU2TableGenerator(Config, n_near=14, p_fac=1)
    
    # Configure table settings (match your 3_generate_table.py settings)
    Tgen.SetEquivalenceRatioLimits(phi_min=0.80, phi_max=0.80)
    Tgen.SetNTableLevels(1)
    Tgen.SetRefinementFields(["ProdRateTot_PV", "Y_dot_net-CO", "Y_dot_pos-CO","Y_dot_neg-CO"])
    
    # Set mesh resolution
    Tgen.SetBaseCellSize(5e-3)
    Tgen.SetRefinedCellSize(5e-3)
    Tgen.SetRefinementRadius(5e-3)
    Tgen.SetRefinementMethod("gradient")
    Tgen.SetMaxRefinementSeeds(500)
    Tgen.SetHullCellSize(5.0e-3)
    
    # Set table axes
    Tgen.SetTableAxes(level_cv_name="MixtureFraction",
                      plane_cv_names=["ProgressVariable", "EnthalpyTot"])
    
    print(f"\nGenerating table level visualization for Z = {z_value}...")
    print("This will create two plots:")
    print("  1. Actual data points with convex hull overlay")
    print("  2. Mesh nodes with actual data overlay")
    print("\nLook for data points that fall OUTSIDE the red hull boundary.")
    print("These are being excluded from the table!\n")
    
    # Visualize with hull diagnostic enabled
    Tgen.VisualizeTableLevel(z_value, visualize_hull=True, show=True)
    
    print("\n" + "="*70)
    print("Visualization complete!")
    print("="*70)

if __name__ == "__main__":
    main()
