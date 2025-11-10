# -*- coding: utf-8 -*-
"""
viscous_disk_evolution.py

This script simulates the 1D viscous evolution of a thin accretion disk.
It solves the surface density evolution equation using a forward-time,
centered-space (FTCS) finite difference scheme.

Parameters are loaded from the `config.py` file.
This script is structured with a main() function for clarity.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from config import params  # Import the single parameters object

def main(p):
    """
    Runs the main simulation logic.
    
    Args:
        p: The SimulationParameters object (from config.py)
           containing all simulation parameters.
    """
    
    # --- All parameters are now accessed via `p` ---
    # Example: p.N_r, p.alpha, p.M_kg

    # --- Setup the Radial Grid ---
    r = np.linspace(p.r_min, p.r_max, p.N_r)
    dr = r[1] - r[0]

    # --- Define Viscosity Profile (alpha-disk model) ---
    Omega_k = np.sqrt(p.G * p.M_kg / r**3)
    T = p.T_const * (r / p.r_min)**(-0.75)
    c_s = np.sqrt(p.k_B * T / (p.mu * p.m_p))
    H = c_s / Omega_k
    nu = p.alpha * c_s * H

    # --- Set the Initial Surface Density Profile (Gaussian Ring) ---
    Sigma = p.Sigma0 * np.exp(-(r - p.r0)**2 / (2 * p.sigma0**2))
    total_mass_initial = np.sum(2 * np.pi * r * Sigma * dr)
    print(f"Initial disk mass: {total_mass_initial:.2e} kg")

    # --- Prepare for Plotting ---
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_times = np.linspace(0, p.N_t * p.dt, 5) # Timestamps to plot
    plot_indices = (plot_times / p.dt).astype(int)

    # --- Main Evolution Loop (FTCS Scheme) ---
    print("Starting simulation...")

    for t_step in range(p.N_t + 1):
        if t_step in plot_indices:
            time_in_years = (t_step * p.dt) / (365.25 * 24 * 3600)
            ax.plot(r / 1e11, Sigma, label=f't = {time_in_years:.1f} years')
            print(f"  ... plotting at step {t_step}, time = {time_in_years:.1f} years")

        # --- Physics: Calculate the Right-Hand-Side (RHS) of the equation ---
        F = nu * Sigma * np.sqrt(r)
        dF_dr = np.gradient(F, dr, edge_order=2)
        d_dr_term = r**(1/2) * dF_dr
        rhs = (3 / r) * np.gradient(d_dr_term, dr, edge_order=2)
        
        # --- Time-Step: Apply the Change ---
        Sigma_new = Sigma + p.dt * rhs
        
        # --- Boundary Conditions ---
        Sigma_new[0] = 0.0  # Inner edge (accretion)
        Sigma_new[-1] = 0.0 # Outer edge (zero density)
        
        # --- Update for next loop ---
        Sigma = Sigma_new

    print("Simulation finished.")

    # --- Final Mass Check ---
    total_mass_final = np.sum(2 * np.pi * r * Sigma * dr)
    print(f"Final disk mass: {total_mass_final:.2e} kg")
    print(f"Mass change: {((total_mass_final - total_mass_initial) / total_mass_initial) * 100:.2f}%")

    # --- Finalize and Save Plot ---
    ax.set_xlabel('Radius ($10^{11}$ m)')
    ax.set_ylabel('Surface Density $Sigma$ (kg/m$^2$)')
    ax.set_title('Viscous Evolution of an Accretion Disk')
    ax.legend()
    ax.grid(True, alpha=0.2)
    plt.tight_layout()

    # --- Save Plot to `results` folder ---
    script_dir = Path(__file__).parent
    results_dir = script_dir.parent / 'results'
    results_dir.mkdir(exist_ok=True)  # Create the folder if it doesn't exist
    save_path = results_dir / 'surface_density_evolution.png'

    plt.savefig(save_path)
    print(f"\nPlot saved to {save_path}")


# --- This is the standard Python entry point ---
# It checks if the script is being run directly.
if __name__ == "__main__":
    # If so, call the main() function and pass it our imported params.
    main(params)