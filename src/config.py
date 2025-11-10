# -*- coding: utf-8 -*-
"""
config.py

This file defines all simulation parameters using a dataclass.
This provides type hinting, auto-completion, and a central place
for all configurable values.
"""

from dataclasses import dataclass

# --- Physical Constants ---
# These are defined outside the class as true constants
G_CONST = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
M_SUN_KG = 1.989e30   # Mass of the Sun (kg)
M_P_KG = 1.67e-27     # Mass of a proton (kg)
K_B_JK = 1.38e-23     # Boltzmann constant (J/K)


@dataclass
class SimulationParameters:
    """Holds all parameters for the accretion disk simulation."""

    # --- Simulation Grid ---
    N_r: int = 200          # Number of radial grid points
    r_min: float = 0.1e11     # Inner radius of the grid (m)
    r_max: float = 10.0e11    # Outer radius of the grid (m)
    N_t: int = 15000        # Number of time steps
    dt: float = 1.0e8         # Time step size (s)

    # --- Disk & Star ---
    M_star_solar: float = 1.0   # Mass of the central object (in solar masses)
    alpha: float = 0.01       # Shakura-Sunyaev alpha parameter
    T_const: float = 1000.0   # Proportionality constant for T ~ r^(-3/4) (K)
    mu: float = 2.3           # Mean molecular weight (for H2)

    # --- Initial Conditions ---
    r0: float = 2.0e11        # Initial radius of the gas ring (m)
    sigma0: float = 0.5e11    # Width of the initial gas ring (m)
    Sigma0: float = 1e5       # Peak initial surface density (kg/m^2)

    # --- Read-only Physical Constants (for easy access) ---
    G: float = G_CONST
    m_p: float = M_P_KG
    k_B: float = K_B_JK

    @property
    def M_kg(self) -> float:
        """
        Calculates the mass of the central star in kilograms.
        This is a "computed property".
        """
        return self.M_star_solar * M_SUN_KG

# --- Global Instance ---
# Create a single, importable instance of the parameters.
# In other files, you will just: `from config import params`
params = SimulationParameters()