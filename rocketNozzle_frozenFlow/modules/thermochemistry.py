"""
thermochemistry.py
Combustion chemistry / equilibrium properties for LOX/CH4 propellant.
"""

import math

try:
    import cantera as ct
    HAS_CANTERA = True
except ImportError:
    HAS_CANTERA = False

Ru = 8314.462618  # Universal gas constant [J/(kmol K)]


def get_thermochemistry(Pc_pa: float, mixture_ratio: float):
    """
    Compute equilibrium combustion properties for LOX/CH4.

    Parameters
    ----------
    Pc_pa : float
        Chamber pressure [Pa]
    mixture_ratio : float
        Oxidizer/Fuel mass ratio (O/F)

    Returns
    -------
    Tc : float
        Chamber (stagnation) temperature [K]
    gamma : float
        Ratio of specific heats, cp/cv
    MW : float
        Mean molecular weight [kg/kmol]
    chem_source : str
        Description of the method used
    """
    if HAS_CANTERA:
        gas = ct.Solution('gri30.yaml')
        gas.TPX = 300, Pc_pa, f'CH4:1, O2:{mixture_ratio * 16.04 / 31.998}'
        gas.equilibrate('HP')

        Tc = gas.T
        gamma = gas.cp / gas.cv
        MW = gas.mean_molecular_weight
        return Tc, gamma, MW, "Cantera (GRI-30 Equilibrium)"
    else:
        return 3600.0, 1.20, 22.5, "Preliminary Fallback Values"


def specific_gas_constant(MW: float) -> float:
    """Specific gas constant R = Ru / MW [J/(kg K)]"""
    return Ru / MW