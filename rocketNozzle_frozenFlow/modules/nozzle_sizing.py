"""
nozzle_sizing.py
1-D isentropic gas dynamics and thrust chamber sizing.
"""

import math
from scipy.optimize import brentq

from thermochemistry import get_thermochemistry, specific_gas_constant


def area_ratio_from_mach(M: float, g: float) -> float:
    """Area ratio A/A* as a function of Mach number (isentropic flow)."""
    term = (2 / (g + 1)) * (1 + 0.5 * (g - 1) * M**2)
    return (1.0 / M) * (term ** ((g + 1) / (2 * (g - 1))))


def static_to_total_pressure(M: float, g: float) -> float:
    """Static-to-total pressure ratio P/P0 as a function of Mach number."""
    return (1.0 + 0.5 * (g - 1) * M**2) ** (-g / (g - 1))


def compute_nozzle_sizing(F_target: float, Pc: float, Pa: float, OF: float,
                           g0: float = 9.80665) -> dict:
    """
    Perform full 1-D thermochemical + gas dynamic sizing of a rocket nozzle.

    Parameters
    ----------
    F_target : float
        Design thrust [N]
    Pc : float
        Chamber pressure [Pa]
    Pa : float
        Ambient pressure [Pa] (nozzle designed for full expansion, Pe = Pa)
    OF : float
        Oxidizer/Fuel mass ratio
    g0 : float
        Standard gravity [m/s^2]

    Returns
    -------
    dict
        All key thermo/geometric/performance results.
    """
    Tc, gamma, MW, chem_source = get_thermochemistry(Pc, OF)
    R = specific_gas_constant(MW)

    # Characteristic velocity c*
    c_star = math.sqrt(R * Tc / gamma) * ((gamma + 1) / 2) ** ((gamma + 1) / (2 * (gamma - 1)))

    # Solve exit Mach number for Pe = Pa
    pe_pc_target = Pa / Pc
    f_mach = lambda M: static_to_total_pressure(M, gamma) - pe_pc_target
    Me = brentq(f_mach, 1.001, 10.0)

    epsilon = area_ratio_from_mach(Me, gamma)
    Te = Tc / (1.0 + 0.5 * (gamma - 1) * Me**2)
    Ve = Me * math.sqrt(gamma * R * Te)

    # Thrust coefficient & sizing
    Cf = Ve / c_star
    At = F_target / (Pc * Cf)
    Ae = epsilon * At

    mdot = Pc * At / c_star
    mdot_fuel = mdot / (1 + OF)
    mdot_ox = mdot - mdot_fuel

    Rt = math.sqrt(At / math.pi)
    Re = math.sqrt(Ae / math.pi)
    Dt, De = 2 * Rt, 2 * Re
    Isp = F_target / (mdot * g0)

    return {
        "chem_source": chem_source,
        "Tc": Tc, "gamma": gamma, "MW": MW, "R": R,
        "c_star": c_star, "Me": Me, "epsilon": epsilon,
        "Te": Te, "Ve": Ve, "Cf": Cf,
        "At": At, "Ae": Ae,
        "mdot": mdot, "mdot_fuel": mdot_fuel, "mdot_ox": mdot_ox,
        "Rt": Rt, "Re": Re, "Dt": Dt, "De": De,
        "Isp": Isp,
    }