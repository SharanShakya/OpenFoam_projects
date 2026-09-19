"""
contour_generator.py
Rao 80% bell nozzle contour generation, CSV export, and plotting.
"""

import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def generate_rao_contour(Rt: float, Re: float, epsilon: float, num_points: int = 200):
    """
    Generate a Rao 80% bell nozzle contour (parabolic approximation).

    Returns
    -------
    x_coords, r_coords : np.ndarray
        Axial and radial wall coordinates [m]
    L_nozzle : float
        Nozzle divergent length [m]
    """
    # Convergent section
    R_conv = 1.5 * Rt
    x_conv = np.linspace(-1.5 * Rt, 0, int(num_points * 0.3))  # matches R_conv exactly
    r_conv = Rt + (R_conv - np.sqrt(np.maximum(0, R_conv**2 - x_conv**2)))

    # Divergent section (80% bell parabola)
    L_nozzle = 0.8 * ((math.sqrt(epsilon) - 1) * Rt / math.tan(math.radians(15)))
    theta_i = math.radians(20.0)
    theta_e = math.radians(7.0)

    t = np.linspace(0, 1, int(num_points * 0.7))
    P0 = np.array([0.0, Rt])
    P2 = np.array([L_nozzle, Re])

    m1 = math.tan(theta_i)
    m2 = math.tan(theta_e)
    C_x = (Re - Rt - m2 * L_nozzle) / (m1 - m2)
    C_y = Rt + m1 * C_x
    P1 = np.array([C_x, C_y])

    bell_pts = (1 - t)[:, None]**2 * P0 + 2 * (1 - t)[:, None] * t[:, None] * P1 + t[:, None]**2 * P2

    x_coords = np.concatenate([x_conv, bell_pts[:, 0]])
    r_coords = np.concatenate([r_conv, bell_pts[:, 1]])

    return x_coords, r_coords, L_nozzle


def export_contour_csv(x_wall, r_wall, output_folder: str, filename: str = "nozzle_contour.csv") -> str:
    """Save contour coordinates to CSV. Returns full save path."""
    os.makedirs(output_folder, exist_ok=True)
    save_path = os.path.join(output_folder, filename)
    contour_df = pd.DataFrame({"x_m": x_wall, "r_m": r_wall})
    contour_df.to_csv(save_path, index=False)
    return save_path


def plot_contour(x_wall, r_wall, Rt, Re, Dt, De, Ln,
                  output_folder: str, filename: str = "nozzle_geometry.png",
                  title: str = "Rocket Nozzle Geometry (80% Bell)",
                  show: bool = True) -> str:
    """Plot and save the nozzle contour. Returns full save path."""
    os.makedirs(output_folder, exist_ok=True)
    save_path = os.path.join(output_folder, filename)

    plt.figure(figsize=(11, 5))

    x_mm = x_wall * 1000
    r_mm = r_wall * 1000

    plt.plot(x_mm, r_mm, 'b-', linewidth=2.5, label='Nozzle Contour')
    plt.plot(x_mm, -r_mm, 'b-', linewidth=2.5)
    plt.fill_between(x_mm, r_mm, -r_mm, color='lightskyblue', alpha=0.15)

    plt.axhline(0, color='gray', linestyle='--', linewidth=1, label='Centerline')
    plt.axvline(0, color='red', linestyle=':', linewidth=1.2, label='Throat Section')

    plt.plot(0, Rt * 1000, 'ro', markersize=6, label=f'Throat ($D_t$ = {Dt*1000:.1f} mm)')
    plt.plot(Ln * 1000, Re * 1000, 'go', markersize=6, label=f'Exit ($D_e$ = {De*1000:.1f} mm)')

    plt.title(title, fontsize=13, fontweight='bold')
    plt.xlabel("Axial Position, $x$ [mm]", fontsize=11)
    plt.ylabel("Radial Position, $r$ [mm]", fontsize=11)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc='upper left', frameon=True)
    plt.tight_layout()

    plt.savefig(save_path, dpi=300)
    if show:
        plt.show()

    return save_path
