from dataclasses import dataclass, replace
from typing import Any, Callable, List, Sequence
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# 设置中文字体和负号显示（如有需要）
matplotlib.rcParams['font.sans-serif'] = ['Heiti TC', 'STHeiti', 'SimHei', 'Arial Unicode MS']
matplotlib.rcParams['axes.unicode_minus'] = False

# === Physical Constants (SI Units) ===
EPSILON_0 = 8.854187817e-12  # Vacuum permittivity (F/m)
M_E = 9.10938356e-31  # Electron mass (kg)
Q_E = 1.602176634e-19  # Elementary charge (C)
C = 299792458  # Speed of light (m/s)


@dataclass(frozen=True)
class BeamParameters:
    """Container for beam and plasma properties used by the envelope solver."""

    gamma: float
    current: float  # Beam current (A)
    normalized_emittance: float  # Normalized emittance (m-rad)
    initial_angle: float  # Initial divergence angle (rad)
    ion_density: float  # Background ion density (m^-3)

    def __post_init__(self) -> None:
        if self.gamma < 1:
            raise ValueError("gamma must be >= 1")
        if self.current <= 0:
            raise ValueError("Beam current must be positive")
        if self.normalized_emittance <= 0:
            raise ValueError("Normalized emittance must be positive")
        if self.initial_angle <= 0:
            raise ValueError("Initial angle must be positive")
        if self.ion_density < 0:
            raise ValueError("Ion density must be non-negative")

    @property
    def beta(self) -> float:
        """Relativistic beta factor (v/c)."""
        return np.sqrt(1.0 - 1.0 / self.gamma**2)

    @property
    def alfven_current(self) -> float:
        """Alfven current (A)."""
        return 4.0 * np.pi * EPSILON_0 * M_E * C**3 * self.beta * self.gamma / Q_E

    @property
    def geometric_emittance(self) -> float:
        """Geometric emittance (m-rad)."""
        return self.normalized_emittance / (self.beta * self.gamma)

    def neutralization_factor(self, radius: float) -> float:
        """Charge neutralization factor f_e."""
        if radius <= 0:
            raise ValueError("Radius must be positive when computing neutralization factor")
        return Q_E * self.beta * C * np.pi * radius**2 * self.ion_density / self.current

    def x_factor(self, radius: float) -> float:
        """Dimensionless x factor for the beam-plasma system."""
        fe = self.neutralization_factor(radius)
        return (fe * self.gamma**2 - 1.0) / (self.gamma**2 - 1.0)

    def initial_radius(self) -> float:
        """Initial beam radius derived from emittance and divergence angle."""
        return self.geometric_emittance / self.initial_angle


def make_envelope_rhs(params: BeamParameters) -> Callable[[float, Sequence[float]], List[float]]:
    """
    Create the first-order system that represents the beam envelope equation.

    The state vector is [R(z), R'(z)].
    """
    inverse_alfven_current = 1.0 / params.alfven_current
    geometric_emittance = params.geometric_emittance

    def rhs(_: float, state: Sequence[float]) -> List[float]:
        radius, slope = state
        # Avoid numerical issues if the solver probes non-physical radii.
        radius_for_force = max(float(radius), 1.0e-12)
        xfac = params.x_factor(radius_for_force)
        second_derivative = (
            -2.0 * xfac * params.current * inverse_alfven_current / radius_for_force
            + (geometric_emittance**2) / radius_for_force**3
        )
        return [float(slope), second_derivative]

    return rhs


def sweep_final_radius_by_gamma(
    params: BeamParameters,
    gamma_values: Sequence[float],
    z_max: float,
) -> np.ndarray:
    """Compute the final beam radius for each gamma value in gamma_values."""
    final_radii: List[float] = []
    for gamma_value in gamma_values:
        gamma_params = replace(params, gamma=float(gamma_value))
        solution = solve_envelope(gamma_params, z_max)
        final_radii.append(solution.y[0][-1])
    return np.asarray(final_radii)


def solve_envelope(
    params: BeamParameters,
    z_max: float,
    num_points: int = 2000,
) -> Any:
    """Integrate the beam envelope equation over z in [0, z_max]."""
    if z_max <= 0:
        raise ValueError("z_max must be positive")
    if num_points < 2:
        raise ValueError("num_points must be at least 2")

    initial_conditions = [params.initial_radius(), params.initial_angle]
    z_span = (0.0, float(z_max))
    z_eval = np.linspace(z_span[0], z_span[1], num_points)

    solution = solve_ivp(
        fun=make_envelope_rhs(params),
        t_span=z_span,
        y0=initial_conditions,
        t_eval=z_eval,
    )

    if not solution.success:
        raise RuntimeError(f"Envelope integration failed: {solution.message}")

    return solution


def main() -> None:
    """Solve and plot the beam envelope evolution for a sample parameter set."""
    params = BeamParameters(
        # 504
        # gamma=40.0,
        # current=50e-6,
        # normalized_emittance=1e-6,
        # initial_angle=5e-6,
        # ion_density=1e5,
        # QY
        # gamma=70.0,
        # current=1.8,
        # normalized_emittance=3e-5,
        # initial_angle=2e-4,
        # ion_density=1e5,
        gamma=200.0,
        current=1.8,
        normalized_emittance=3e-6,
        initial_angle=2e-6,
        ion_density=1e5,
    )
    z_max = 1e5

    solution = solve_envelope(params, z_max)
    final_radius = solution.y[0][-1]

    para_res_str = "\n".join([
    fr"能量={params.gamma/2:.1f} $\mathrm{{MeV}}$",
    # fr"流强={params.current*1e6:.1f} $\mathrm{{\mu A}}$",
    fr"流强={params.current:.1f} $\mathrm{{A}}$",
    fr"归一化发射度={params.normalized_emittance*1e6:.1f} $\mathrm{{mm-mrad}}$",
    fr"初始发散角={params.initial_angle*1e6:.1f} $\mathrm{{\mu rad}}$",
    # fr"初始发散角={params.initial_angle*1e3:.1f} $\mathrm{{mrad}}$",
    fr"初始半径={params.initial_radius()*1e2:.2f} $\mathrm{{cm}}$",
    fr"背景离子密度={params.ion_density*1e-6:.1f} $\mathrm{{cm}}^{{-3}}$",
    fr"到靶半径={final_radius:.2f} $\mathrm{{m}}$",
])



    plt.figure("电子束包络半径随传播距离变化")
    plt.plot(solution.t*1e-3, solution.y[0])
    # plt.xticks(np.arange(0,z_max,10), np.linspace(0,int(z_max/1000),10))
    plt.xlabel("传输距离(km)")
    plt.ylabel("包络半径(m)")
    plt.title("电子束包络半径随传输距离变化")
    plt.text(x=0.5,y=0.60*final_radius,
             s=para_res_str)
    plt.grid(False)

    # E_values = np.linspace(10.0, 100.0, 100)
    # # gamma_values = np.linspace(20.0, 200.0, 100)
    # final_radii = sweep_final_radius_by_gamma(params, 2*E_values, z_max)

    # plt.figure("到靶半径与能量关系")
    # plt.plot(E_values, final_radii)
    # plt.xlabel("能量 (MeV)")
    # plt.ylabel("到靶半径(m)")
    # plt.title("到靶半径与能量关系")
    # plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()
