"""Numerical scout for the archimedean energy in the zeta_2(5) holonomy bound.

This is deliberately a scout, not a rigorous certificate.  It evaluates the
specific lune uniformization printed in Calegari--Dimitrov--Tang's ICM survey
and compares the raw Haar log-energy with the rearrangement upper bound used
there.  Multiple resolutions and radial approaches are retained as hostile
controls against diagonal and boundary-branch quadrature artifacts.
"""

from __future__ import annotations

import math

import numpy as np


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lune_q(z: np.ndarray, a: float = 2 / 5) -> np.ndarray:
    """The normalized lune map psi_a, scaled into |q|<2/3.

    The published choice is a=2/5 and has psi'(0)=14/29.  The square root is
    the branch equal to 1 at z=0.  On |z|<1 the principal branch realizes it
    for this quadratic; radial limits are used below.
    """

    a2 = a * a
    radicand = 1 - 2 * (1 - 6 * a2 + a2 * a2) / (1 + a2) ** 2 * z + z * z
    return (2 / 3) * (1 + a2) / (2 * (1 - a2)) * (
        1 + z - np.sqrt(radicand)
    )


def lune_derivative(a: float) -> float:
    return (2 / 3) * (1 - a * a) / (1 + a * a)


def hauptmodul_x(q: np.ndarray, terms: int = 120) -> np.ndarray:
    """x(q)=q prod_{n>=1}(1+q^n)^24, truncated safely for |q|<=2/3."""

    ans = q.astype(np.complex128, copy=True)
    qn = q.astype(np.complex128, copy=True)
    for _ in range(1, terms + 1):
        ans *= (1 + qn) ** 24
        qn *= q
    return ans


def boundary_values(
    n: int, radius: float, a: float = 2 / 5
) -> tuple[np.ndarray, np.ndarray]:
    theta = 2 * math.pi * (np.arange(n) + 0.371) / n
    z = radius * np.exp(1j * theta)
    q = lune_q(z, a)
    return z, hauptmodul_x(q)


def off_diagonal_energy(values: np.ndarray, block: int = 256) -> float:
    """N^-2 sum_{i != j} log|v_i-v_j|, accumulated in blocks."""

    n = len(values)
    total = 0.0
    for lo in range(0, n, block):
        hi = min(n, lo + block)
        d = np.abs(values[lo:hi, None] - values[None, :])
        rows = np.arange(lo, hi)
        d[np.arange(hi - lo), rows] = 1.0
        total += float(np.log(d).sum())
    return total / (n * n)


def local_derivative_mean(z: np.ndarray, values: np.ndarray) -> float:
    """Centered cyclic estimate of mean log|d phi/dz| on the sample circle."""

    dz = np.roll(z, -1) - np.roll(z, 1)
    dv = np.roll(values, -1) - np.roll(values, 1)
    return float(np.log(np.abs(dv / dz)).mean())


def diagonal_corrected_energy(z: np.ndarray, values: np.ndarray) -> float:
    """Remove the universal N-point diagonal defect.

    The local spacing is |z|*|d(values)/dz|*2*pi/N.  Including log|z|
    matters for the radial hostile controls, even though it vanishes in the
    boundary limit.  For a linear circle values=c*z this correction is exact.
    """

    n = len(values)
    radius = float(np.abs(z).mean())
    derivative = local_derivative_mean(z, values)
    return off_diagonal_energy(values) - (
        math.log(n) - derivative - math.log(radius)
    ) / n


def corrected_from_parts(
    z: np.ndarray, raw: float, derivative_mean: float
) -> float:
    n = len(z)
    radius = float(np.abs(z).mean())
    return raw - (math.log(n) - derivative_mean - math.log(radius)) / n


def kappa_threshold(energy: float, derivative: float = 14 / 29) -> float:
    """Solve equality in the single-place specialization of Proposition 2.4."""

    gamma = 1 / 6
    alpha = 12 * math.log(2)
    base = math.log(derivative) + alpha - 175 / 36
    numerator = energy + alpha

    # 6 = numerator / (base-alpha*(2(1-gamma)k-(1-gamma)^2)/k^2)
    # gives a quadratic in k.  The larger positive root is the threshold.
    a = base - numerator / 6
    b = -2 * alpha * (1 - gamma)
    c = alpha * (1 - gamma) ** 2
    roots = np.roots([a, b, c])
    positive = [float(r.real) for r in roots if abs(r.imag) < 1e-10 and r.real > 0]
    return max(positive)


def energy_required(kappa: float, derivative: float = 14 / 29) -> float:
    """Energy which would make Proposition 2.4 an equality at kappa."""

    gamma = 1 / 6
    alpha = 12 * math.log(2)
    base = math.log(derivative) + alpha - 175 / 36
    correction = alpha * (
        2 * (1 - gamma) * kappa - (1 - gamma) ** 2
    ) / (kappa * kappa)
    return 6 * (base - correction) - alpha


def modular_preimage_lower(sample_count: int) -> tuple[float, float, float]:
    """Two explicit Gamma_0(2) collision contributions to Jensen energy.

    On the published lune boundary q=psi(z), write q=exp(2*pi*i*tau).
    The matrices [[1,0],[+/-2,1]] preserve the Hauptmodul

        x(q) = q product_n (1+q^n)^24.

    Pulling their q-images back through the exact rational inverse

        psi^{-1}(q) = q(63q-58)/(2(29q-14))

    produces two explicit interior preimages whenever their modulus is below
    one.  Jensen counts max(0,-log|z|) for each.  This routine evaluates the
    resulting one-dimensional integrals numerically; an interval enclosure is
    still needed to turn the numerical lower device into a theorem.
    """

    pi = math.pi
    angle0 = math.acos((82 / 841) / 2)
    t = 2 * pi * (np.arange(sample_count) + 0.371) / sample_count
    z = np.exp(1j * t)
    analytic_root = np.sqrt(1 - z * np.exp(-1j * angle0)) * np.sqrt(
        1 - z * np.exp(1j * angle0)
    )
    q = (29 / 63) * (1 + z - analytic_root)
    tau = np.log(q) / (2j * pi)

    contributions: list[float] = []
    for sign in (+1, -1):
        transformed_q = np.exp(2j * pi * tau / (1 + 2 * sign * tau))
        pulled_back = transformed_q * (63 * transformed_q - 58) / (
            2 * (29 * transformed_q - 14)
        )
        contributions.append(
            float(np.maximum(0, -np.log(np.abs(pulled_back))).mean())
        )
    lower = math.log(14 / 29) + sum(contributions)
    return contributions[0], contributions[1], lower


def main() -> None:
    published_rearrangement = 3.92881
    candidate = 13.806686
    print("published rearrangement energy", published_rearrangement)
    print("published threshold replay", f"{kappa_threshold(published_rearrangement):.9f}")
    required = energy_required(candidate)
    print(f"candidate energy for kappa={candidate}", f"{required:.12f}")
    print(
        "candidate collision excess I-log|phi'(0)|",
        f"{required - math.log(14 / 29):.12f}",
    )
    require(abs(required - 2.2100946773550234) < 1e-12, "candidate inversion drift")
    require(
        abs(kappa_threshold(published_rearrangement) - 19.74389750921136) < 1e-10,
        "published threshold replay drift",
    )
    print(
        "forbidden cross-template mix: circle energy 2.13322 + lune derivative ->",
        f"{kappa_threshold(2.13322, 14 / 29):.9f}",
    )
    print(
        "legitimate circle packet: circle energy 2.13322 + derivative 1/3 ->",
        f"{kappa_threshold(2.13322, 1 / 3):.9f}",
    )
    jensen_plus, jensen_minus, jensen_lower = modular_preimage_lower(2**18)
    require(abs(jensen_plus - jensen_minus) < 3e-5, "modular branches lost symmetry")
    require(jensen_lower > required + 0.28, "modular Jensen hostile lost its gap")
    print(
        "Gamma_0(2) modular-preimage contributions (VERIFIED-NUMERIC)",
        f"{jensen_plus:.10f}",
        f"{jensen_minus:.10f}",
    )
    print(
        "Jensen lower device",
        f"{jensen_lower:.10f}",
        "already exceeds candidate energy by",
        f"{jensen_lower-required:.10f}",
    )
    print("formal exclusion requires a validated interval enclosure of these integrals")
    print()

    # Exact discrete capacity control for a linear image of a circle.
    n_control = 1024
    r_control = 0.93
    c_control = 2.7
    angles = 2 * math.pi * (np.arange(n_control) + 0.371) / n_control
    z_control = r_control * np.exp(1j * angles)
    linear_control = diagonal_corrected_energy(z_control, c_control * z_control)
    require(
        abs(linear_control - math.log(c_control * r_control)) < 1e-12,
        "linear-circle diagonal correction failed",
    )
    print("linear-circle diagonal control PASS", f"{linear_control:.12f}")
    print()

    for radius in (0.99, 0.997, 0.999, 0.9997, 0.9999):
        for n in (512, 1024, 2048, 4096):
            z, values = boundary_values(n, radius)
            raw = off_diagonal_energy(values)
            deriv = local_derivative_mean(z, values)
            corrected = corrected_from_parts(z, raw, deriv)
            print(
                f"r={radius:.4f} n={n:4d} raw={raw:.10f} "
                f"dmean={deriv:.8f} corrected={corrected:.10f} "
                f"kappa={kappa_threshold(corrected):.8f}"
            )
        print()

    print("parameter scan (FINITE-NUMERIC, not a conformal-domain proof)")
    for a in np.linspace(0.0, 0.85, 18):
        n = 2048
        z, values = boundary_values(n, 0.9997, float(a))
        energy = diagonal_corrected_energy(z, values)
        derivative = lune_derivative(float(a))
        print(
            f"a={a:.3f} psi'={derivative:.9f} energy={energy:.9f} "
            f"collision_excess={energy-math.log(derivative):.9f} "
            f"kappa={kappa_threshold(energy, derivative):.8f}"
        )


if __name__ == "__main__":
    main()
