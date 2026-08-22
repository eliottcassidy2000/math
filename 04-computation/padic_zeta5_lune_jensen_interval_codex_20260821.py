"""Validated Jensen-energy hostile for the published zeta_2(5) lune.

Status: FINITE-INTERVAL certificate plus an exact modular/Jensen reduction.

Uses python-flint complex balls.  A cell is discarded (given lower bound zero)
whenever any branch operation becomes indeterminate or the enclosure for |Z|
reaches one.  The two retained integrals count explicit Gamma_0(2) companion
preimages of the published Hauptmodul.  Their lower sum therefore gives a
rigorous lower bound for the Jensen collision energy, conditional only on the
exact algebraic/modular reduction documented in the session reflection.

This script excludes 13.806686 for the *unchanged published lune/arithmetic
packet*.  It does not exclude a different conformal domain, local system, or
arithmetic denominator improvement.
"""

from __future__ import annotations

from flint import acb, arb, ctx


ctx.dps = 50
PI = arb.pi()
I = acb(0, 1)
A = arb(82) / 841
THETA0 = (A / 2).acos()
CANDIDATE_REQUIRED_ENERGY = arb("2.2100946773550234")


def invpsi(q: acb) -> acb:
    """Exact inverse z=q(63q-58)/(2(29q-14)) of the published lune map."""

    return q * (63 * q - 58) / (2 * (29 * q - 14))


def integrand_cell(mid: arb, radius: arb, sign: int) -> arb:
    """Lower enclosure of max(0,-log|z_sign(t)|) on one t-ball."""

    t = arb(mid, radius)
    z = (I * t).exp()
    root = (1 - z * (-I * THETA0).exp()).sqrt() * (
        1 - z * (I * THETA0).exp()
    ).sqrt()
    q = (arb(29) / 63) * (1 + z - root)
    tau = q.log() / (2 * PI * I)
    transformed_q = (2 * PI * I * tau / (1 + sign * 2 * tau)).exp()
    pulled_back = invpsi(transformed_q)
    modulus = abs(pulled_back)
    upper = modulus.upper()
    if not upper < 1:
        return arb(0)

    # The inverse equation is
    #   63q^2-(58+58z)q+28z=0,
    # so its two q-roots multiply to 4z/9.  The branch normalized by
    # psi(0)=0 is the smaller one.  Discard the cell unless that sheet is
    # certified throughout the ball.
    q_upper = abs(transformed_q).upper()
    z_lower = modulus.lower()
    if not q_upper * q_upper < arb(4) * z_lower / 9:
        return arb(0)
    value = -arb(upper).log()
    lower = value.lower()
    return arb(lower) if lower > 0 else arb(0)


def lower_sum(sign: int, cells: int, left: str, right: str) -> arb:
    aa = arb(left)
    bb = arb(right)
    width = (bb - aa) / cells
    radius = width / 2
    total = arb(0)
    discarded = 0
    for j in range(cells):
        mid = aa + (arb(2 * j + 1) * width) / 2
        try:
            value = integrand_cell(mid, radius, sign)
        except (ValueError, ZeroDivisionError):
            value = arb(0)
        if value == 0:
            discarded += 1
        total += width * value
    normalized = total / (2 * PI)
    print(
        f"sign={sign:+d} interval=[{left},{right}] cells={cells} "
        f"discarded={discarded} raw_lower={total} normalized_lower={normalized}"
    )
    return normalized


def main() -> None:
    # Each branch's mass is concentrated on one side of pi.  We deliberately
    # omit the branch cut and every other interval; the integrand is
    # nonnegative, so discarded territory only weakens the lower bound.
    minus = lower_sum(-1, 12000, "0.5", "3.1415")
    plus = lower_sum(+1, 12000, "3.1417", "5.8")
    modular_lower = minus + plus
    energy_lower = (arb(14) / 29).log() + modular_lower
    print("sum modular lower", modular_lower)
    print("energy lower", energy_lower)
    print("candidate-required energy", CANDIDATE_REQUIRED_ENERGY)
    print("certified margin", energy_lower - CANDIDATE_REQUIRED_ENERGY)
    if not energy_lower > CANDIDATE_REQUIRED_ENERGY:
        raise RuntimeError("candidate energy is not excluded")
    print("published-lune candidate exclusion PASS")


if __name__ == "__main__":
    main()
