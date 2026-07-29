#!/usr/bin/env python3
"""Exact companion for THM-2857 endpoint Galois carry torsor."""

from __future__ import annotations

from hashlib import sha256
from math import comb, gcd
from pathlib import Path


N = 1183
P = 13
RELATIVE_STEP = 91
ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    ROOT / "04-computation/lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    ROOT / "05-knowledge/results/lrc14_q3_q11_transverse_endpoint_horn_thm2847.out":
        "155fce129c750a9505fdda3c71a250ff3a57fcd4044bb1df941da83c08baee1d",
    ROOT / "04-computation/lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py":
        "2227f59c717095da0f2042096ada145de4e3661530c9aa2cc9020f42c8237a8b",
    ROOT / "05-knowledge/results/lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out":
        "424fd2e83a618f862a5ee1b5f073a93fe236d10cdc5412eab1b54dec5e537eac",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


def transform(values, omega, modulus):
    return tuple(
        sum(
            values[r] * pow(omega, (-frequency * r) % P, modulus)
            for r in range(P)
        )
        % modulus
        for frequency in range(P)
    )


def predicted_power_support(power):
    return tuple(sorted(3 * j % P for j in range(power + 1)))


def main():
    # Gal(Q(zeta_1183)/Q(zeta_91)) is the kernel of reduction mod 91.
    units = tuple((1 + RELATIVE_STEP * r) % N for r in range(P))
    require(
        len(set(units)) == P
        and all(gcd(unit, N) == 1 for unit in units)
        and {
            left * right % N
            for left in units
            for right in units
        }
        == set(units)
        and all((13 * (unit - 1)) % N == 0 for unit in units),
        "relative Galois group failed",
    )
    require(
        (6 * 156) // (6 * 12) == P,
        "relative cyclotomic degree changed",
    )

    # c=zeta^624-zeta^510 and sigma_r(zeta)=zeta^(1+91r).
    orbit_exponents = tuple(
        ((624 * unit) % N, (510 * unit) % N)
        for unit in units
    )
    formula_exponents = tuple(
        (624, (510 + RELATIVE_STEP * (3 * r % P)) % N)
        for r in range(P)
    )
    require(orbit_exponents == formula_exponents, "Galois orbit formula failed")
    require(
        len({right for _, right in orbit_exponents}) == P,
        "endpoint orbit is not free",
    )
    canonical_scalar_exponents = orbit_exponents[0]
    sigma_one_exponents = orbit_exponents[1]
    require(
        canonical_scalar_exponents == (624, 510)
        and sigma_one_exponents == (624, 783)
        and canonical_scalar_exponents != sigma_one_exponents,
        "canonical-linear versus semilinear exponent boundary failed",
    )

    # The centered orbit is exactly the faithful character 3.  Its forward
    # difference is the same character multiplied by omega^3-1.
    linear_support = (0, 3)
    centered_support = (3,)
    derivative_support = (3,)
    power_four_support = predicted_power_support(4)
    power_five_support = predicted_power_support(5)
    power_ten_support = predicted_power_support(10)
    require(
        all(
            len(predicted_power_support(power)) == power + 1
            and all(comb(power, j) != 0 for j in range(power + 1))
            for power in range(1, P)
        )
        and power_four_support == (0, 3, 6, 9, 12)
        and power_five_support == (0, 2, 3, 6, 9, 12)
        and power_ten_support == (0, 1, 2, 3, 4, 5, 6, 8, 9, 11, 12),
        "power Fourier support law failed",
    )

    # A single exact finite-field realization independently certifies all
    # orbit, power, transform, and difference claims.  Since the values are
    # algebraic integers, distinct reduction implies distinct characteristic-
    # zero values.
    prime = 4733
    zeta = 625
    require(
        pow(zeta, N, prime) == 1
        and pow(zeta, N // 7, prime) != 1
        and pow(zeta, N // 13, prime) != 1,
        "finite-field element lacks exact order 1183",
    )
    omega = pow(zeta, RELATIVE_STEP, prime)
    orbit = tuple(
        (
            pow(zeta, 624 * unit, prime)
            - pow(zeta, 510 * unit, prime)
        )
        % prime
        for unit in units
    )
    mean = pow(zeta, 624, prime)
    centered = tuple((value - mean) % prime for value in orbit)
    derivative = tuple(
        (orbit[(r + 1) % P] - orbit[r]) % prime
        for r in range(P)
    )
    power_four_orbit = tuple(pow(value, 4, prime) for value in orbit)
    power_five_orbit = tuple(pow(value, 5, prime) for value in orbit)
    power_ten_orbit = tuple(pow(value, 10, prime) for value in orbit)
    require(
        len(set(orbit)) == P
        and len(set(power_four_orbit)) == P
        and len(set(power_five_orbit)) == P
        and len(set(power_ten_orbit)) == P,
        "endpoint values fail to separate thirteen sections",
    )
    require(
        tuple(i for i, value in enumerate(transform(orbit, omega, prime)) if value)
        == linear_support
        and tuple(
            i
            for i, value in enumerate(transform(centered, omega, prime))
            if value
        )
        == centered_support
        and tuple(
            i
            for i, value in enumerate(transform(derivative, omega, prime))
            if value
        )
        == derivative_support
        and tuple(
            i
            for i, value in enumerate(transform(power_four_orbit, omega, prime))
            if value
        )
        == power_four_support
        and tuple(
            i
            for i, value in enumerate(transform(power_five_orbit, omega, prime))
            if value
        )
        == power_five_support
        and tuple(
            i
            for i, value in enumerate(transform(power_ten_orbit, omega, prime))
            if value
        )
        == power_ten_support,
        "finite-field Fourier certificate failed",
    )

    # Galois norm and the q7-to-q11 group-ring alignment.
    norm_orbit = 1
    for value in orbit:
        norm_orbit = norm_orbit * value % prime
    norm_formula = (
        pow(pow(zeta, 624, prime), P, prime)
        - pow(pow(zeta, 510, prime), P, prime)
    ) % prime
    require(norm_orbit == norm_formula != 0, "relative norm formula failed")
    require((3 * 10 + 7) % P == 11, "q7-to-q11 alignment failed")
    require(
        comb(8, 4) == 70 and comb(10, 5) == 252,
        "Gaussian lift coefficient changed",
    )

    horn_output = (
        ROOT / "05-knowledge/results/lrc14_q3_q11_transverse_endpoint_horn_thm2847.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    bockstein_output = (
        ROOT / "05-knowledge/results/lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    require(
        "zeta1183_exponents=(624, 510)" in horn_output
        and "group_ring_over_inherited_K=c*z^3" in horn_output
        and "inverse_monomial_exponent=10" in horn_output,
        "pinned endpoint scalar or group-ring edge changed",
    )
    require(
        "identity=L9L8=T*L4; carry=1" in bockstein_output
        and "semantic_leg=THM2835_PINNED:"
        "449_QA(q11,a)_to_QAB(q7,a+1)" in bockstein_output,
        "pinned carry triangle or semantic leg changed",
    )

    print("ENDPOINT GALOIS CARRY TORSOR AND PHASE ALIGNMENT SIDECAR")
    print(
        "relative_field=Q(zeta1183)/Q(zeta91); degree=13; "
        "Gal={sigma_r:zeta->zeta^(1+91r)}"
    )
    print(
        "endpoint_orbit=c_r=A-B*omega^(3r); orbit_size=13; "
        "mean=A; centered_character=3; derivative_character=3"
    )
    print(
        "endpoint_alignment=(c_r*x^3)^10*(449*x^7)="
        "c_r^10*(449*x^11); tenth_power_orbit_size=13"
    )
    print(
        f"tenth_power_Fourier_support={power_ten_support}; "
        "missing=(7,10); point_separating=1"
    )
    print(
        "factorial_exit_tensor=first_exit_m<13_gives_m+1_carry_channels; "
        f"m4_support={power_four_support}; m4_point_separating=1"
    )
    print(
        f"m5_support={power_five_support}; "
        "Gaussian_exit4=70*c_r^4*L(H^4); "
        "Gaussian_exit5=252*c_r^5*L(H^5)"
    )
    print(
        f"exact_specialization=(prime={prime},zeta1183={zeta}); "
        "relative_norm=A^13-B^13_nonzero"
    )
    print(
        "canonical_current_T=K0_linear; scalar_exponents_stay=(624,510); "
        "sigma1_exponents=(624,783); mismatch=1"
    )
    print(
        "missing_bridge=semilinear_clutch_T_on_scalar_equals_sigma1; "
        "future_test=one_new_carried_endpoint_coefficient_then_T2380_pair_twist"
    )
    print(
        "scope=coefficient_Galois_torsor_only; "
        "no_physical_ancestry_action=1; no_response_section=1; "
        "no_positivity_or_LRC14_conclusion=1"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
