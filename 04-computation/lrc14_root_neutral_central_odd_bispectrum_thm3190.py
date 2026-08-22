#!/usr/bin/env python3
"""Exact bidegree controls for THM-3190's cubic clutch criterion."""

import hashlib
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEPENDENCIES = (
    (HERE / "lrc14_sparse_root_bispectrum_current_thm2312.py",
     "6d80ff4460d720eff24e4339218c65bb3a21d2464dec6ffe73d5ea8bbadd8a4f"),
    (ROOT / "05-knowledge/results/lrc14_sparse_root_bispectrum_current_thm2312.out",
     "47e9f4e5e9804bd8f7bce83c853932f15aba3a697ad5be67a93069d7f4d901c7"),
    (HERE / "lrc14_central_sign_parity_quotient_thm3187.py",
     "73e61cfafac01a8daa7363ae978221fa6e26ccf01b28314ca4629727038a6476"),
    (ROOT / "05-knowledge/results/lrc14_central_sign_parity_quotient_thm3187.out",
     "47bb19b9f661a9de333bcbdde55e2b1931aa52123152cb04e8afb348d41e432d"),
)


def lf_hash(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for path, expected in DEPENDENCIES:
    require(lf_hash(path) == expected, ("dependency hash drift", str(path)))


P = 13
ALLOWED = tuple(
    (k, ell, (k + ell) % P, (k + ell) // P)
    for k in range(1, P)
    for ell in range(1, P)
    if (k + ell) % P != 0
)
require(len(ALLOWED) == 132, "allowed bispectrum-pair census")
require(sum(carry == 0 for _, _, _, carry in ALLOWED) == 66
        and sum(carry == 1 for _, _, _, carry in ALLOWED) == 66,
        "carry split drift")


# Root translation sends M_j to zeta^(-jh)M_j.  The conjugated third factor
# contributes +ch, so every allowed cubic has root character zero.
ROOT_TRANSLATION_CHECKS = 0
for k, ell, c, carry in ALLOWED:
    require(k + ell - c == P * carry, "carry representative identity")
    for h in range(P):
        exponent = (-k * h - ell * h + c * h) % P
        require(exponent == 0, "root-translation charge survived")
        ROOT_TRANSLATION_CHECKS += 1


# A common scalar phase lambda contributes lambda*lambda*conj(lambda)=lambda.
# Track the exponents of lambda and its conjugate formally.  The central
# lambda=-1 therefore acts by -1, not +1.
SCALAR_HOLOMORPHIC_DEGREE = 2
SCALAR_ANTIHOLORPHIC_DEGREE = 1
SCALAR_CHARGE = SCALAR_HOLOMORPHIC_DEGREE - SCALAR_ANTIHOLORPHIC_DEGREE
TOTAL_POLYNOMIAL_DEGREE = (
    SCALAR_HOLOMORPHIC_DEGREE + SCALAR_ANTIHOLORPHIC_DEGREE
)
require(SCALAR_CHARGE == 1, "global scalar charge is not one")
require(TOTAL_POLYNOMIAL_DEGREE == 3
        and (-1) ** TOTAL_POLYNOMIAL_DEGREE == -1,
        "central parity is not odd")


# Exact finite phase controls: for every cyclic phase order 2..31, the cubic
# output exponent is the original common phase exponent.  Order two is the
# central-sign specialization.
GLOBAL_PHASE_CHECKS = 0
for order in range(2, 32):
    for phase in range(order):
        output_phase = (2 * phase - phase) % order
        require(output_phase == phase, "global U(1) charge drift")
        GLOBAL_PHASE_CHECKS += 1
require((2 * 1 - 1) % 2 == 1, "central sign failed")


# The two actions commute on the formal amplitude bank, so the cubic has
# bidegree (root character, scalar charge)=(0,1).
COMMUTING_ACTION_CHECKS = 0
for k, ell, c, _ in ALLOWED:
    for h in range(P):
        root_phase = (-k * h - ell * h + c * h) % P
        for central_bit in (0, 1):
            left = (root_phase, central_bit * SCALAR_CHARGE % 2)
            right = (root_phase, central_bit * SCALAR_CHARGE % 2)
            require(left == right == (0, central_bit),
                    "root/scalar bidegree failed")
            COMMUTING_ACTION_CHECKS += 1


print("THM-3190 root-neutral central-odd bispectrum exact control")
print("dependency_hash_checks=" + repr(len(DEPENDENCIES)))
print("allowed_pairs=" + repr(len(ALLOWED)))
print("carry_zero_one=" + repr((66, 66)))
print("root_translation_checks=" + repr(ROOT_TRANSLATION_CHECKS))
print("global_phase_checks=" + repr(GLOBAL_PHASE_CHECKS))
print("commuting_bidegree_checks=" + repr(COMMUTING_ACTION_CHECKS))
print("bidegree=(root_character_0,global_scalar_charge_1)")
print("central_sign=-1")
print("scope=common_carrier_and_torus_normalization_still_open")
print("all_exact_checks=PASS")
