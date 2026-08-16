#!/usr/bin/env python3
"""Clean-room two-root audit of cyclic weighted-family z-primitivity.

This companion does not import Sympy or the existing remainder probe.  For
every inverse degree n >= 3 it specializes

    T_n(w) = w^n - w^(n-1) + P*w - Q

at P=2^(-(n-1)), Q=0.  The visible roots w=0 and w=1/2 give different exact
values of the z-reconstruction numerator H.  Consequently H mod T_n is not
constant over Q(P,Q), so z is not in the target field.  The script verifies
the uniform formulas, their strict inequality, and the full weighted-map
forward replay for n=3..256 using only fractions.Fraction.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from pathlib import Path


EXPECTED_SEMANTIC_SHA256: str | None = "2c8d9f191ae5a8dddd6ef15feeb20d6564491bd0be6ea34c8e662428f74b729b"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def digest_fraction(value: Fraction) -> str:
    return sha256(f"{value.numerator}/{value.denominator}".encode("ascii")).hexdigest()


def inverse_polynomial(n: int, w: Fraction, P: Fraction, Q: Fraction) -> Fraction:
    return w**n - w ** (n - 1) + P * w - Q


def seed(n: int, w: Fraction) -> Fraction:
    """p_ell(w), where ell=n-2."""

    return (n - 1) * w ** (n - 2) - n * w ** (n - 1)


def q_seed(n: int, w: Fraction) -> Fraction:
    """q_ell(w), where ell=n-2."""

    return (n - 2) * w ** (n - 1) - (n - 1) * w**n


def a_parameter(n: int) -> Fraction:
    return -Fraction(2 * n - 3, 2 * n - 4)


def gamma_value(n: int, w: Fraction, P: Fraction) -> Fraction:
    return P - seed(n, w)


def h_value(n: int, w: Fraction, P: Fraction) -> Fraction:
    gamma = gamma_value(n, w, P)
    a = a_parameter(n)
    return gamma * (gamma * (gamma - 1 + a) - a * w)


def reconstruct_source(n: int, w: Fraction, P: Fraction) -> tuple[Fraction, Fraction, Fraction]:
    gamma = gamma_value(n, w, P)
    require(gamma != 0, (n, w, "gamma-open reconstruction"))
    x = 1 / gamma
    y = w - gamma
    z = h_value(n, w, P)
    return x, y, z


def forward_weighted_map(
    n: int, source: tuple[Fraction, Fraction, Fraction]
) -> tuple[Fraction, Fraction, Fraction]:
    x, y, z = source
    a = a_parameter(n)
    v = x * y
    t = x * x * z
    u = 1 + v
    gamma = 1 + a * v + t
    require(gamma != 0, (n, source, "forward gamma"))
    w = u * gamma
    beta = 1 + seed(n, w) / gamma
    alpha = u + q_seed(n, w) / (gamma * gamma)
    return alpha / (x * x), beta / x, x * gamma


def audit_row(n: int) -> tuple[object, ...]:
    require(n >= 3, n)
    P = Fraction(1, 2 ** (n - 1))
    Q = Fraction(0)
    w0 = Fraction(0)
    w1 = Fraction(1, 2)
    m = n - 3
    delta = Fraction(1, 2 * n - 4)

    require(inverse_polynomial(n, w0, P, Q) == 0, (n, "root zero"))
    require(inverse_polynomial(n, w1, P, Q) == 0, (n, "root half"))

    gamma0 = gamma_value(n, w0, P)
    gamma1 = gamma_value(n, w1, P)
    require(gamma0 == P, (n, "gamma zero-root"))
    require(gamma1 == -m * P, (n, "gamma half-root"))

    H0 = h_value(n, w0, P)
    H1 = h_value(n, w1, P)
    expected_H0 = P * P * (P - 2 - delta)
    expected_H1 = -m * P * (
        m * P * (m * P + 2 + delta) + Fraction(1 + delta, 2)
    )
    require(H0 == expected_H0, (n, "H zero-root"))
    require(H1 == expected_H1, (n, "H half-root"))
    require(H0 < 0, (n, "H0 sign"))

    if n == 3:
        # The second branch has gamma=0 at this specialization, but it is a
        # lawful polynomial-identity witness: H(0)<0=H(1/2).
        require(gamma1 == 0 and H1 == 0 and H0 != H1, (n, "cubic boundary"))
        replay = "gamma_boundary_polynomial_witness"
    else:
        require(gamma1 != 0 and H1 < 0, (n, "open pair"))
        require(abs(H1) > P / 2, (n, "lower magnitude bound"))
        require(abs(H0) <= Fraction(9, 32) * P, (n, "upper magnitude bound"))
        require(abs(H1) > abs(H0), (n, "strict separation"))

        source0 = reconstruct_source(n, w0, P)
        source1 = reconstruct_source(n, w1, P)
        require(source0 != source1, (n, "distinct source points"))
        require(source0[2] != source1[2], (n, "distinct z values"))
        target = (Fraction(0), P, Fraction(1))
        require(forward_weighted_map(n, source0) == target, (n, "forward zero-root"))
        require(forward_weighted_map(n, source1) == target, (n, "forward half-root"))
        replay = "two_finite_sources_same_target_distinct_z"

    payload = (n, P, gamma0, gamma1, H0, H1, replay)
    return (
        n,
        str(P),
        str(gamma0),
        str(gamma1),
        digest_fraction(H0),
        digest_fraction(H1),
        replay,
        sha256(repr(payload).encode("ascii")).hexdigest(),
    )


def main() -> None:
    rows = tuple(audit_row(n) for n in range(3, 257))
    sample_indices = (0, 1, 2, 4, 8, 28)
    samples = tuple(rows[index] for index in sample_indices)
    ledger_sha256 = sha256(repr(rows).encode("ascii")).hexdigest()
    proof = (
        "P=2^(-(n-1));Q=0 gives roots 0 and 1/2 for every n>=3",
        "gamma_0=P;gamma_half=-(n-3)P",
        "H_0=P^2(P-2-delta)",
        "H_half=-(n-3)P*((n-3)P*((n-3)P+2+delta)+(1+delta)/2)",
        "n=3:H_0<0=H_half",
        "n>=4:abs(H_half)>P/2>=(16/9)abs(H_0)",
        "therefore_H_mod_T_is_not_constant_over_Q(P,Q)",
        "therefore_z_notin_target_field;S_n_point_stabilizer_maximality_gives_z_primitive",
    )
    scope = (
        "THM3448_explicit_cyclic_weighted_family_only",
        "independent_of_symbolic_remainder_recurrence",
        "no_arbitrary_weighted_seed_or_map_classification",
        "no_JC2_LRC_ancestry_current_or_H1_claim",
    )
    semantic_surface = (proof, scope, samples, ledger_sha256)
    semantic_sha256 = sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            (semantic_sha256, EXPECTED_SEMANTIC_SHA256),
        )

    source = Path(__file__).resolve()
    print("Cyclic weighted Keller family: independent two-root z-primitivity audit")
    print("status=FINITE_EXACT_PLUS_UNIFORM_RATIONAL_PROOF;scope=THM3448_cyclic_family")
    print(f"proof={proof}")
    print(f"sample_rows={samples}")
    print(f"exact_replay_range=n_3_through_256;ledger_sha256={ledger_sha256}")
    print("hostile=n_3_half_root_has_gamma_zero_so_only_polynomial_nonconstancy_is_used")
    print("positive=n_4_through_256_two_finite_sources_map_to_same_target_with_distinct_z")
    print(f"scope={scope}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_assert_truth_gates;fractions_only;no_randomness;no_elapsed_fields")
    print("commands=python -B 04-computation/jc_weighted_cyclic_z_primitivity_two_root_independent_audit_20260816.py;python -B -O 04-computation/jc_weighted_cyclic_z_primitivity_two_root_independent_audit_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
