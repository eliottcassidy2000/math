#!/usr/bin/env python3
"""Dependency-free algebra audit for THM-4268.

This path deliberately imports neither SymPy nor the primary companion.  It
checks the load-bearing polynomial identities by evaluation in a deterministic
finite box and by their independently entered coefficient formulae.  Degree
bounds make the box test exact over Q for each univariate slice used below.
"""

from __future__ import annotations

from fractions import Fraction


def face(U: int, W: int, Z: int, S: int, P: int) -> int:
    return 1 - U * P**6 - W * S**2 * P**5 - Z * S**4 * P**4


def audit_box() -> int:
    checks = 0
    values = (-3, -2, -1, 1, 2, 3)
    s_values = tuple(range(-7, 8))
    for U in values:
        for W in values:
            for Z in values:
                D = W * W - 4 * U * Z
                Lam = U + W + Z
                for S in s_values:
                    for P in values:
                        sbr = W * P + 2 * Z * S * S
                        pbr = 6 * U * P * P + 5 * W * S * S * P + 4 * Z * S**4
                        lhs = 2 * Z * pbr - (4 * Z * S * S + 3 * W * P) * sbr
                        assert lhs == -3 * D * P * P
                        assert face(U, W, Z, S, S * S) == 1 - Lam * S**12
                        checks += 2
    return checks


def audit_kummer() -> int:
    checks = 0
    for T in range(-9, 10):
        for R in range(-9, 10):
            if T == 0 and R == 0:
                continue
            Uh = 4 * T * T * R * R
            Zh = (T * T - R * R) ** 2
            assert Uh + Zh == (T * T + R * R) ** 2
            if T and R:
                q = Fraction(Zh, Uh)
                q_neg = Fraction(((-T) ** 2 - R * R) ** 2, 4 * (-T) ** 2 * R * R)
                q_inv = Fraction((R * R - T * T) ** 2, 4 * R * R * T * T)
                assert q == q_neg == q_inv
            checks += 1
    return checks


def main() -> None:
    box_checks = audit_box()
    kummer_checks = audit_kummer()
    # Direct edge-discriminant constants, independently entered from
    # disc(a*x^n+b)=(-1)^(n(n-1)/2)n^n*a^(n-1)*b^(n-1).
    assert -(4**4) == -256
    assert 6**6 == 46656
    print("theorem=THM-4268-independent")
    print(f"critical_and_contact_box_checks={box_checks}")
    print(f"kummer_projective_checks={kummer_checks}")
    print("critical_identity=pass")
    print("contact_substitution=pass")
    print("kummer_V4_invariance=pass")
    print("hostile=Spec(R[X]/(W*X-1))_special_empty_generic_nonempty")
    print("scope=algebra_audit_only_properness_is_theorem_proof")


if __name__ == "__main__":
    main()
