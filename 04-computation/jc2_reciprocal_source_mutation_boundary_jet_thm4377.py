#!/usr/bin/env python3
"""Exact audit for THM-4377.

The program deliberately reconstructs every source fibre twice: first by
scanning the source exponent ``e`` and testing all four exponents, and then
from the closed interval formulas.  It checks the reciprocal exponent
mutation, its three coordinate-wall defects, the matched-core clock, and the
eventual boundary-jet cokernel.
"""

from __future__ import annotations

from collections import Counter


ELL_MAX = 90
W_MAX = 180


class Audit:
    def __init__(self) -> None:
        self.assertions = 0

    def check(self, condition: bool, message: object) -> None:
        self.assertions += 1
        if not condition:
            raise RuntimeError(f"audit failure: {message}")


AUDIT = Audit()


def source_tuple(ell: int, u: int, v: int, e: int) -> tuple[int, int, int, int]:
    s = (ell + 1) // 2
    return (2 * v + e - ell, s + u - v - 1 - e, v - e, e)


def direct_fibre(ell: int, u: int, v: int) -> tuple[tuple[int, int, int, int], ...]:
    # The equation c+e=v forces 0<=e<=v for a polynomial source monomial.
    return tuple(
        source_tuple(ell, u, v, e)
        for e in range(v + 1)
        if min(source_tuple(ell, u, v, e)) >= 0
    )


def interval_fibre(ell: int, u: int, v: int) -> tuple[tuple[int, int, int, int], ...]:
    s = (ell + 1) // 2
    lo = max(0, ell - 2 * v)
    hi = min(v, s + u - v - 1)
    return tuple(source_tuple(ell, u, v, e) for e in range(lo, hi + 1))


def mutate(monomial: tuple[int, int, int, int], d: int) -> tuple[int, int, int, int]:
    a, b, c, e = monomial
    return (a + 2 * d, b - 2 * d, c + d, e)


def r_word(exponents: set[int]) -> str:
    if not exponents:
        return "0"
    words = []
    for e in sorted(exponents, reverse=True):
        words.append("1" if e == 0 else ("r" if e == 1 else f"r^{e}"))
    return "+".join(words)


def signed_r_word(positive: set[int], negative: set[int]) -> str:
    pieces = []
    for sign, exponents in (("+", positive), ("-", negative)):
        for e in sorted(exponents, reverse=True):
            atom = "1" if e == 0 else ("r" if e == 1 else f"r^{e}")
            if not pieces and sign == "+":
                pieces.append(atom)
            else:
                pieces.append(sign + atom)
    return "".join(pieces) if pieces else "0"


def analyse_case(ell: int, w: int, d: int) -> dict[str, object]:
    s = (ell + 1) // 2
    plus = direct_fibre(ell, w + d, w)
    minus = direct_fibre(ell, w, w + d)
    plus_closed = interval_fibre(ell, w + d, w)
    minus_closed = interval_fibre(ell, w, w + d)
    AUDIT.check(plus == plus_closed, ("plus interval", ell, w, d))
    AUDIT.check(minus == minus_closed, ("minus interval", ell, w, d))
    AUDIT.check(bool(plus) and bool(minus), ("bilateral nonempty", ell, w, d))

    plus_set = set(plus)
    minus_set = set(minus)
    matched_plus = {m for m in plus if m[1] >= 2 * d}
    b_defect = plus_set - matched_plus
    matched_minus = {mutate(m, d) for m in matched_plus}
    AUDIT.check(matched_minus <= minus_set, ("forward mutation", ell, w, d))

    inverse_matched = {
        m for m in minus if m[0] >= 2 * d and m[2] >= d
    }
    AUDIT.check(matched_minus == inverse_matched, ("inverse mutation", ell, w, d))
    for m in matched_plus:
        AUDIT.check(mutate(mutate(m, d), -d) == m, ("inverse", ell, w, d, m))

    a_defect = {m for m in minus if m[0] < 2 * d}
    c_defect = {m for m in minus if m[2] < d}
    AUDIT.check(not (a_defect & c_defect), ("disjoint target walls", ell, w, d))
    AUDIT.check(
        minus_set == matched_minus | a_defect | c_defect,
        ("target decomposition", ell, w, d),
    )
    AUDIT.check(
        len(minus_set) == len(matched_minus) + len(a_defect) + len(c_defect),
        ("target count", ell, w, d),
    )

    a_plus = max(0, ell - 2 * w)
    b_plus = min(w, s + d - 1)
    a_minus = max(0, ell - 2 * w - 2 * d)
    b_minus = min(w + d, s - d - 1)
    nu = max(0, min(w, s - d - 1) - a_plus + 1)
    beta = max(0, b_plus - max(a_plus, s - d) + 1)
    alpha = max(0, min(b_minus, a_plus - 1) - a_minus + 1)
    gamma = max(0, b_minus - max(a_minus, w + 1) + 1)
    AUDIT.check(len(matched_plus) == nu, ("nu", ell, w, d))
    AUDIT.check(len(b_defect) == beta, ("beta", ell, w, d))
    AUDIT.check(len(a_defect) == alpha, ("alpha", ell, w, d))
    AUDIT.check(len(c_defect) == gamma, ("gamma", ell, w, d))
    AUDIT.check(len(plus) == nu + beta, ("plus count", ell, w, d))
    AUDIT.check(len(minus) == nu + alpha + gamma, ("minus count", ell, w, d))
    AUDIT.check(nu <= s - d, ("matched ceiling", ell, w, d))
    AUDIT.check(not (b_defect and c_defect), ("opposite high walls", ell, w, d))

    if d == 0:
        AUDIT.check(plus == minus, ("fixed fibre", ell, w))
        AUDIT.check(nu == len(plus), ("fixed matching", ell, w))
    else:
        AUDIT.check((nu == s - d) == (w >= s), ("sharp match clock", ell, w, d))
        if w < s:
            AUDIT.check(bool(a_defect), ("preclock a wall", ell, w, d))
        else:
            expected_beta = min(w - s + d + 1, 2 * d)
            AUDIT.check(not a_defect and not c_defect, ("postclock surjection", ell, w, d))
            AUDIT.check(len(b_defect) == expected_beta, ("postclock jet", ell, w, d))
            AUDIT.check(
                {m[3] for m in minus} <= {m[3] for m in plus},
                ("postclock e inclusion", ell, w, d),
            )
            if w >= s + d - 1:
                AUDIT.check(len(b_defect) == 2 * d, ("full jet rank", ell, w, d))
                AUDIT.check(
                    {m[1] for m in b_defect} == set(range(2 * d)),
                    ("full u jet", ell, w, d),
                )

    plus_e = {m[3] for m in plus}
    minus_e = {m[3] for m in minus}
    matched_e = {m[3] for m in matched_plus}
    AUDIT.check(
        plus_e - minus_e == {m[3] for m in b_defect},
        ("graded plus defect", ell, w, d),
    )
    AUDIT.check(
        minus_e - plus_e == {m[3] for m in a_defect | c_defect},
        ("graded minus defect", ell, w, d),
    )
    AUDIT.check(matched_e == plus_e & minus_e, ("graded match", ell, w, d))

    return {
        "plus": plus,
        "minus": minus,
        "matched": matched_plus,
        "b_defect": b_defect,
        "a_defect": a_defect,
        "c_defect": c_defect,
        "symmetric_defect": len(plus_e ^ minus_e),
    }


def e_values(monomials: object) -> tuple[int, ...]:
    return tuple(sorted(m[3] for m in monomials))  # type: ignore[union-attr]


def main() -> None:
    cases = 0
    patterns: Counter[str] = Counter()
    max_symmetric_defect = 0
    max_controls: list[tuple[int, int, int]] = []

    for ell in range(2, ELL_MAX + 1):
        s = (ell + 1) // 2
        rho = (ell + 2) // 3
        for d in range(s):
            for w in range(rho, W_MAX + 1):
                data = analyse_case(ell, w, d)
                key = "".join(
                    name
                    for name, field in (("a", "a_defect"), ("b", "b_defect"), ("c", "c_defect"))
                    if data[field]
                ) or "none"
                patterns[key] += 1
                defect = int(data["symmetric_defect"])
                if defect > max_symmetric_defect:
                    max_symmetric_defect = defect
                    max_controls = [(ell, w, d)]
                elif defect == max_symmetric_defect:
                    max_controls.append((ell, w, d))
                AUDIT.check(defect <= 2 * s - 2, ("sharp defect ceiling", ell, w, d))
                cases += 1

    no_match = analyse_case(3, 1, 1)
    transient = analyse_case(10, 4, 1)
    two_target_walls = analyse_case(15, 5, 1)
    stable = analyse_case(10, 8, 4)

    print("THM-4377 reciprocal source mutation and boundary-jet audit")
    print(f"universe: ell=2..{ELL_MAX}, rho<=w<={W_MAX}, 0<=d<s")
    print(f"bilateral cases checked: {cases}")
    print(f"explicit assertions: {AUDIT.assertions}")
    print("observed wall-pattern counts: " + ", ".join(f"{k}={patterns[k]}" for k in sorted(patterns)))
    print(f"largest symmetric e-fibre defect: {max_symmetric_defect}")
    print(f"first maximum control (ell,w,d): {max_controls[0]}")
    print()
    print("no-match boundary control ell=3, (u,v)=(2,1)<->(1,2)")
    print(f"  plus e={e_values(no_match['plus'])}, minus e={e_values(no_match['minus'])}")
    print(f"  matched e={e_values(no_match['matched'])}")
    print()
    print("THM-4375 first hostile ell=10, (u,v)=(5,4)<->(4,5)")
    print(f"  plus e={e_values(transient['plus'])}, minus e={e_values(transient['minus'])}")
    print(f"  matched e={e_values(transient['matched'])}")
    print(f"  b-wall e={e_values(transient['b_defect'])}")
    print(f"  a-wall e={e_values(transient['a_defect'])}, c-wall e={e_values(transient['c_defect'])}")
    print(
        "  graded difference H_+-H_-="
        + signed_r_word(
            {m[3] for m in transient["b_defect"]},
            {m[3] for m in transient["a_defect"] | transient["c_defect"]},
        )
    )
    print()
    print("two-target-wall control ell=15, (u,v)=(6,5)<->(5,6)")
    print(f"  plus e={e_values(two_target_walls['plus'])}, minus e={e_values(two_target_walls['minus'])}")
    print(f"  a-wall e={e_values(two_target_walls['a_defect'])}, c-wall e={e_values(two_target_walls['c_defect'])}")
    print()
    print("full boundary-jet control ell=10, (u,v)=(12,8)<->(8,12)")
    print(f"  plus e={e_values(stable['plus'])}, minus e={e_values(stable['minus'])}")
    print(f"  matched e={e_values(stable['matched'])}")
    print(f"  cokernel b-degrees={tuple(sorted(m[1] for m in stable['b_defect']))}")
    print(f"  stable graded defect={r_word({m[3] for m in stable['b_defect']})}")
    print()
    print("PASS: mutation, three-wall decomposition, sharp clocks, and graded cokernel agree")


if __name__ == "__main__":
    main()
