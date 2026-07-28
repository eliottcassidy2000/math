"""Exact finite/algebraic checks for THM-2685.

The codimension-one DVR and class-group arguments are mathematical proof,
not delegated to this companion.  This file checks the finite F2 code and
the polynomial identities used by the theorem.
"""

from itertools import permutations, product

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rank_f2(rows, width):
    work = [sum((bit & 1) << j for j, bit in enumerate(row)) for row in rows]
    rank = 0
    for j in range(width):
        pivot = next((i for i in range(rank, len(work)) if (work[i] >> j) & 1), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for i in range(len(work)):
            if i != rank and ((work[i] >> j) & 1):
                work[i] ^= work[rank]
        rank += 1
    return rank


def main():
    # Values of a functional (a,b) on the three nonzero classes
    # e1,e2,e1+e2 are exactly the even-weight [3,2,2] code.
    code = tuple(sorted((a, b, a ^ b) for a, b in product(range(2), repeat=2)))
    expected_code = ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0))
    require(code == expected_code, "quartic parity code changed")
    require(tuple(sum(word) for word in code) == (0, 2, 2, 2),
            "nonzero inertia weights changed")

    nonzero = ((1, 0), (0, 1), (1, 1))
    cycle = {nonzero[i]: nonzero[(i + 1) % 3] for i in range(3)}
    invariant_lines = tuple(v for v in nonzero if cycle[v] == v)
    require(not invariant_lines, "standard F2[C3] plane became reducible")
    require(len(tuple(permutations(nonzero))) == 6,
            "GL2(F2)=S3 permutation count changed")

    # THM-2681 divisor rows and maximal unramified quotient.
    boundary_rows = ((1, 0), (0, 1), (1, 1))
    boundary_rank = rank_f2(boundary_rows, 2)
    kernel = tuple(
        (x, y) for x, y in product(range(2), repeat=2)
        if all((a * x + b * y) % 2 == 0 for a, b in boundary_rows)
    )
    require(boundary_rank == 2 and kernel == ((0, 0),),
            "ordered-root boundary no longer kills every subcover")

    # Derive the squared-pair-sum resolvent from four roots with zero sum.
    x1, x2, x3 = sp.symbols("x1 x2 x3")
    x4 = -x1 - x2 - x3
    roots = (x1, x2, x3, x4)
    p = sp.expand(sum(roots[i] * roots[j] for i in range(4) for j in range(i + 1, 4)))
    e3 = sp.expand(sum(roots[i] * roots[j] * roots[k]
                       for i in range(4) for j in range(i + 1, 4)
                       for k in range(j + 1, 4)))
    q = -e3
    r = sp.expand(x1 * x2 * x3 * x4)
    u = tuple(sp.expand((x1 + roots[j]) ** 2) for j in (1, 2, 3))
    require(sp.expand(sum(u) + 2 * p) == 0,
            "quartic resolvent U^2 coefficient changed")
    require(sp.expand(sum(u[i] * u[j] for i in range(3) for j in range(i + 1, 3))
                      - (p**2 - 4 * r)) == 0,
            "quartic resolvent U coefficient changed")
    require(sp.expand(u[0] * u[1] * u[2] - q**2) == 0,
            "quartic resolvent constant term changed")

    # Integral-resolvent coefficient gate and its translation law.
    A, B, C, h, W = sp.symbols("A B C h W")
    cubic = W**3 + A * W**2 + B * W + C
    shifted = sp.Poly(sp.expand(cubic.subs(W, W + h)), W)
    Ah, Bh, Ch = shifted.all_coeffs()[1:]
    require(sp.expand(Ah * Bh - Ch + cubic.subs(W, -A - 2 * h)) == 0,
            "translated integral-resolvent gate changed")

    aa, bb, cc = sp.symbols("a b c")
    A0, B0, C0 = -1, bb * cc / 4, -aa * cc**2 / 4
    normalized_gate = sp.factor(A0 * B0 - C0)
    require(normalized_gate == cc * (aa * cc - bb) / 4,
            "THM-2681 literal-resolvent gate changed")

    # Ramified/even-base-change hostile: t has parity 110, s^2 has 000.
    ramified_word = (1, 1, 0)
    even_pullback_word = tuple((2 * bit) % 2 for bit in ramified_word)
    require(ramified_word in code and even_pullback_word == (0, 0, 0),
            "even base-change parity hostile changed")

    # Riemann-Hurwitz parity sanity checks.
    curve_genus = tuple((g, b, 2 * g - 1 + b // 2)
                        for g in range(3) for b in (0, 2, 4, 6)
                        if not (g == 0 and b == 0))
    require(curve_genus[0] == (0, 2, 0) and curve_genus[-1] == (2, 6, 6),
            "curve parity/genus table changed")

    print("THM-2685 EQUIVARIANT KUMMER BOUNDARY PARITY AUDIT")
    print(f"even_weight_code={code}; weights={tuple(sum(w) for w in code)}")
    print(f"C3_invariant_lines={invariant_lines}; GL2F2_permutations=6")
    print(f"thm2681_boundary_rows={boundary_rows}; rank={boundary_rank}; kernel={kernel}")
    print("quartic_squared_pair_resolvent=(1,2p,p^2-4r,-q^2): PASS")
    print(f"translated_gate=Ah*Bh-Ch=-R(-A-2h): {sp.factor(Ah * Bh - Ch)}")
    print(f"thm2681_literal_gate={normalized_gate}; nonsquare_odd_factors=(c,ac-b)")
    print(f"ramified_word={ramified_word}; even_base_change={even_pullback_word}")
    print(f"curve_genus_controls={curve_genus}")
    print("scope=finite code and polynomial identities only; DVR quasi-etale and class-group exactness are proved in text")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
