#!/usr/bin/env python3
"""Exact referee for the local-ruler deletion-surrogate obstruction.

The row

    W = {14,28,29,56,70,84,98,112,126,153,154,168}
    V = W union {182}

is primitive, compact, and Cover14.  At every shallow point j/13, W has a
two-owner local maximum of height 1/13, the two owners sum to the common ruler
182, and the deleted runner 182 is exactly at zero.  Nevertheless V is loose:
M(V)=2/17, while M(W)=2/13.  Thus local shallow witnesses plus their common
coupling ruler and the deleted-runner blocking data do not imply a strict
cover, an AP deletion, or crown collapse.

All arithmetic is exact.  The maximum evaluator uses the standard
piecewise-linear candidate set: opposite-slope pair intersections have
denominator v+w, same-slope changes use |v-w|, and self cusps use 2v.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import gcd


W = (14, 28, 29, 56, 70, 84, 98, 112, 126, 153, 154, 168)
Q = 182
V = tuple(sorted(W + (Q,)))
THRESHOLD = F(1, 13)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def dist_num(v: int, p: int, q: int) -> int:
    r = (v * p) % q
    return min(r, q - r)


def margin(values: tuple[int, ...], t: F) -> F:
    return min(F(dist_num(v, t.numerator, t.denominator), t.denominator) for v in values)


def exact_max(values: tuple[int, ...]) -> tuple[F, tuple[F, ...], int]:
    denoms: set[int] = {2 * v for v in values}
    for i, a in enumerate(values):
        for b in values[i + 1 :]:
            denoms.add(a + b)
            denoms.add(abs(a - b))
    denoms.discard(0)

    best = F(0)
    args: set[F] = set()
    candidates = 0
    for q in sorted(denoms):
        for p in range(q):
            candidates += 1
            t = F(p, q)
            value = margin(values, t)
            if value > best:
                best = value
                args = {t}
            elif value == best:
                args.add(t)
    return best, tuple(sorted(args)), candidates


def cover_carriers(values: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple((q, min(v for v in values if v % q == 0)) for q in range(2, 15))


def shallow_rows() -> tuple[tuple[int, int, int, int], ...]:
    rows = []
    epsilon = F(1, 13 * 1000 * max(V))
    for j in range(1, 13):
        residues = {w: (w * j) % 13 for w in W}
        plus = [w for w, r in residues.items() if r == 1]
        minus = [w for w, r in residues.items() if r == 12]
        require(len(plus) == len(minus) == 1, f"bad active-owner count at j={j}")
        wp, wm = plus[0], minus[0]
        t = F(j, 13)
        require(margin(W, t) == THRESHOLD, f"bad shallow margin at j={j}")
        require(wp + wm == Q, f"bad complementary sum at j={j}")
        require(
            margin(W, t - epsilon) == THRESHOLD - wp * epsilon,
            f"bad left local germ at j={j}",
        )
        require(
            margin(W, t + epsilon) == THRESHOLD - wm * epsilon,
            f"bad right local germ at j={j}",
        )
        require(margin((Q,), t) == 0, f"deleted runner does not kill j={j}")
        rows.append((j, wp, wm, wp + wm))
    return tuple(rows)


def main() -> None:
    require(len(V) == 13 and len(set(V)) == 13, "V is not a thirteen-set")
    require(gcd(*V) == 1, "V is not primitive")
    require(
        tuple(sorted(v % 13 for v in W)) == tuple(range(1, 13)),
        "W is not a full nonzero residue packet",
    )
    require(all(W[r - 1] % 13 == r for r in range(1, 13)), "residue labels are wrong")
    require(
        all(W[r - 1] + W[12 - r] == Q for r in range(1, 7)),
        "complementary pair sums are not constant",
    )
    require(len(set(W[i + 1] - W[i] for i in range(11))) > 1, "W is an AP")

    carriers = cover_carriers(V)
    require(tuple(q for q, _ in carriers) == tuple(range(2, 15)), "Cover14 failed")
    ratio = F(V[-1], V[-2])
    require(ratio == F(91, 84) < 13, "compact ratio failed")

    local = shallow_rows()
    mw, aw, cw = exact_max(W)
    mv, av, cv = exact_max(V)
    require(mw == F(2, 13), "wrong exact value for W")
    require(mv == F(2, 17), "wrong exact value for V")
    require(F(3, 119) in av, "named global escape is not an argmax")
    require(margin(V, F(3, 119)) == F(2, 17), "wrong named escape margin")

    print("LOCAL-RULER DELETION-SURROGATE EXACT REFEREE")
    print(f"W={W}")
    print(f"V={V}")
    print(f"primitive={gcd(*V) == 1} distinct={len(set(V)) == 13}")
    print(f"Cover14 carriers={carriers}")
    print(f"compact top ratio={ratio} < 13")
    print("W is non-AP: sorted gaps are not constant")
    print("shallow rows (j, +1 owner, -1 owner, common sum):")
    for row in local:
        print(f"  {row}")
    print("at every j/13: M_local(W)=1/13, owners cross on ruler 182,")
    print("and the deleted 182-runner has phase zero")
    print(f"exact M(W)={mw}; argmax_count={len(aw)}; candidates={cw}")
    print(f"exact M(V)={mv}; argmax_count={len(av)}; candidates={cv}")
    print(f"explicit global escape t=3/119 has margin {margin(V, F(3, 119))}")
    print("CONCLUSION: the shallow local-ruler/deleted-runner data are not an")
    print("AP-extraction or crown-collapse certificate; global safe-set confinement")
    print("and the full edge/curtain atlas remain necessary.")
    print("ALL EXACT CHECKS PASS")


if __name__ == "__main__":
    main()
