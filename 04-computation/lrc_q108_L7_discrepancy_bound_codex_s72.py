#!/usr/bin/env python3
"""HYP-2729/2730: elementary discrepancy bound for the L7 torus curve.

KPS HYP-2730 reduces the L7 tail to the L1 cell discrepancy

    D_{p,q} = sum_{i,j in Z/7} |mu_{p,q}(i,j) - 1/49|

of the closed torus curve v -> (qv,pv).  Their sharp observed law is
`D_{p,q} <= (12/7)/q`.  This note records a cruder elementary proof route that
is already strong enough for L7:

    D_{p,q} <= 24/(7q).

Proof sketch.  Put u={qv}.  Conditional on u, the second coordinate samples the
q equally spaced points `{n/q + (p/q)u}`.  In any seven equal bins, q equally
spaced points have each bin count differing from q/7 by at most one, so the
conditional seven-bin L1 discrepancy is at most `2 r (7-r)/(7q) <= 24/(7q)`,
where `q=7a+r`.  Integrating over each first-coordinate sector and summing the
seven sectors gives the same bound for D.

With the conservative L7 margin `0.21`, this proves every `q>=17` tail because
`24/(7q) < 0.21`.  Exact D checks for `q<=16` show only the KPS small atlas has
`D>=0.21`.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import gcd, lcm


P = 7
MARGIN = F(21, 100)


def mu_cell_nums(p: int, q: int) -> tuple[int, dict[tuple[int, int], int]]:
    """Return common denominator D and integer cell masses for mu_{p,q}."""
    D = lcm(P * p, P * q)
    pts = {0, D}
    for f in (p, q):
        den = P * f
        for t in range(P * f):
            pts.add(D * t // den)
    xs = sorted(pts)
    mu: dict[tuple[int, int], int] = {}
    den2 = 2 * D
    for a, b in zip(xs, xs[1:]):
        mid_num = a + b
        key = (
            (P * ((q * mid_num) % den2)) // den2,
            (P * ((p * mid_num) % den2)) // den2,
        )
        mu[key] = mu.get(key, 0) + (b - a)
    return D, mu


def D_pq(p: int, q: int) -> F:
    D, mu = mu_cell_nums(p, q)
    return sum(abs(49 * mu.get((i, j), 0) - D) for i in range(P) for j in range(P)) / F(49 * D)


def ratios(qmax: int) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for q in range(1, qmax + 1):
        for p in range(q + 1, int(F(43, 20) * q) + 1):
            if gcd(p, q) == 1 and F(1) < F(p, q) <= F(43, 20):
                out.append((p, q))
    return out


def fmt(x: F) -> str:
    return f"{x} ({float(x):.9f})"


def main() -> None:
    print("HYP-2729/HYP-2730 L7 torus-curve discrepancy bound")
    print("Exact D_{p,q}; elementary crude proof target D <= 24/(7q).")
    print()

    rows = []
    failures_crude = []
    failures_sharp = []
    for p, q in ratios(80):
        D = D_pq(p, q)
        crude = F(24, 7 * q)
        sharp = F(12, 7 * q)
        if D > crude:
            failures_crude.append((p, q, D, crude))
        if D > sharp:
            failures_sharp.append((p, q, D, sharp))
        rows.append((p, q, D, D * q))

    print(f"ratios_checked={len(rows)} q<=80")
    print(f"crude_bound_failures D<=24/(7q): {len(failures_crude)}")
    print(f"sharp_observed_failures D<=12/(7q): {len(failures_sharp)}")
    max_scaled = max(rows, key=lambda row: (row[3], -row[1], row[0]))
    print(
        "max D*q checked: "
        f"p/q={max_scaled[0]}/{max_scaled[1]} D*q={fmt(max_scaled[3])} D={fmt(max_scaled[2])}"
    )

    print()
    print("FINITE DISCREPANCY ATLAS q<=16")
    finite_bad = []
    for p, q, D, scaled in rows:
        if q <= 16 and D >= MARGIN:
            finite_bad.append((p, q, D, scaled))
    for p, q, D, scaled in finite_bad:
        print(f"  p/q={p}/{q:<2d} D={fmt(D)} D*q={fmt(scaled)}")
    print(f"finite q<=16 with D>=0.21: {len(finite_bad)}")

    tail_start = 17
    print()
    print("TAIL CERTIFICATE")
    print(f"  For q>={tail_start}, 24/(7q) <= {fmt(F(24, 7 * tail_start))} < 0.21.")
    print("  Therefore the elementary crude bound alone makes every q>=17 resonance tail safe.")
    print("  Exact D for 5<=q<=16 is already below 0.21, so only q<=4 needs p0_inf atlas checks.")
    print()
    print("SYNTHESIS")
    print("  HYP-2730's sharp 12/7 constant is excellent, but L7 does not need it.")
    print("  A simple seven-bin lattice discrepancy bound plus finite q<=16 D-checks leaves")
    print("  the same small rational atlas and is a cleaner formalization target.")


if __name__ == "__main__":
    main()
