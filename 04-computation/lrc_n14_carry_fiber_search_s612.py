#!/usr/bin/env python3
"""
lrc_n14_carry_fiber_search_s612.py

monad-compute-2026-06-03-S3

Bounded exhaustive carry-fiber search for the n=14 LRC lift/CRT
conservativity problem.  Direct compute follow-on to S611 / HYP-2167.

Context
-------
S611 (lrc_n14_carry_conservativity_s611.py, HYP-2167) probed the carry fiber
of v = r + 27k (27 == -1 mod 14) over the two primitive least-positive floor
shadows AP and V*.  It tested ONLY isolated wraps: local carry perturbations of
Hamming weight 1 and 2, each a single +27 wrap.  All were strict (M > 1/14),
and it concluded:

  "Any new floor lift must use a globally coherent carry pattern, not an
   isolated wrap."

This left the obvious computational question open: is there ANY nonzero carry
vector k -- not just weight 1 or 2 single wraps -- that drives the lifted row
back to the floor M = 1/14 (or below)?  HYP-2167's proof target conjectures
NO, except for carry patterns cohomologous to the AP / V* scalar carries.

Method
------
We search two complementary bounded neighbourhoods of each floor shadow R in
{AP, V*}, lifting to the row R + 27*k and computing the EXACT loneliness
radius M(R + 27k):

  1. L1-budget search: ALL carry vectors k in Z_{>=0}^{13} with sum(k) <= L1MAX.
     This covers large-magnitude wraps on a few coordinates (the "coherent"
     patterns the S611 weight-1/2 sweep could not reach) and is the natural
     generalisation of bounded carry budget.

  2. Magnitude-1 lattice: ALL k in {0,1}^13 (8192 rows per shadow).  This is the
     full Boolean carry lattice; it subsumes S611's weight-1 and weight-2 single
     wraps and extends to every Hamming weight 0..13 at magnitude 1.

For each lifted row we record M (exact Fraction).  We report:
  * how many nonzero-k rows land AT the floor (the conservativity question),
  * how many land BELOW the floor (would be an LRC(14) counterexample if any),
  * the minimum M over nonzero k and the witness,
  * the per-Hamming-weight min-M profile.

exact_maximin is imported unchanged from the S611 module so the maximin oracle
is the same VERIFIED code path (no reintroduced arithmetic).  We re-assert that
the two base shadows AP, V* and the scalar 2*AP are at the floor as a sanity
gate before trusting any "strict" verdict.

This is a compute-node diagnostic: it produces data on the carry fiber.  It
does NOT prove the carry-conservativity theorem; it bounds the radius in which
a counterexample would have to live.
"""

from __future__ import annotations

import sys
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

# Reuse the VERIFIED maximin oracle and shadows from the S611 module.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lrc_n14_carry_conservativity_s611 import (  # noqa: E402
    AP,
    VSTAR,
    C,
    FLOOR,
    N,
    exact_maximin,
    fmt_frac,
    shell_orbits,
)

BASE_ROWS = (("AP", AP), ("V*", VSTAR))
L1MAX = 6  # L1 carry budget for search (1)


@dataclass(frozen=True)
class FiberHit:
    base_name: str
    carry: tuple[int, ...]
    row: tuple[int, ...]
    m: Fraction
    t: Fraction

    @property
    def weight(self) -> int:
        return sum(1 for k in self.carry if k)

    @property
    def l1(self) -> int:
        return sum(self.carry)

    @property
    def label(self) -> str:
        moved = [(i, k) for i, k in enumerate(self.carry) if k]
        body = ",".join(f"i{i}:+{k}" for i, k in moved) if moved else "ZERO"
        return f"{self.base_name}:carry[{body}]"


def lift_row(base: tuple[int, ...], carry: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(b + C * k for b, k in zip(base, carry)))


FLOOR_F = float(FLOOR)


def _fdist(x: float) -> float:
    """Distance from x to the nearest integer (the circle metric)."""
    r = x - int(x)
    if r < 0:
        r += 1.0
    return r if r <= 0.5 else 1.0 - r


def approx_maximin(row: tuple[int, ...]) -> float:
    """
    Fast float maximin: M = max_t min_i ||v_i t|| over the same critical-time
    candidate set used by the exact oracle (single-runner midpoints (2m+1)/(2a)
    and pairwise fractions m/(a+b), m/|a-b|).  No Fraction construction, no set
    dedup -- candidates are streamed.  Used only to SCREEN rows; any row whose
    approximate M lands near or below the floor is re-checked exactly.
    """
    best = 0.0
    n = len(row)
    for i in range(n):
        a = row[i]
        inv2a = 1.0 / (2 * a)
        for m in range(a):
            t = (2 * m + 1) * inv2a
            if 0.0 < t < 1.0:
                s = 1.0
                for v in row:
                    d = _fdist(v * t)
                    if d < s:
                        s = d
                        if s <= best:
                            break
                if s > best:
                    best = s
        for j in range(i + 1, n):
            b = row[j]
            for den in (a + b, a - b if a > b else b - a):
                if den <= 0:
                    continue
                invden = 1.0 / den
                for m in range(1, den):
                    t = m * invden
                    s = 1.0
                    for v in row:
                        d = _fdist(v * t)
                        if d < s:
                            s = d
                            if s <= best:
                                break
                    if s > best:
                        best = s
    return best


def l1_carry_vectors(dim: int, budget: int):
    """Yield every nonneg integer vector of length dim with sum <= budget."""
    vec = [0] * dim

    def rec(pos: int, remaining: int):
        if pos == dim:
            yield tuple(vec)
            return
        for k in range(remaining + 1):
            vec[pos] = k
            yield from rec(pos + 1, remaining - k)
        vec[pos] = 0

    yield from rec(0, budget)


def boolean_carry_vectors(dim: int):
    for mask in range(1 << dim):
        yield tuple((mask >> i) & 1 for i in range(dim))


def sanity_gate() -> None:
    print("Sanity gate (must all be floor=1/14):")
    for name, row in BASE_ROWS:
        m, t = exact_maximin(row)
        ok = "OK" if m == FLOOR else "*** FAIL ***"
        print(f"  {name:3s} M={fmt_frac(m)} t={fmt_frac(t)} {ok}")
    two_ap = tuple(sorted(2 * v for v in AP))
    m, t = exact_maximin(two_ap)
    ok = "OK" if m == FLOOR else "*** FAIL ***"
    print(f"  2*AP M={fmt_frac(m)} t={fmt_frac(t)} {ok}  (nonprimitive)")
    print()


TOL = 1e-6  # float screening tolerance; distinct rationals differ from 1/14 by >> this


def _exact_hit(name: str, base, carry) -> FiberHit:
    row = lift_row(base, carry)
    m, t = exact_maximin(row)
    return FiberHit(base_name=name, carry=carry, row=row, m=m, t=t)


def run_search(name: str, base: tuple[int, ...], vectors, label: str) -> dict:
    total = 0
    near_floor: list[tuple] = []   # carry vectors with float M <= floor+TOL -> exact check
    weight_min_f: dict[int, float] = {}     # per-weight min float M
    weight_min_carry: dict[int, tuple] = {}
    min_f = float("inf")
    min_f_carry: tuple | None = None
    zero_m = None

    for carry in vectors:
        total += 1
        row = lift_row(base, carry)
        is_zero = all(k == 0 for k in carry)
        mf = approx_maximin(row)
        w = sum(1 for k in carry if k)
        if w not in weight_min_f or mf < weight_min_f[w]:
            weight_min_f[w] = mf
            weight_min_carry[w] = carry
        if is_zero:
            zm, _ = exact_maximin(row)
            zero_m = zm
            continue
        if mf < min_f:
            min_f = mf
            min_f_carry = carry
        if mf <= FLOOR_F + TOL:
            near_floor.append(carry)

    # Exact verification of every near-floor row (catches all floor / below-floor).
    floor_nonzero: list[FiberHit] = []
    below: list[FiberHit] = []
    for carry in near_floor:
        hit = _exact_hit(name, base, carry)
        if hit.m == FLOOR:
            floor_nonzero.append(hit)
        elif hit.m < FLOOR:
            below.append(hit)
    # Exact value for the global-min and per-weight-min winners.
    min_hit = _exact_hit(name, base, min_f_carry) if min_f_carry is not None else None
    weight_min_exact = {
        w: _exact_hit(name, base, weight_min_carry[w]).m for w in sorted(weight_min_carry)
    }

    print(f"[{label}] base={name} rows={total}")
    if zero_m is not None:
        print(f"  k=0 baseline M={fmt_frac(zero_m)} (expect 1/14)")
    print(f"  nonzero-k rows AT floor (M=1/14): {len(floor_nonzero)}")
    print(f"  nonzero-k rows BELOW floor (M<1/14): {len(below)}")
    print(f"  rows exact-verified near floor (float<=1/14+{TOL:g}): {len(near_floor)}")
    if min_hit is not None:
        print(
            f"  min M over nonzero k: {fmt_frac(min_hit.m)} "
            f"(weight={min_hit.weight}, L1={min_hit.l1}) "
            f"via {min_hit.label} t={fmt_frac(min_hit.t)}"
        )
    if floor_nonzero:
        print("  *** NONZERO FLOOR CARRIES FOUND ***")
        for hit in floor_nonzero[:20]:
            print(f"    {hit.label} row={hit.row} M={fmt_frac(hit.m)}")
    if below:
        print("  *** BELOW-FLOOR CARRIES FOUND (would refute LRC(14)) ***")
        for hit in below[:20]:
            print(f"    {hit.label} row={hit.row} M={fmt_frac(hit.m)}")
    prof = ", ".join(f"w{w}:{fmt_frac(weight_min_exact[w])}" for w in sorted(weight_min_exact))
    print(f"  min-M by Hamming weight (exact at the min witness): {prof}")
    print()

    return {
        "total": total,
        "floor_nonzero": len(floor_nonzero),
        "below": len(below),
        "near_floor": len(near_floor),
        "min_nonzero_m": None if min_hit is None else min_hit.m,
    }


def main() -> None:
    print(f"==== LRC n={N} carry-fiber bounded search (S612) ====")
    print(f"C=2n-1={C}; floor=1/{N}; lift row = shadow + {C}*k")
    print(f"L1 carry budget L1MAX={L1MAX}")
    print()

    print("THM-407 shell orbits under G=<2,-1> modulo 27 (context):")
    for orbit in shell_orbits():
        from math import gcd

        g = gcd(orbit[0], C)
        print(f"  gcd={g:2d}: {orbit}")
    print()

    sanity_gate()

    dim = len(AP)  # = N-1 = 13
    summary = {}

    # Count L1 vectors once for reporting.
    print("Search 1: L1-budget carry vectors (sum(k) <= L1MAX), all magnitudes")
    for name, base in BASE_ROWS:
        res = run_search(
            name,
            base,
            l1_carry_vectors(dim, L1MAX),
            f"L1<= {L1MAX}",
        )
        summary[(name, "L1")] = res

    print("Search 2: Boolean carry lattice {0,1}^13 (all Hamming weights, mag 1)")
    for name, base in BASE_ROWS:
        res = run_search(
            name,
            base,
            boolean_carry_vectors(dim),
            "bool{0,1}^13",
        )
        summary[(name, "bool")] = res

    print("Conclusion:")
    any_floor = any(r["floor_nonzero"] for r in summary.values())
    any_below = any(r["below"] for r in summary.values())
    if any_below:
        print("  *** A below-floor carry lift was found -- CHECK IMMEDIATELY. ***")
    elif any_floor:
        print(
            "  Nonzero carry vectors that stay AT the floor exist within the "
            "searched radius -- inspect them as candidate 'globally coherent' "
            "patterns (HYP-2167 case 1)."
        )
    else:
        print(
            "  No nonzero carry vector within the searched radius (L1<=%d over "
            "AP and V*, plus the full {0,1}^13 lattice) returns the lifted row "
            "to the floor M=1/14." % L1MAX
        )
        print(
            "  This strengthens HYP-2167: in this bounded fiber neighbourhood the "
            "least-positive section AP/V* is the ONLY floor representative; every "
            "isolated OR coherent bounded carry is strictly above the floor.  Any "
            "new floor lift must live outside this radius (large coherent carry) "
            "or route through a multiple/owner certificate."
        )


if __name__ == "__main__":
    main()
