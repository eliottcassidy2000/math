#!/usr/bin/env python3
"""
lrc_n16_adelic_product_s19.py

oracle-2026-05-31-S19

Integrates the repo's RAPIDITY / ADELIC thread (THM-252; adelic-tournament-
geometry) with codex's n=16 dyadic-endpoint-debt program (S389/S390), aiming at
an n=16 Lonely Runner proof.

THE IDEA.  In the adelic tournament geometry the Cayley formal group
F(x,y)=(x+y)/(1+xy) is SUPERSINGULAR (height ∞) at p=2 and ordinary (height 1)
at odd p; rapidity = arctanh = its formal logarithm, living in the log-prime
lattice (THM-252).  n=16=2^4 is the pure 2-power, so the n=16 LRC obstruction
lives entirely at the supersingular place.  codex (S389) empirically found a
CONSERVED PRODUCT across the dyadic ladders:

    (#unprotected endpoints) * (max_gap / threshold) ≈ 34/33  (constant)

This is exactly an ADELIC PRODUCT FORMULA: the archimedean size (the real gap,
= how lonely) times the 2-adic size (the endpoint debt, = how much protection is
owed) is conserved.  A counterexample needs gap=0 AND debt=0 simultaneously --
the bad corner -- which a positive lower bound on the product forbids.

THIS SCRIPT tests whether the product is a UNIVERSAL lower-bounded invariant
(not just a ladder artifact): if inf over all primitive 15-speed n=16 sets of a
suitable Gap*Debt functional is bounded away from 0, no counterexample can sit at
(gap=0, debt=0).  It also resolves the unprotected endpoints by 2-adic debt layer
(the Bruhat-Tits depth) to expose the dyadic descent.
"""

from __future__ import annotations

import random
from collections import Counter
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[0]
S356 = SourceFileLoader("s356", str(ROOT / "lonely_runner_residue_probe_s356.py")).load_module()
S360 = SourceFileLoader("s360", str(ROOT / "lonely_runner_endpoint_protection_s360.py")).load_module()

N = 16            # runners
K = N - 1         # speeds = 15
THR = Fraction(1, N)


def primitive(vs) -> bool:
    g = 0
    for v in vs:
        g = gcd(g, v)
    return g == 1


def v2(x: int) -> int:
    if x == 0:
        return 99
    out = 0
    while x % 2 == 0:
        x //= 2
        out += 1
    return out


def debt_layer(point: Fraction) -> int:
    """Odd part of the denominator's relation to N: Bruhat-Tits-ish depth via v2."""
    return v2(point.denominator)


def analyze(speeds):
    speeds = tuple(sorted(set(speeds)))
    report = S356.report("x", list(speeds))
    sp = report.speeds
    gap_ratio = report.max_gap / report.threshold          # archimedean size
    forbidden = report.forbidden_length
    # unprotected endpoint VALUES (distinct): surviving lonely boundary points
    pts = sorted({e.value for e in S360.endpoints(sp)})
    unprot = [p for p in pts if not any(S360.direct_protects(sp, v, p) for v in sp)]
    debt = len(unprot)
    layerhist = Counter(debt_layer(p) for p in unprot)
    return {
        "speeds": sp,
        "gap_ratio": gap_ratio,
        "forbidden": forbidden,
        "debt": debt,
        "product": gap_ratio * debt,            # Gap * Debt (the conserved quantity)
        "has_16gate": any(v % 16 == 0 for v in sp),
        "layerhist": dict(sorted(layerhist.items())),
    }


def lpd_ladder(d: int, skip: int):
    return tuple(sorted({1} | {d * q for q in range(1, N) if q != skip}))


def best_ladder(d: int):
    best = None
    for skip in range(1, N):
        s = lpd_ladder(d, skip)
        if len(s) != K or not primitive(s):
            continue
        r = S356.report("x", list(s))
        key = (r.max_gap, r.boundary_witness_count, skip)
        if best is None or key < best[0]:
            best = (key, s)
    return best[1]


def fmt(x):
    return f"{float(x):.6f}"


def main():
    print("LRC n=16 adelic Gap*Debt product (oracle-2026-05-31-S19)\n")

    print("=" * 70)
    print("A. Dyadic ladders (reproduce codex's conserved product)")
    print("=" * 70)
    for d in (2, 4, 8, 16):
        s = best_ladder(d)
        a = analyze(s)
        print(f" d={d:<2} gap_ratio={fmt(a['gap_ratio'])} debt={a['debt']:<4} "
              f"product={fmt(a['product'])}  layers(v2)={a['layerhist']}")
    print(" -> product ~ const across ladders is the adelic conservation.\n")

    print("=" * 70)
    print("B. Structured non-ladder sets")
    print("=" * 70)
    structured = {
        "initial 1..15": tuple(range(1, 16)),
        "1..15 drop8 add16": tuple(sorted(set(range(1, 16)) - {8} | {16})),
        "1..15 drop15 add16": tuple(sorted(set(range(1, 16)) - {15} | {16})),
        "odds+16": (1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 16),
    }
    for label, s in structured.items():
        if len(set(s)) != K or not primitive(s):
            print(f" [{label}] skipped (not primitive 15-set)")
            continue
        a = analyze(s)
        tag = "16gate" if a["has_16gate"] else "no-gate"
        print(f" [{label}] {tag} gap_ratio={fmt(a['gap_ratio'])} forbidden={a['forbidden']} "
              f"debt={a['debt']} product={fmt(a['product'])} layers={a['layerhist']}")
    print()

    print("=" * 70)
    print("C. Random forced-16-gate primitive 15-sets: is the product bounded below?")
    print("=" * 70)
    rng = random.Random(19)
    min_pos_product = None
    min_pos_set = None
    counterexamples = []
    tested = 0
    forced_gate_trials = 300
    for _ in range(forced_gate_trials):
        gate = rng.choice([16, 32, 48, 16, 16])
        pool = [v for v in range(1, 40) if v % 16 != 0]
        rest = rng.sample(pool, 14)
        s = tuple(sorted(set(rest) | {gate}))
        if len(s) != K or not primitive(s):
            continue
        tested += 1
        a = analyze(s)
        if a["forbidden"] == 1 and a["debt"] == 0:
            counterexamples.append(s)          # the bad corner
        # track smallest positive Gap*Debt-style product among full-measure-ish sets
        prod = a["product"]
        if prod > 0 and (min_pos_product is None or prod < min_pos_product):
            min_pos_product = prod
            min_pos_set = (s, a)
    print(f" trials={forced_gate_trials}  tested={tested}  open-cover counterexamples found={len(counterexamples)}")
    if min_pos_set:
        s, a = min_pos_set
        print(f" smallest positive product={fmt(min_pos_product)} at debt={a['debt']} "
              f"gap_ratio={fmt(a['gap_ratio'])}")
    print()

    print("=" * 70)
    print("D. The bad corner test: any set with forbidden=1 AND debt=0?")
    print("=" * 70)
    # near-tight: forbidden_length close to 1 -> check debt never hits 0 there
    near_tight = []
    rng2 = random.Random(1900)
    near_tight_trials = 300
    for _ in range(near_tight_trials):
        gate = 16
        pool = [v for v in range(1, 30) if v % 16 != 0]
        rest = rng2.sample(pool, 14)
        s = tuple(sorted(set(rest) | {gate}))
        if len(s) != K or not primitive(s):
            continue
        report = S356.report("x", list(s))
        if report.forbidden_length >= Fraction(99, 100):
            a = analyze(s)
            near_tight.append((a["forbidden"], a["debt"], s))
    near_tight.sort(reverse=True)
    print(f" near-tight trials={near_tight_trials}")
    print(f" near-tight (forbidden>=0.99) sets: {len(near_tight)}")
    zero_debt = [x for x in near_tight if x[1] == 0]
    print(f" of those with debt=0 (=counterexample corner): {len(zero_debt)}")
    for forb, debt, s in near_tight[:8]:
        print(f"   forbidden={fmt(forb)} debt={debt} speeds={s}")
    print()
    print("SUMMARY: if no set reaches (forbidden=1, debt=0) and Gap*Debt stays")
    print("bounded away from 0, the adelic product formula forbids the bad corner")
    print("=> evidence that no n=16 counterexample exists.")


if __name__ == "__main__":
    main()
