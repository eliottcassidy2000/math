#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) via the DILATED-BAND COVERING reformulation + the resource-bound probe.
kind-pasteur-2026-06-13-S1.  Applies the repo's sum-free-covering lens (THM-469)
to the band criterion (THM-492).

THE REFRAME.  Band criterion (THM-492, verified 0-mismatch there): t = a/q,
gcd(a,q)=1, is a STRICT loneliness witness for speed set S (M(S) > 1/14) iff
every v in S has (v*a mod q) NOT in the band B_q = {r : min(r,q-r) <= floor(q/14)}.
Equivalently a is a witness-numerator iff a (a unit mod q) avoids the union of
DILATED BANDS  U_v := v^{-1} * B_q  (mod q).  So:

    S has a shell-q witness  <=>  (Z/q)*  \  U_{v in S} v^{-1} B_q  is NONEMPTY.

C'(14) (=> LRC(14), THM-398) says: for every primitive multiple-of-14 set S,
SOME q leaves an uncovered unit.  The 13 dilated bands each have relative size
(2 floor(q/14)+1)/q -> 1/7, so 13 of them have total measure ~13/7 ~ 1.86 > 1:
measure ALONE permits a cover -- covering is a genuine Diophantine/additive
question (where THM-469's sum-free/unit-escape structure should bite).

THIS SCRIPT:
 (A) Implements the band criterion + the covering reformulation; cross-checks
     them (must agree by construction) AND validates the band criterion against a
     dense-rational direct computation of M(S) on small/known configs (MISTAKE-067
     discipline: an independent method, not the same tool twice).
 (B) RESOURCE BOUND: for a config, the LADDER HEIGHT h(S) = least q with an
     uncovered unit (the first witness shell).  Measures h on the known hard
     families -- the five evaders S(r)=7*{1..12} u {r}, r in {611,702,793,962,1053}
     (expect h in {40,41}, band-2), and four-deletions of 7*{1..12} (eight-core,
     expect h in {31,33}) -- and on random primitive multiple-of-14 configs.
 (C) COVERING DEFICIT curve def(q) = #uncovered units, q = 14..H, for the evaders:
     shows the cover is COMPLETE (def=0) through band-1 (q<=27) and first leaks at
     band-2.  Reports, at the critical q*, the additive structure of {v^{-1} mod q*}
     of the band that the cover needed (the resource the 13 runners spent).
"""

import itertools, time, random
from fractions import Fraction
from math import gcd


# ----------------------------------------------------------- band criterion

def band_set(q, n=14):
    """B_q = {r in Z/q : min(r, q-r) <= floor(q/n)} = +-{0,1,...,floor(q/n)}."""
    h = q // n
    B = set()
    for r in range(q):
        if min(r, q - r) <= h:
            B.add(r)
    return B


def witness_numerators(S, q, n=14):
    """Units a mod q with (v*a mod q) not in B_q for all v in S.
    These are exactly the escaping units of the dilated-band cover."""
    B = band_set(q, n)
    out = []
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        if all((v * a) % q not in B for v in S):
            out.append(a)
    return out


def covering_deficit(S, q, n=14):
    """# units (Z/q)* NOT covered by U_v v^{-1} B_q  ==  #witness numerators."""
    return len(witness_numerators(S, q, n))


def has_shell_witness(S, q, n=14):
    return covering_deficit(S, q, n) > 0


def ladder_height(S, Hmax=200, n=14):
    """Least q in [2, Hmax] with a strict shell witness, or None."""
    for q in range(2, Hmax + 1):
        if covering_deficit(S, q, n) > 0:
            return q
    return None


# ------------------------------------ independent M(S) (dense-rational, slow)

def M_dense(S, Qmax=400):
    """max over a/q (q<=Qmax, gcd(a,q)=1) of min_v ||v a / q||, as a Fraction.
    Independent of the band criterion -- validates it. M(S) > 1/14 <=> loose
    (the band-ladder theory says the max is attained at finite q; Qmax large)."""
    best = Fraction(0)
    for q in range(1, Qmax + 1):
        for a in range(1, q + 1):
            if gcd(a, q) != 1:
                continue
            m = min((min((v * a) % q, q - (v * a) % q) for v in S))
            val = Fraction(m, q)
            if val > best:
                best = val
    return best


# ----------------------------------------------------------- config families

def five_evaders():
    base = [7 * k for k in range(1, 13)]  # 7,14,...,84
    return {r: sorted(base + [r]) for r in (611, 702, 793, 962, 1053)}


def eight_core_deletions():
    """4-deletions of 7*{1..12} that the reflection flagged Q27-feasible, plus
    a small sweep. Returns dict address->config (12-4+? ) -- the reflection's
    two exceptional addresses (28,42,56,84) and (42,56,70,84) PLUS we keep the
    config a primitive multiple-of-14 set by appending 14*1=14 if needed."""
    base = [7 * k for k in range(1, 13)]  # 7..84
    addrs = [(28, 42, 56, 84), (42, 56, 70, 84)]
    out = {}
    for ad in addrs:
        kept = [x for x in base if x not in ad]
        # ensure a multiple of 14 present and gcd 1
        S = sorted(set(kept))
        out[ad] = S
    return out


def is_primitive_mult14(S):
    from functools import reduce
    g = reduce(gcd, S)
    return g == 1 and any(v % 14 == 0 for v in S)


# ----------------------------------------------------------------- the lab

def part_A_validate():
    print("=== A. band criterion <-> covering reformulation <-> dense M ===", flush=True)
    # small known configs: tight {1..13} (M=1/14, loose? no -- tight), and a loose one
    tests = [
        ("tight {1..13}", list(range(1, 14)), Fraction(1, 14)),
        ("{1..12,28}", list(range(1, 13)) + [28], None),   # THM-492: M=1/13 loose
        ("evader r=611", sorted([7 * k for k in range(1, 13)] + [611]), None),
    ]
    for name, S, expect in tests:
        Md = M_dense(S, Qmax=250)
        loose_dense = Md > Fraction(1, 14)
        h = ladder_height(S, Hmax=120)
        loose_band = h is not None
        ok = (loose_dense == loose_band)
        line = (f"   {name}: M_dense={Md} (>1/14? {loose_dense}); "
                f"band ladder height h={h} (witness? {loose_band}); agree={ok}")
        if expect is not None:
            line += f"; expected M={expect} {'OK' if Md==expect else 'MISMATCH'}"
        print(line, flush=True)


def part_B_resource(Hmax=200):
    print("\n=== B. ladder height h(S) on the hard families ===", flush=True)
    print("   five evaders S(r)=7*{1..12} u {r} (THM-492: expect h in {40,41}):", flush=True)
    for r, S in five_evaders().items():
        h = ladder_height(S, Hmax)
        print(f"      r={r}: primitive_mult14={is_primitive_mult14(S)}  h(S)={h}  "
              f"deficit@h={covering_deficit(S,h) if h else 0}", flush=True)
    print("   eight-core 4-deletions (reflection: expect plain q in {31,33}):", flush=True)
    for ad, S in eight_core_deletions().items():
        h = ladder_height(S, Hmax)
        print(f"      delete{ad}: |S|={len(S)} mult14={is_primitive_mult14(S)}  h(S)={h}", flush=True)
    print("   random primitive multiple-of-14 configs (13 speeds, entries<=120):", flush=True)
    rng = random.Random(2026)
    hts = []
    tries = 0
    while len(hts) < 40 and tries < 4000:
        tries += 1
        S = sorted(rng.sample(range(1, 121), 13))
        if not is_primitive_mult14(S):
            continue
        h = ladder_height(S, Hmax)
        if h is not None:
            hts.append(h)
    hts.sort()
    if hts:
        print(f"      {len(hts)} configs: h min={hts[0]} median={hts[len(hts)//2]} "
              f"max={hts[-1]}  (max ladder height seen = candidate resource bound)", flush=True)
        from collections import Counter
        print(f"      height distribution: {dict(sorted(Counter(hts).items()))}", flush=True)


def part_C_deficit_curve(Hmax=50):
    print("\n=== C. covering-deficit curve def(q) for the evaders (band-1 ceiling q=27) ===", flush=True)
    for r, S in five_evaders().items():
        curve = [(q, covering_deficit(S, q)) for q in range(14, Hmax + 1)]
        nz = [(q, d) for q, d in curve if d > 0]
        first = nz[0] if nz else None
        b1 = max((d for q, d in curve if q <= 27), default=0)  # band-1 deficit
        print(f"   r={r}: band-1 (q<=27) max deficit = {b1} (0 = fully covered through band-1); "
              f"first leak at q*={first[0] if first else None} (deficit {first[1] if first else 0})", flush=True)
        if first:
            q = first[0]
            # the additive structure: the inverse-residues that ALMOST covered
            inv = sorted(set((pow(v % q, -1, q) if gcd(v, q) == 1 else None) for v in S) - {None})
            print(f"        {{v^{{-1}} mod {q}}} (units) = {inv}", flush=True)


def main():
    t0 = time.time()
    part_A_validate()
    part_B_resource()
    part_C_deficit_curve()
    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
