#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S9b -- HYP-4352: THE RESIDUE-RIGIDITY toward transversal
AP-rigidity (the pinned-modulus crux).

From S9 (HYP-4342): a gap member (M in (1/13,2/25)) is NON-CLEARING mod every
q; for q in {13,...,25} that means the residues are PAIR-HITTING ({0,+-1}-hitting
under every rotation).  Question: is the residue structure FORCED further -- to
the AP's complete nonzero system?

TESTS:
 (A) COMPLETE-RESIDUE mod 13: among covering families (M > 1/13) non-clearing
     mod 13, are the residues mod 13 the COMPLETE nonzero system {1..12}?
     Construct families with prescribed mod-13 residue multisets and measure M.
 (B) THE INTERSECTION: pair-hitting mod EVERY q in 13..25 simultaneously -- how
     constrained is the residue structure?  What families (covering, primitive)
     survive besides the AP?
 (C) THE AP FINGERPRINT: characterize the AP's residue signature across moduli
     and test whether a gap member must match it.
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def dist_int(x, q):
    r = x % q
    return min(r, q - r)

def exact_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0); seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def float_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = 0.0
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = j / s
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def nonclearing_mod(W, q):
    thr = F(2, 25) * q
    residues = [v % q for v in W]
    for a in range(1, q):
        if all(dist_int(r * a, q) >= thr for r in residues):
            return False
    return True

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

RHO = F(2, 25); FLOOR = F(1, 13)

def part_A():
    print("=" * 78)
    print("PART A: does covering (M>1/13) + non-clearing mod 13 force COMPLETE")
    print("        residues {1..12} mod 13?")
    print("=" * 78)
    random.seed(91)
    # build families with controlled mod-13 residue completeness, measure M
    complete_stats = {"gap": 0, "floor": 0, "clear": 0, "below": 0}
    incomplete_stats = {"gap": 0, "floor": 0, "clear": 0, "below": 0}
    inc_covering = []   # incomplete-mod-13 families that are covering (M>1/13)
    tested = 0
    fl = float(FLOOR); rh = float(RHO)
    for _ in range(40000):
        W = sorted(random.sample(range(1, 70), 12))
        W = primitive(W)
        if len(W) != 12 or not nonclearing_mod(W, 13):
            continue
        tested += 1
        res13 = set(v % 13 for v in W)
        complete = (res13 == set(range(1, 13)))   # all nonzero residues
        fm = float_M(W)
        if fm > fl + 1e-6 and fm < rh - 1e-6:       # candidate gap or between
            M = exact_M(W)                          # exact only near/in gap
            bucket = "gap" if FLOOR < M < RHO else ("floor" if M == FLOOR else ("clear" if M >= RHO else "below"))
        elif abs(fm - fl) <= 1e-6:
            bucket = "floor"
        elif fm >= rh - 1e-6:
            bucket = "clear"
        else:
            bucket = "below"
        (complete_stats if complete else incomplete_stats)[bucket] += 1
        if not complete and bucket == "gap":
            inc_covering.append((W, bucket, sorted(res13)))
    print(f"  non-clearing-mod-13 families tested: {tested}")
    print(f"  COMPLETE residues {{1..12}} mod 13: {complete_stats}")
    print(f"  INCOMPLETE (misses a residue) mod 13: {incomplete_stats}")
    print(f"  => INCOMPLETE-mod-13 GAP members (M in (1/13,2/25)): {len(inc_covering)}")
    if inc_covering:
        print("    HYPOTHESIS 'gap member => complete mod 13' is FALSE; examples:")
        for W, b, r in inc_covering[:6]:
            print(f"      res13={r} W={list(W)}")
    else:
        print("    HYPOTHESIS SUPPORTED: every GAP member non-clearing mod 13 has")
        print("    COMPLETE nonzero residues mod 13 (the AP's mod-13 fingerprint).")
        print(f"    (incomplete-mod-13 families are only floor/clear/below, never in-gap)")

def part_B():
    print()
    print("=" * 78)
    print("PART B: families non-clearing mod EVERY q in 13..25 -- residue structure")
    print("=" * 78)
    random.seed(92)
    MODS = list(range(13, 26))
    survivors = []
    tested = 0
    # search AP-perturbations (the near-tight regime where survivors live)
    cands = set()
    cands.add(tuple(range(1, 13)))
    for _ in range(60000):
        W = list(range(1, 13))
        for _i in range(random.randint(1, 5)):
            idx = random.randrange(12)
            W[idx] += random.choice([13, 25, 26, 39, 50, 65]) * random.randint(1, 4)
        W = primitive(tuple(W))
        if len(W) == 12:
            cands.add(W)
    for W in cands:
        if len(set(W)) != 12:
            continue
        tested += 1
        if all(nonclearing_mod(W, q) for q in MODS):
            M = exact_M(W)
            res13 = tuple(sorted(set(v % 13 for v in W)))
            res25 = tuple(sorted(set(v % 25 for v in W)))
            survivors.append((W, M, res13, res25))
    print(f"  candidates: {tested};  non-clearing mod all of 13..25: {len(survivors)}")
    gap = [s for s in survivors if FLOOR < s[1] < RHO]
    print(f"  survivors in-gap: {len(gap)};  at floor 1/13: {sum(1 for s in survivors if s[1]==FLOOR)};  "
          f">=2/25: {sum(1 for s in survivors if s[1]>=RHO)}")
    # residue completeness among survivors
    complete13 = sum(1 for s in survivors if set(s[2]) == set(range(1,13)))
    print(f"  survivors with COMPLETE residues mod 13: {complete13}/{len(survivors)}")
    if gap:
        print("  *** IN-GAP survivors:")
        for W, M, r13, r25 in gap[:6]:
            print(f"    M={M} W={list(W)}")
    for W, M, r13, r25 in survivors[:8]:
        print(f"    M={M} res13={list(r13)} res25(units)={[x for x in r25 if gcd(x,25)==1]} W={list(W)}")

if __name__ == "__main__":
    part_A()
    part_B()
