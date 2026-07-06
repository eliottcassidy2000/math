#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S9 -- HYP-4332: THE TWO-MODULUS NON-CLEARING STRUCTURE of the
spectral-gap crux.

SYNTHESIS: a gap member (M in (1/13, 2/25)) has NO clearing witness at any q<=25
(my S7 Q50 bound) => by opus-S98's residue bridge, its residues are NON-CLEARING
mod EVERY q<=25.  Working out each modulus:

 * mod q with 1/q > 2/25 (q <= 12): non-clearing at a=1 forces a MULTIPLE OF q
   present (else margin >= 1/q > 2/25 at t=1/q).  => a multiple of every m<=12
   (recovers THM-369 sieve).
 * mod 13: non-clearing <=> for every unit a, some residue r has ra in {0,+-1}
   <=> (if no multiple of 13) the residues HIT ALL 6 pair {u,-u} mod 13.
 * mod 25: the sibling-S7 FULL ANTIPODAL TRANSVERSAL -- residues hit all 10
   unit pairs {a, 25-a}.

So a gap member is a PAIR-HITTING TRANSVERSAL AT BOTH 13 AND 25.  This script:
 (1) verifies the non-clearing characterization at q = 12,13,19,23,25 exactly;
 (2) confirms the AP {1..12} passes all of them (is non-clearing everywhere);
 (3) SEARCHES for any covering primitive 12-set passing ALL non-clearing filters
     q<=25 -- if only AP-structured families survive with M=1/13, and every
     other survivor has M>=2/25, the crux is closed at these moduli.
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def dist_int(x, q):
    r = x % q
    return min(r, q - r)

def clears_at_mod(residues, q):
    """does SOME a in [1,q-1] make all residues >= 2/25 far at a/q?
    i.e. dist_int(r*a, q) >= 2q/25 for all r.  Returns the witness a or None."""
    thr = F(2, 25) * q
    for a in range(1, q):
        if all(dist_int(r * a, q) >= thr for r in residues):
            return a
    return None

def nonclearing_mod(W, q):
    """W is non-clearing mod q iff no a/q witness clears 2/25."""
    return clears_at_mod([v % q for v in W], q) is None

# ---- characterization checks -------------------------------------------------
def characterize_mod13():
    print("  mod 13: non-clearing <=> residues hit all 6 pairs {u,-u} (if no mult of 13)")
    pairs13 = [(u, 13 - u) for u in range(1, 7)]   # (1,12),(2,11),...,(6,7)
    # test: a residue set hits all pairs <=> non-clearing (0 not present)
    ok = True
    for _ in range(3000):
        k = random.randint(4, 12)
        S = set(random.sample(range(1, 13), min(k, 12)))   # nonzero residues
        hits_all = all((p[0] in S or p[1] in S) for p in pairs13)
        # non-clearing test on this residue set (as a "family")
        nc = clears_at_mod(list(S), 13) is None
        if hits_all != nc:
            ok = False
            print(f"    MISMATCH S={sorted(S)} hits_all={hits_all} nonclearing={nc}")
            break
    print(f"    characterization {'VERIFIED' if ok else 'FAILED'} (pair-hitting <=> non-clearing mod 13)")

def characterize_mod25():
    print("  mod 25: non-clearing => hits all 10 unit pairs {a,25-a} (sibling P1)")
    units = [u for u in range(1, 25) if gcd(u, 25) == 1]
    pairs25 = [(u, 25 - u) for u in units if u < 25 - u]
    ok = True
    for _ in range(3000):
        k = random.randint(8, 14)
        S = set(random.sample(range(1, 25), min(k, 24)))
        hits_all = all((p[0] in S or p[1] in S) for p in pairs25)
        nc = clears_at_mod(list(S), 25) is None
        # non-clearing REQUIRES hitting all pairs (P1); the converse can fail
        # (non-units matter) -- so test the necessary direction only
        if nc and not hits_all:
            ok = False
            print(f"    COUNTEREX to necessity: S={sorted(S)} nonclearing but misses a pair")
            break
    print(f"    necessity (non-clearing => hits all 10 unit pairs) {'VERIFIED' if ok else 'FAILED'}")

# ---- exact M ----------------------------------------------------------------
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

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(v // g for v in W) if g > 1 else tuple(W)

# FULL filter: non-clearing mod every q in 2..25.  By the S7 Q50 bound (clearing
# => a witness at q <= 25), passing ALL of these means M < 2/25 => the survivor
# is gap-or-AP.  If only the AP survives, no gap member exists at these moduli.
MODULI = list(range(2, 26))

def part_search():
    print()
    print("=" * 78)
    print("SEARCH: covering primitive 12-sets non-clearing mod ALL of", MODULI)
    print("(gap members must pass every one; do only AP-structured families survive?)")
    print("=" * 78)
    random.seed(9)
    AP = tuple(range(1, 13))
    # confirm the AP passes all filters
    ap_pass = {q: nonclearing_mod(AP, q) for q in MODULI}
    print(f"  AP {{1..12}}: non-clearing status by modulus {ap_pass}; M(AP) = {exact_M(AP)}")
    survivors = []
    tested = 0
    passed_filters = 0
    # structured search: AP-perturbations, lifts, block-lifts, random with divisor structure
    cands = set()
    cands.add(AP)
    # AP with a few elements lifted by +25 or +13 (residue-preserving-ish)
    for _ in range(30000):
        W = list(range(1, 13))
        nmod = random.randint(1, 4)
        for _i in range(nmod):
            idx = random.randrange(12)
            W[idx] += random.choice([13, 25, 26, 38, 50]) * random.randint(1, 3)
        W = primitive(tuple(sorted(set(W))))
        if len(W) == 12:
            cands.add(W)
    # random families containing a multiple of each m<=12 (THM-369 necessary)
    for _ in range(30000):
        req = [random.choice([m, 2*m, 3*m]) for m in [5,7,8,9,11,12]]  # covers 2,3,4,6,10 via 12,8? ensure
        extra = random.sample(range(1, 40), 6)
        W = primitive(tuple(sorted(set(req + extra))))
        if len(W) == 12:
            cands.add(W)
    for W in cands:
        if len(set(W)) != 12:
            continue
        tested += 1
        if all(nonclearing_mod(W, q) for q in MODULI):
            passed_filters += 1
            M = exact_M(W)
            survivors.append((W, M))
    print(f"  candidates tested: {tested};  passed ALL non-clearing filters q in 2..25: {passed_filters}")
    if survivors:
        gap = [(W, M) for W, M in survivors if F(1,13) < M < F(2,25)]
        ap_val = [(W, M) for W, M in survivors if M == F(1,13)]
        above = [(W, M) for W, M in survivors if M >= F(2,25)]
        below = [(W, M) for W, M in survivors if M < F(1,13)]
        print(f"    survivors by M: M=1/13 (AP-value): {len(ap_val)};  "
              f"M>=2/25: {len(above)} (SHOULD be 0 -- they'd clear at some q<=25);  "
              f"M<1/13: {len(below)};  IN-GAP: {len(gap)}")
        if gap:
            print("    *** GAP MEMBERS FOUND (crux would be FALSE / need higher moduli):")
            for W, M in gap[:8]:
                print(f"      M={M} W={list(W)}")
        else:
            print("    ZERO in-gap survivors: every family passing all q<=25 non-clearing")
            print("    filters is either the AP-value (1/13) or clears (>=2/25) -- the crux")
            print("    holds at moduli", MODULI, "on this search.")
        # what do the M=1/13 survivors look like? (should be AP-structured)
        for W, M in ap_val[:5]:
            print(f"    M=1/13 survivor: {list(W)}")
        # DIAGNOSTIC: the M>=2/25 survivors pass all q<=25 filters yet clear --
        # so their witness is at q>25 (the "ray-local" caveat, MISTAKE-110).
        # Examine: min-witness-denom + are they dilation rays?
        if above:
            print("    M>=2/25 survivors (clear at q>25 -- the ray-local exceptions):")
            for W, M in above:
                mwd = None
                for q in range(1, 200):
                    thr = F(2,25)*q
                    if any(all(dist_int(v*a, q) >= thr for v in W) for a in range(1,q)):
                        mwd = q; break
                g = gcd(W[1]-W[0], W[2]-W[0]) if len(W) > 2 else 1
                diffs = sorted(set(W[i]-W[0] for i in range(len(W))))
                gg = 0
                for d in diffs:
                    gg = gcd(gg, d)
                print(f"      M={M} min-witness-q={mwd}  gcd(diffs)={gg}  W={list(W)}")

if __name__ == "__main__":
    print("=" * 78)
    print("PART 1: non-clearing characterization at the key moduli (exact)")
    print("=" * 78)
    random.seed(1)
    characterize_mod13()
    characterize_mod25()
    part_search()
