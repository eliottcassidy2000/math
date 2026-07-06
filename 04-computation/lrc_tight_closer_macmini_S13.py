#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S13 -- HYP-4392 (the closer test): do the TWO PROVEN
necessary conditions for tightness force the AP?

A tight (M=1/13) primitive 12-family satisfies BOTH (each a green corpus lemma):
 (RP) residue_pinning_13: residues {1..12} mod 13, bijectively;
 (SV) divisor-protection (sieve): a multiple of each m in 2..12
      (at t=1/m, margin <= M = 1/13 < 1/m forces some m | v_i).

QUESTION: does (RP) + (SV) + primitive + 12 distinct elements => the AP {1..12}
(uniquely)?  If YES, the tight-side rigidity CLOSES: tight => (RP)+(SV) => AP,
a finite combinatorial fact composed from two proven lemmas.  This searches for
any NON-AP family satisfying both.
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def dist_int(x, q):
    r = x % q
    return min(r, q - r)

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
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

def residue_pinned(W):
    """residues {1..12} mod 13, bijectively (no mult of 13)."""
    res = [v % 13 for v in W]
    return set(res) == set(range(1, 13)) and len(set(res)) == 12

def sieve_ok(W):
    """a multiple of each m in 2..12."""
    for m in range(2, 13):
        if not any(v % m == 0 for v in W):
            return False
    return True

FLOOR = F(1, 13)

def is_dilated_AP(W):
    """W == c*{1..12} for some c."""
    W = sorted(W)
    d = W[0]
    return all(W[i] == d * (i + 1) for i in range(12))

def main():
    print("=" * 78)
    print("THE TIGHT-SIDE CLOSER: does (RP residue-pinning) + (SV sieve) => AP?")
    print("=" * 78)
    AP = tuple(range(1, 13))
    print(f"  AP {{1..12}}: RP={residue_pinned(AP)}, SV={sieve_ok(AP)}, "
          f"primitive={primitive(AP)==AP}, M={exact_M(AP)}")
    random.seed(13)
    # SEARCH for non-dilated-AP families satisfying RP + SV + primitive.
    # Build by assigning each residue class mod 13 an element (residue + 13k),
    # then filter by the sieve.
    found = []
    found_M = {}
    tested = 0
    kmax = 12   # allow lifts up to residue + 13*12 (height ~ 160) -- generous
    for _ in range(400000):
        perm = list(range(1, 13)); random.shuffle(perm)
        W = tuple(sorted(perm[i] + 13 * random.randint(0, kmax) for i in range(12)))
        if len(set(W)) != 12:
            continue
        Wp = primitive(W)
        if len(set(Wp)) != 12 or not residue_pinned(Wp) or not sieve_ok(Wp):
            continue
        tested += 1
        if is_dilated_AP(Wp) or Wp == AP:
            continue
        # a NON-AP family satisfying both conditions
        found.append(Wp)
        M = exact_M(Wp)
        found_M[M] = found_M.get(M, 0) + 1
    print(f"  families satisfying RP + SV + primitive found: {tested}")
    print(f"  of those, NON-(dilated-)AP: {len(found)}")
    if found:
        print("    (RP)+(SV) do NOT force the AP -- non-AP families satisfy both:")
        # show whether they're tight or not
        tight = [W for W in found if exact_M(W) == FLOOR]
        print(f"    of the non-AP RP+SV families: TIGHT (M=1/13): {len(tight)}, non-tight: {len(found)-len(tight)}")
        for W in found[:8]:
            print(f"      M={exact_M(W)} W={list(W)}  (multiples: "
                  f"{[m for m in range(2,13) if any(v%m==0 for v in W)]})")
        if tight:
            print(f"    *** NON-AP TIGHT families satisfying RP+SV: {len(tight)} -- these would")
            print(f"        REFUTE 'primitive tight => AP' UNLESS excluded by another condition!")
            for W in tight[:6]:
                print(f"        TIGHT NON-AP: {list(W)}")
        else:
            print("    BUT none are tight (all M > 1/13) -- so (RP)+(SV) are necessary but")
            print("    NOT sufficient; tightness picks the AP among RP+SV families.  The")
            print("    tight-side closer needs a THIRD condition beyond RP+SV.")
    else:
        print("    (RP) + (SV) + primitive => AP UNIQUELY (on this search, height ~160)!")
        print("    => tight => (RP)+(SV) [both proven] => AP: the tight-side rigidity")
        print("    CLOSES as a finite combinatorial fact from two green lemmas.")

if __name__ == "__main__":
    main()
