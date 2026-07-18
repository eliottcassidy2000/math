#!/usr/bin/env python3
"""
Farey resonance profile of the HARD covering families (boxeph-S81).

THM-998 recast the covering/non-covering dichotomy in Farey terms:
  - a family has a deep circle at denominator b iff it has speeds divisible by b;
  - NON-covering  = some small circle is EMPTY  => live witness at its center;
  - COVERING      = every small circle non-empty => no easy witness.

But "non-empty" is not "full".  The census/B5 route needs a modulus where enough
LIVE multipliers survive.  Here we compute, for the crux covering families, the
resonance profile: for each denominator b, how many speeds are divisible by b
(the "fill" of circle b), and where the family FIRST goes live (a rational
loneliness witness).  Goal: see whether an UNDER-FILLED small circle is the
census target the fleet's adaptive-modulus search keeps rediscovering.
"""
from math import gcd
from fractions import Fraction

N = 14
THRESH = Fraction(1, 14)

def band_count(V, q, p):
    c = 0
    for vi in V:
        r = (vi * p) % q
        if not (q <= 14 * r <= 13 * q):
            c += 1
    return c

def first_live(V, qmax=400):
    for q in range(2, qmax+1):
        for p in range(1, q):
            if band_count(V, q, p) == 0:
                return (q, p, Fraction(p, q))
    return None

def resonance_fill(V, bmax=28):
    """circle b 'fill' = #speeds divisible by b (these nest at 0 near a/b)."""
    return {b: sum(1 for v in V if v % b == 0) for b in range(1, bmax+1)}

def is_covering(V, upto=14):
    """covering = a multiple of every b in 2..upto present."""
    return all(any(v % b == 0 for v in V) for b in range(2, upto+1))

def continuous_M(V):
    """Exact M = max_t min_i ||v_i t|| via candidate rationals a/(v_i+v_j) and a/v_i
    (the min-of-sawtooths is piecewise linear; optima at these kinks/crossings)."""
    cands = set()
    vs = sorted(set(V))
    for v in vs:
        for a in range(1, v):
            cands.add(Fraction(a, v))
    for i in range(len(vs)):
        for j in range(len(vs)):
            s = vs[i] + vs[j]
            for a in range(1, s):
                cands.add(Fraction(a, s))
            d = abs(vs[i] - vs[j])
            if d:
                for a in range(1, d):
                    cands.add(Fraction(a, d))
    best = Fraction(0)
    for t in cands:
        m = min(min((v*t) % 1, 1-((v*t) % 1)) for v in V)
        if m > best:
            best = m
    return best

FAMILIES = {
    "tight prefix {1..13}     ": list(range(1,14)),
    "GW doubling {1..11,13,24}": list(range(1,12))+[13,24],
    "deep well {1..12,182}    ": list(range(1,13))+[182],
    "residue {1..11,13,84}    ": list(range(1,12))+[13,84],
    "sporadic {1..11,13,24}   ": list(range(1,12))+[13,24],
    "missing-2 {1,3,4,..,14}  ": [1]+list(range(3,15)),
}

print("="*78)
print("FAREY RESONANCE PROFILE of LRC(14) crux families")
print("="*78)
for name, V in FAMILIES.items():
    cov = is_covering(V)
    M = continuous_M(V)
    fl = first_live(V, 300)
    fill = resonance_fill(V, 16)
    empty_b = [b for b in range(2,15) if fill[b]==0]
    print(f"\n{name}")
    print(f"   covering={cov}   M={M} (={float(M):.5f})   first-live (q,p,t)={fl}")
    print(f"   empty small circles b in 2..14: {empty_b if empty_b else 'NONE (covering)'}")
    print(f"   circle fill (b:count divisible), b=1..14: "
          + " ".join(f"{b}:{fill[b]}" for b in range(1,15)))
print()
print("READING TO WATCH FOR:")
print(" - non-covering families have an empty small circle b0; first-live t = 1/b0.")
print(" - covering families: no empty small circle; first-live sits at a modulus")
print("   tied to the LEAST-FILLED circle (the census's natural target).")
