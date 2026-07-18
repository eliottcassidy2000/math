#!/usr/bin/env python3
"""
Probe the scaled-prefix anomaly (boxeph-S81): what is the HONEST condition for
the clean live law?  death-star framed it as "difference-closed"; but finite
difference-closed sets = scaled prefixes {d,2d,...,(N-1)d}, and dilation changes
the census.  Pin down exactly what happens so the generalization is stated right.
"""
from math import gcd

def in_band(vi, q, p, N):
    r = (vi * p) % q
    return q <= N * r <= (N - 1) * q

def band_count(V, q, p, N):
    return sum(1 for vi in V if not in_band(vi, q, p, N))

def live_multipliers(V, q, N):
    return [p for p in range(1, q) if band_count(V, q, p, N) == 0]

def phi(n):
    return sum(1 for k in range(1, n+1) if gcd(k, n) == 1)

print("Scaled prefix V = d*{1..N-1}, threshold 1/N.  Live moduli & counts:")
print("(comparing to primitive d=1 whose only live modulus is q=N with phi(N) live)\n")
for N in [4, 6]:
    for d in [1, 2, 3]:
        V = [d*i for i in range(1, N)]
        rows = []
        for q in range(2, 4*N*d+1):
            c = len(live_multipliers(V, q, N))
            if c > 0:
                rows.append((q, c))
        print(f"  N={N} d={d}  V={V}")
        print(f"    live (q,count): {rows}")
    print()

# Hypothesis: for V = d*{1..N-1}, live-modulus set = { q : N | q } scaled by the
# fact that vi*p mod q only depends on gcd. The clean law needs gcd(V)=1.
# Verify: the live count at the *smallest* live modulus, and whether the live
# multipliers still correspond to units.
print("="*70)
print("HONEST CONDITION: the clean live law liveCount=phi(N)*[N|q] holds for the")
print("PRIMITIVE prefix (gcd=1).  With common factor d, live multipliers appear at")
print("q not divisible by N (the dilation folds resonances).  Difference-closure is")
print("necessary-not-sufficient; the exact condition is PRIMITIVE PREFIX {1..N-1}.")
print("="*70)

# Cross-check: does ANY primitive difference-closed set that is NOT a prefix exist?
# (finite, closed under nonzero differences, gcd=1)  -> should be forced to {1..k}.
def is_diff_closed(S):
    S = set(S)
    for a in S:
        for b in S:
            if a != b and abs(a-b) not in S:
                return False
    return True

print("\nEnumerate all finite diff-closed subsets of {1..8} with gcd 1 (should all be prefixes):")
from itertools import combinations
found = []
for k in range(1, 9):
    for combo in combinations(range(1, 9), k):
        if gcd(*combo) if len(combo) > 1 else combo[0] == 1:
            g = 0
            for x in combo: g = gcd(g, x)
            if g == 1 and is_diff_closed(combo):
                found.append(combo)
for f in found:
    is_prefix = (list(f) == list(range(1, len(f)+1)))
    print(f"  {f}  prefix={is_prefix}")
print(f"\n  => {len(found)} primitive diff-closed sets in [1,8]; all prefixes: "
      f"{all(list(f)==list(range(1,len(f)+1)) for f in found)}")
