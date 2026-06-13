#!/usr/bin/env python3
"""
neutral_arcs_s20dj.py — Classify neutral arcs by mechanism
kind-pasteur-2026-03-23-S20dj

A neutral arc (u,v) in tournament T means: flipping (u,v) gives T' ~ T.
The isomorphism tau with tau(T') = T can be:
  1. TWIN: tau = transposition (u,v). Works iff u,v have same out-set
     among other vertices. Simplest mechanism.
  2. COMPLEX: tau involves moving other vertices. Harder to characterize.

For each class C and each neutral arc orbit, classify the mechanism.
Then see if twin-neutral has a formula.

TWIN COUNT FORMULA:
  For a specific arc (u,v) in a random labeled tournament:
  P(twin) = 2^{-(n-2)} (each other vertex independently matches).
  Expected twin arcs per tournament = C(n,2) * 2^{-(n-2)}.
  If this accounts for ALL neutral arcs, then:
  SL_mine = T_n * twin_fraction (computable from Burnside!)
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  NEUTRAL ARC CLASSIFICATION: TWIN vs COMPLEX")
print("  kind-pasteur-2026-03-23-S20dj")
print("=" * 80)

for n in [3, 4, 5, 6]:
    t0 = time.time()
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    perms = list(permutations(range(n)))

    def make_adj(bits):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(ALL_ARCS):
            if bits & (1 << k): A[i][j] = 1
            else: A[j][i] = 1
        return A

    def canon(A):
        best = None
        for p in perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best

    def Hdp(A):
        dp = {}
        for v in range(n): dp[(1<<v, v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S & (1<<v)): continue
                val = dp.get((S,v), 0)
                if val == 0: continue
                for u in range(n):
                    if S & (1<<u): continue
                    if A[v][u]:
                        dp[(S|(1<<u), u)] = dp.get((S|(1<<u), u), 0) + val
        return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

    def is_twin(A, u, v):
        """Check if u,v are twins: same out-set among other vertices."""
        for w in range(n):
            if w == u or w == v: continue
            if A[u][w] != A[v][w]: return False
        return True

    # Build classes
    classes = {}
    for bits in range(1 << m):
        A = make_adj(bits)
        cn = canon(A)
        if cn not in classes:
            classes[cn] = {'rep_bits': bits, 'rep_adj': A, 'H': Hdp(A), 'members': []}
        classes[cn]['members'].append(bits)

    # For each class: count neutral arcs, decomposed by twin vs complex
    print(f"\nn={n}: {len(classes)} classes, {time.time()-t0:.1f}s build")

    total_twin_neutral = 0
    total_complex_neutral = 0
    total_neutral = 0

    print(f"  {'H':>3} {'|Aut|':>5} {'orb':>5} {'neut':>5} {'twin':>5} {'cplx':>5} {'base_n':>6} {'tile_n':>6}")

    for cn in sorted(classes.keys(), key=lambda c: classes[c]['H']):
        ci = classes[cn]
        A = ci['rep_adj']
        H = ci['H']
        bits = ci['rep_bits']
        orbit_size = len(ci['members'])
        aut_size = factorial(n) // orbit_size

        # Automorphisms
        auts = [p for p in perms if all(A[p[i]][p[j]] == A[i][j]
                for i in range(n) for j in range(n))]

        # Arc orbits under Aut
        arc_parent = list(range(m))
        def find(x):
            while arc_parent[x] != x: arc_parent[x] = arc_parent[arc_parent[x]]; x = arc_parent[x]
            return x
        def unite(x, y): arc_parent[find(x)] = find(y)
        for sigma in auts:
            for k, (i,j) in enumerate(ALL_ARCS):
                si, sj = sigma[i], sigma[j]
                unite(k, ALL_ARCS.index((min(si,sj), max(si,sj))))

        orbits = defaultdict(list)
        for k in range(m): orbits[find(k)].append(k)

        # Classify each orbit
        neutral_count = 0
        twin_count = 0
        complex_count = 0
        base_neutral = 0
        tile_neutral = 0

        for root, members in orbits.items():
            k = members[0]
            u, v = ALL_ARCS[k]
            A2 = make_adj(bits ^ (1 << k))
            cn2 = canon(A2)

            if cn2 == cn:  # neutral!
                neutral_count += 1
                is_base = (abs(u - v) == 1)

                if is_base:
                    base_neutral += 1
                else:
                    tile_neutral += 1

                # Check if twin mechanism works
                if is_twin(A, u, v):
                    twin_count += 1
                else:
                    complex_count += 1

        total_twin_neutral += twin_count
        total_complex_neutral += complex_count
        total_neutral += neutral_count

        if neutral_count > 0 or n <= 5:
            print(f"  {H:3d} {aut_size:5d} {len(orbits):5d} {neutral_count:5d} "
                  f"{twin_count:5d} {complex_count:5d} {base_neutral:6d} {tile_neutral:6d}")

    print(f"  TOTALS: neutral={total_neutral} twin={total_twin_neutral} "
          f"complex={total_complex_neutral}")
    print(f"  SL_mine = {total_neutral}")
    print(f"  Twin fraction of neutral: {total_twin_neutral/total_neutral:.4f}" if total_neutral > 0 else "")

    # Expected twin arcs per tournament
    expected_twins = comb(n,2) * 2**(-(n-2))
    print(f"  Expected twin arcs per random T: {expected_twins:.3f}")
    print(f"  Twin neutral orbits / V_n: {total_twin_neutral/len(classes):.3f}")
    print(f"  Complex neutral orbits / V_n: {total_complex_neutral/len(classes):.3f}")

print(f"\n{'='*60}")
print("SYNTHESIS")
print(f"{'='*60}")
print("""
The TWIN mechanism: arc (u,v) is neutral because swapping u,v
is an isomorphism (they have identical out-neighborhoods).

If ALL neutral arcs are twin-neutral, then:
  SL_mine = sum_C (# twin arc orbits in C)
  = (1/n!) * sum_T (#twin pairs in T) * C(n,2)/C(n,2)
  Hmm, this needs more careful counting.

The question: does TWIN account for ALL neutral arcs at large n?
Or are there always COMPLEX neutral arcs too?
""")

print("DONE.")
print("=" * 80)
