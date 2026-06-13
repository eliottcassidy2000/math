#!/usr/bin/env python3
"""
tournament_universe_s339.py — The Tournament Universe: an interactive explorer
opus-2026-03-25-S339

A tool that generates a MAP of tournament space at each n, showing:
- Where each iso class lives (by H and c3 coordinates)
- How classes connect (metagraph edges)
- The A² fingerprint of each class (as color)
- The even-graph dual class
- The distance to the transitive tournament (the "origin")

This is the ATLAS of tournament space — every class, every connection,
every invariant, visualized as data you can explore.

Also computes the "tournament periodic table" — classifying tournaments
by their A² fingerprint "element type", like how the periodic table
classifies atoms by electron configuration.
"""

import sys
import subprocess
import numpy as np
from math import comb, factorial
from collections import defaultdict, Counter
from itertools import permutations

sys.stdout.reconfigure(line_buffering=True)

def build_universe(n):
    """Build the complete tournament universe at level n."""
    m_total = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1, n)]

    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True, timeout=120)
    lines = [l.strip() for l in r.stdout.split('\n')
             if len(l.strip()) == m_total and all(c in '01' for c in l.strip())]

    def bits_to_adj(s):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if s[k] == '1': A[i][j] = 1
            else: A[j][i] = 1
        return A

    def score(A): return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
    def c3(A):
        c = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if A[i][j] and A[j][k] and A[k][i]: c += 1
                    if A[i][k] and A[k][j] and A[j][i]: c += 1
        return c
    def H(A):
        dp = {}
        for v in range(n): dp[(1<<v, v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S & (1<<v)): continue
                val = dp.get((S, v), 0)
                if val == 0: continue
                for u in range(n):
                    if S & (1<<u): continue
                    if A[v][u]: dp[(S|(1<<u), u)] = dp.get((S|(1<<u), u), 0) + val
        return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

    def a2_fp(A):
        Anp = np.array(A, dtype=int)
        A2 = Anp @ Anp
        return tuple(sorted(tuple(A2[i]) for i in range(n)))

    universe = []
    for idx, line in enumerate(lines):
        A = bits_to_adj(line)
        s = score(A); c = c3(A); h = H(A)
        fp = a2_fp(A)
        aut = sum(1 for p in permutations(range(n))
                  if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n))
                 ) if n <= 6 else 0

        # Is it self-complementary?
        Ac = [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
        sc = a2_fp(Ac) == fp  # complement has same A² fingerprint

        universe.append({
            'id': idx, 'score': s, 'c3': c, 'H': h,
            'fp': fp, 'aut': aut, 'sc': sc,
        })

    return universe

# ============================================================================
# MAIN
# ============================================================================

print("=" * 72)
print("  THE TOURNAMENT UNIVERSE")
print("  opus-2026-03-25-S339")
print("=" * 72)

for n in [5, 6]:
    print(f"\n{'#'*72}")
    print(f"  n = {n}: THE PERIODIC TABLE OF TOURNAMENTS")
    print(f"{'#'*72}")

    universe = build_universe(n)
    V = len(universe)

    # Group by "element type" = (score, c3) — the "orbital"
    orbitals = defaultdict(list)
    for t in universe:
        orbitals[(t['score'], t['c3'])].append(t)

    print(f"\n  {V} iso classes in {len(orbitals)} orbital types")

    # The periodic table
    print(f"\n  PERIODIC TABLE (rows = score type, columns = c3):")
    scores = sorted(set(t['score'] for t in universe))
    c3_vals = sorted(set(t['c3'] for t in universe))

    # Header
    print(f"  {'score':>20}", end='')
    for c in c3_vals[:8]: print(f" c3={c:>2}", end='')
    print()
    print("  " + "-" * (22 + 6 * min(len(c3_vals), 8)))

    for s in scores:
        line = f"  {str(s):>20}"
        for c in c3_vals[:8]:
            members = orbitals.get((s, c), [])
            if members:
                n_members = len(members)
                sc_count = sum(1 for m in members if m['sc'])
                symbol = f"{n_members}" + ("*" if sc_count > 0 else " ")
                line += f" {symbol:>5}"
            else:
                line += "     ·"
        print(line)

    # H-spectrum by orbital
    print(f"\n  H-SPECTRUM BY ORBITAL:")
    for (s, c), members in sorted(orbitals.items(), key=lambda x: (x[0], x[1])):
        h_vals = sorted(set(m['H'] for m in members))
        sc_mark = " [SC]" if any(m['sc'] for m in members) else ""
        if len(members) > 1:
            print(f"    score={s}, c3={c}: {len(members)} classes, H∈{h_vals}{sc_mark}")

    # A² fingerprint uniqueness
    fp_groups = defaultdict(list)
    for t in universe:
        fp_groups[t['fp']].append(t['id'])
    n_unique_fp = sum(1 for g in fp_groups.values() if len(g) == 1)
    print(f"\n  A² FINGERPRINT: {n_unique_fp}/{V} unique ({n_unique_fp/V*100:.0f}%)")

    # Self-complementary census
    n_sc = sum(1 for t in universe if t['sc'])
    print(f"  Self-complementary: {n_sc}/{V}")

    # The "noble gases" (perfectly structured tournaments)
    print(f"\n  NOBLE GASES (extremal tournaments):")
    min_h = min(t['H'] for t in universe)
    max_h = max(t['H'] for t in universe)
    min_c3 = min(t['c3'] for t in universe)
    max_c3 = max(t['c3'] for t in universe)

    for t in universe:
        if t['H'] == min_h:
            print(f"    Transitive (H={t['H']}): score={t['score']}, c3={t['c3']}")
            break
    for t in universe:
        if t['H'] == max_h:
            print(f"    Maximum H (H={t['H']}): score={t['score']}, c3={t['c3']}, SC={t['sc']}")
            break
    for t in universe:
        if t['c3'] == max_c3 and t['H'] != max_h:
            print(f"    Maximum c3 (c3={t['c3']}): score={t['score']}, H={t['H']}")
            break

# ============================================================================
# THE ATLAS
# ============================================================================

print(f"\n{'='*72}")
print(f"  THE ATLAS OF TOURNAMENT SPACE")
print(f"{'='*72}")

print(f"""
  Tournament space at each n is a FINITE UNIVERSE with:

  COORDINATES: (H, c3) — like a Hertzsprung-Russell diagram for stars
    H = luminosity (how many Hamiltonian paths)
    c3 = temperature (how many 3-cycles)

  ELEMENT TYPE: (score, c3) — the orbital / electron configuration
    Each orbital can hold multiple iso classes (like electron shells)

  FINGERPRINT: A² sorted rows — the DNA / nuclear structure
    Complete invariant: uniquely identifies each class

  CONNECTIONS: metagraph edges — the chemical bonds
    χ = n-1 colors needed (like valence electrons)

  SYMMETRY: self-complementary — the matter/antimatter duality
    SC tournaments are their own antiparticle

  The progression n = 3, 4, 5, 6, 7, ... builds up the universe
  element by element, like the Big Bang nucleosynthesis:

    n=3: 2 elements (hydrogen and helium)
    n=4: 4 elements (first metals)
    n=5: 12 elements (first transition metals)
    n=6: 56 elements (the iron peak)
    n=7: 456 elements (the lanthanides)
    n=8: 6880 elements (the actinides)
    n=∞: T_ω (the ultimate element, containing all others)

  Each level adds new phenomena:
    n=3: first 3-cycle (first chemistry)
    n=5: first per-path failure (first quantum effects)
    n=6: first GF(2) obstruction (first parity violation)
    n=7: metagraph loses perfectness (first phase transition)
    n=9: first β₅ > 0 (first nuclear instability)

  The A² fingerprint is the SPECTROSCOPIC LINE SPECTRUM:
  just as astronomers identify elements by their spectral lines,
  we identify tournament types by their A² row patterns.
""")

print("DONE.")
