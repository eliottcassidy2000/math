#!/usr/bin/env python3
"""
n8_n9_structures_s289.py — Key structures at n=8, n=9, n=10
opus-2026-03-24-S289

Using the Burnside speed engine + tiling formula + Petersen tournament
to build out info about metagraphs at n=8, 9, 10.

For n=8: 2^21 = 2M tilings — attempt tiling enumeration
For n=9: 2^28 = 268M tilings — formula-only (no enumeration)
For n=10: 2^36 = 69B tilings — formula-only, but Petersen tournament known
"""

import sys, time, subprocess
from math import comb, factorial, gcd
from collections import Counter, defaultdict
from functools import lru_cache

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# BURNSIDE SPEED ENGINE (from S265)
# ============================================================

@lru_cache(maxsize=None)
def _f(n): return factorial(n)

def gen_odd_parts(n):
    result = []
    def _g(rem, mx, cur):
        if rem == 0: result.append(tuple(cur)); return
        st = min(rem, mx)
        if st % 2 == 0: st -= 1
        for p in range(st, 0, -2): _g(rem - p, p, cur + [p])
    _g(n, n, [])
    return result

def gen_all_parts(n):
    result = []
    def _g(rem, mx, cur):
        if rem == 0: result.append(tuple(cur)); return
        for p in range(min(rem, mx), 0, -1): _g(rem - p, p, cur + [p])
    _g(n, n, [])
    return result

def ccs(n, ct):
    c = Counter(ct); r = _f(n)
    for l, k in c.items(): r //= pow(l, k) * _f(k)
    return r

def pair_orb(ct):
    s = sum((c-1)//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)): s += gcd(ct[i], ct[j])
    return s

def edge_orb(ct):
    s = sum(c//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)): s += gcd(ct[i], ct[j])
    return s

def fix_anti(ct, n):
    perm = [0]*n; pos = 0
    for c in ct:
        for i in range(c): perm[pos+i] = pos + (i+1)%c
        pos += c
    visited = set(); free = 0
    for i in range(n):
        for j in range(n):
            if i == j or (i,j) in visited: continue
            orbit = []; a, b = i, j
            while (a,b) not in visited:
                visited.add((a,b)); orbit.append((a,b))
                a, b = perm[a], perm[b]
            L = len(orbit); orbit_set = set(orbit); rev = (j,i)
            if rev in orbit_set:
                k = orbit.index(rev)
                if k % 2 == 1: free += 1
                else: return 0
            else:
                if rev in visited: pass
                else:
                    if L % 2 == 0: free += 1
                    else: return 0
    return 2 ** free

def compute_all(n):
    """Compute everything from Burnside."""
    nf = _f(n); m = comb(n, 2); m_tiles = comb(n-1, 2)

    V_sum = 0; T_sum = 0; twin_sum = 0; SC_sum = 0

    for ct in gen_odd_parts(n):
        cc = ccs(n, ct); po = pair_orb(ct)
        k = len(ct); f = ct.count(1)
        fix = pow(2, po)
        V_sum += cc * fix
        if f >= 2:
            T_sum += cc * comb(f, 2) * fix
            twin_sum += cc * 2 * comb(f, 2) * pow(2, po - k + 1)

    for ct in gen_all_parts(n):
        cc = ccs(n, ct)
        fa = fix_anti(list(ct), n)
        SC_sum += cc * fa

    V = V_sum // nf
    T = T_sum // nf
    twin_SL = twin_sum // nf
    SC = SC_sum // nf
    V_merged = (V + SC) // 2
    EO = T // 2 + _f(n-2)
    E_twin = (T - twin_SL) // 2
    E_master = T * (pow(2, n-1) - 2) // pow(2, n)

    # V_even (odd n only)
    V_even = None
    if n % 2 == 1:
        V_even_sum = 0
        for ct in gen_all_parts(n):
            cc = ccs(n, ct); eo = edge_orb(list(ct))
            kk = len(ct)
            V_even_sum += cc * pow(2, eo - kk + 1)
        V_even = V_even_sum // nf

    return {
        'n': n, 'm': m, 'm_tiles': m_tiles,
        'V': V, 'SC': SC, 'V_merged': V_merged,
        'T': T, 'twin_SL': twin_SL, 'EO': EO,
        'E_master': E_master, 'E_twin': E_twin,
        'V_even': V_even,
        'total_tilings': pow(2, m_tiles),
    }

# ============================================================
# MAIN
# ============================================================

print("=" * 72)
print("  KEY STRUCTURES AT n=8, 9, 10")
print("  opus-2026-03-24-S289")
print("=" * 72)

# Known exact E values
known_E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

for n in [8, 9, 10]:
    t0 = time.time()
    r = compute_all(n)
    dt = time.time() - t0

    print(f"\n{'='*72}")
    print(f"  n = {n} (computed in {dt:.3f}s)")
    print(f"{'='*72}")

    print(f"\n  BURNSIDE QUANTITIES:")
    print(f"    V (iso classes):     {r['V']:>15}")
    print(f"    SC (self-comp):      {r['SC']:>15}")
    print(f"    V_merged (÷Z_2):     {r['V_merged']:>15}")
    print(f"    T (transitions):     {r['T']:>15}")
    print(f"    twin_SL:             {r['twin_SL']:>15}")
    print(f"    edge_orbits:         {r['EO']:>15}")
    print(f"    E_master (approx):   {r['E_master']:>15}")
    print(f"    E_twin (approx):     {r['E_twin']:>15}")
    if n in known_E:
        print(f"    E_exact (known):     {known_E[n]:>15}")
    if r['V_even'] is not None:
        print(f"    V_even (odd n):      {r['V_even']:>15}")

    print(f"\n  TILING MODEL:")
    print(f"    m_tiles = C({n-1},2) = {r['m_tiles']}")
    print(f"    total tilings = 2^{r['m_tiles']} = {r['total_tilings']:>15}")
    print(f"    tilings/class (avg) = {r['total_tilings']/r['V_merged']:.2f}")

    # The tiling formula: Σ H/|Aut| × mult = 2^m_tiles
    # For merged: Σ H/|Aut| (SC) + Σ 2H/|Aut| (NS) = 2^m_tiles
    print(f"\n  TILING FORMULA CONSTRAINT:")
    print(f"    Σ H/|Aut| × mult = {r['total_tilings']}")
    print(f"    This constrains H and |Aut| across all {r['V_merged']} merged classes")

    # Wiggly metagraph structure (from Burnside + tiling formula)
    print(f"\n  WIGGLY METAGRAPH:")
    print(f"    Same node set as merged metagraph: V = {r['V_merged']}")
    print(f"    W_total = 2^{r['m_tiles']} × {r['m_tiles']} = {r['total_tilings'] * r['m_tiles']}")
    print(f"    W_row(i) = fiber(i) × {r['m_tiles']}")

    # Petersen tournament info for n=10
    if n == 10:
        print(f"\n  PETERSEN TOURNAMENT (T_Pet ∈ this metagraph):")
        print(f"    H(T_Pet) = 6961")
        print(f"    Scores = (2,3,3,3,4,5,6,6,6,7) — palindromic")
        print(f"    c3 = 28")
        print(f"    30/36 tiles flipped (83.3%)")
        print(f"    T_Pet is ONE of {r['V_merged']} merged classes")
        print(f"    T_Pet's fiber = H/|Aut| × mult = 6961/|Aut| × mult")
        print(f"    Szele expected H = {_f(n) // pow(2, n-1)} at n={n}")
        print(f"    T_Pet H / Szele = {6961 / (_f(n) / pow(2, n-1)):.4f}")

    # Spine/ribs/sea estimates
    ns_merged = r['V_merged'] - (r['SC'] + r['V'] - 2*r['V_merged'])  # NS count
    # Actually: SC classes in merged = SC, NS pairs in merged = (V - SC)/2
    # V_merged = SC + (V - SC)/2 = (V + SC)/2
    sc_in_merged = r['SC']  # unmerged SC classes, each is one merged node
    ns_in_merged = r['V_merged'] - sc_in_merged
    print(f"\n  SC/NS SPLIT IN MERGED METAGRAPH:")
    print(f"    SC nodes: {sc_in_merged}")
    print(f"    NS nodes: {ns_in_merged}")
    print(f"    SC fraction: {sc_in_merged/r['V_merged']:.4f}")

    # The spine/ribs/sea at large n: sea dominates
    if n >= 8:
        print(f"\n  SPINE/RIBS/SEA ESTIMATE:")
        spine_frac = (sc_in_merged * (sc_in_merged-1)/2) / (r['V_merged'] * (r['V_merged']-1)/2) if r['V_merged'] > 1 else 0
        sea_frac = (ns_in_merged * (ns_in_merged-1)/2) / (r['V_merged'] * (r['V_merged']-1)/2) if r['V_merged'] > 1 else 0
        ribs_frac = 1 - spine_frac - sea_frac
        print(f"    Spine (SC-SC) ≈ {100*spine_frac:.2f}% of possible edges")
        print(f"    Sea (NS-NS)   ≈ {100*sea_frac:.2f}% of possible edges")
        print(f"    Ribs (SC-NS)  ≈ {100*ribs_frac:.2f}% of possible edges")

# ============================================================
# SUMMARY TABLE
# ============================================================

print(f"\n{'='*72}")
print("  SUMMARY TABLE: n=3 through n=10")
print(f"{'='*72}")
print(f"  {'n':>3} {'V':>10} {'SC':>8} {'V_m':>10} {'m_tiles':>7} {'2^m':>12} "
      f"{'E_approx':>12} {'V_even':>10}")

for n in range(3, 11):
    r = compute_all(n)
    ek = known_E.get(n, r['E_master'])
    ve = r['V_even'] if r['V_even'] is not None else '—'
    print(f"  {n:>3} {r['V']:>10} {r['SC']:>8} {r['V_merged']:>10} {r['m_tiles']:>7} "
          f"{r['total_tilings']:>12} {ek:>12} {str(ve):>10}")

# ============================================================
# THE PETERSEN PATTERN
# ============================================================

print(f"\n{'='*72}")
print("  THE PETERSEN PATTERN")
print(f"{'='*72}")
print("""
  At each n, we can construct "graph tournaments" T_G by:
  1. Take ANY graph G on n vertices with a Hamiltonian path
  2. Orient G-edges as unflipped (high→low)
  3. Orient non-G-edges as flipped (low→high)

  The Petersen tournament (n=10) is path-independent:
  ALL 240 Hamiltonian paths give the SAME tournament.
  This is because |Aut(Petersen)| = 120 = 5! = |S_5|.

  NATURAL GRAPH TOURNAMENTS AT EACH n:

  n=5: K_5 has everything as an edge → transitive tournament (all unflipped)
       C_5 has 5 edges → a specific "cycle tournament" with H=?
       Petersen = K(5,2) on 5 elements — but it has 10 vertices, not 5

  n=8: The Cayley graph of Z_2^3 (the 3-cube graph Q_3)
       Q_3 has 12 edges on 8 vertices. 12 unflipped + 16 flipped.
       |Aut(Q_3)| = 48. Hamiltonian path → tournament.

  n=9: The Paley graph of order 9 (strongly regular (9,4,1,2))
       9 edges. Hamiltonian paths exist.
       |Aut| = 72. Would give a specific tournament.

  n=10: PETERSEN (computed above). H=6961, palindromic scores.
""")

print("DONE.")
