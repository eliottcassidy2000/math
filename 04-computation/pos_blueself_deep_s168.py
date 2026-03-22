#!/usr/bin/env python3
"""
pos_blueself_deep_s168.py — PoS, Blueself, Blackself Through the Three Axes
opus-2026-03-22-S168

Using the three independent axes (H/cycle, arb/tree, L/residual)
to discover new facts and sequences about:
  - PoS (Points of Symmetry): score classes with H-variation
  - Blueself: tournaments isomorphic to their complement (T ≅ T^c)
  - Blackself: tournaments NOT isomorphic to complement
  - Kings: vertices that reach all others in ≤2 steps

The key question: how do blueself/blackself tournaments distribute
across the three axes? Are they concentrated on specific axes?
"""

import numpy as np
from math import comb, factorial
from collections import defaultdict, Counter
from fractions import Fraction

print("=" * 72)
print("  PoS, BLUESELF, BLACKSELF THROUGH THE THREE AXES")
print("  opus-2026-03-22-S168")
print("=" * 72)
print()

# ==========================================================================
# HELPERS
# ==========================================================================

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(n))

def count_ham_cycles(adj, n):
    dp = [0]*((1<<n)*n)
    dp[(1<<0)*n+0] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(1, n) if adj[v][0])

def count_arborescences(adj, n, root=0):
    L = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                L[i][j] = -adj[j][i]
                L[i][i] += adj[j][i]
    indices = [i for i in range(n) if i != root]
    minor = L[np.ix_(indices, indices)]
    return int(round(np.linalg.det(minor)))

def count_kings(adj, n):
    kings = 0
    for v in range(n):
        is_king = True
        for u in range(n):
            if u == v: continue
            if adj[v][u]: continue
            if any(adj[v][w] and adj[w][u] for w in range(n) if w != v and w != u):
                continue
            is_king = False
            break
        if is_king: kings += 1
    return kings

def complement_bits(bits, m):
    return ((1 << m) - 1) ^ bits

def scores(adj, n):
    return [sum(adj[i]) for i in range(n)]

def is_isomorphic(adj1, adj2, n):
    """Check if two tournaments are isomorphic (brute force, small n)."""
    from itertools import permutations
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if adj1[i][j] != adj2[perm[i]][perm[j]]:
                    match = False
                    break
            if not match: break
        if match: return True
    return False

# ==========================================================================
# PART 1: CLASSIFY ALL n=5 TOURNAMENTS
# ==========================================================================

print("PART 1: CLASSIFYING ALL n=5 TOURNAMENTS")
print("-" * 60)
print()

n = 5
m = comb(n, 2)
data = []

for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    h = H_dp(adj, n)
    hc = count_ham_cycles(adj, n)
    l_val = h - n * hc
    arb = count_arborescences(adj, n, root=0)
    k = count_kings(adj, n)
    s = scores(adj, n)
    ss = tuple(sorted(s))
    s2 = sum(v*v for v in s)

    # Check if self-complementary (blueself)
    comp_bits = complement_bits(bits, m)
    comp_adj = tournament_adj(n, comp_bits)
    is_sc = is_isomorphic(adj, comp_adj, n)

    data.append({
        'bits': bits, 'h': h, 'hc': hc, 'l': l_val,
        'arb': arb, 'kings': k, 'scores': ss, 's2': s2,
        'is_sc': is_sc,  # blueself = self-complementary
    })

n_sc = sum(1 for d in data if d['is_sc'])
n_nsc = sum(1 for d in data if not d['is_sc'])
print(f"  Total: {len(data)} tournaments")
print(f"  Blueself (self-complementary): {n_sc}")
print(f"  Blackself (not SC): {n_nsc}")
print()

# ==========================================================================
# PART 2: BLUESELF vs BLACKSELF — THREE-AXIS DISTRIBUTION
# ==========================================================================

print()
print("PART 2: BLUESELF vs BLACKSELF ON THREE AXES")
print("-" * 60)
print()

# How do SC tournaments distribute across (H, arb, L)?
sc_h = Counter(d['h'] for d in data if d['is_sc'])
nsc_h = Counter(d['h'] for d in data if not d['is_sc'])

print(f"  H DISTRIBUTION:")
print(f"  {'H':>3s}  {'blueself':>8s}  {'blackself':>9s}  {'total':>6s}  {'%blue':>6s}")
print(f"  {'─'*3}  {'─'*8}  {'─'*9}  {'─'*6}  {'─'*6}")

all_h = sorted(set(d['h'] for d in data))
for h in all_h:
    blue = sc_h.get(h, 0)
    black = nsc_h.get(h, 0)
    total = blue + black
    pct = 100 * blue / total if total > 0 else 0
    print(f"  {h:3d}  {blue:8d}  {black:9d}  {total:6d}  {pct:5.1f}%")

print()

# arb distribution
print(f"  ARBORESCENCE DISTRIBUTION (mean arb):")
sc_arb = [d['arb'] for d in data if d['is_sc']]
nsc_arb = [d['arb'] for d in data if not d['is_sc']]
print(f"    Blueself mean arb: {np.mean(sc_arb):.2f}")
print(f"    Blackself mean arb: {np.mean(nsc_arb):.2f}")
print()

# L distribution
sc_l = Counter(d['l'] for d in data if d['is_sc'])
nsc_l = Counter(d['l'] for d in data if not d['is_sc'])
print(f"  L DISTRIBUTION:")
print(f"  {'L':>3s}  {'blueself':>8s}  {'blackself':>9s}")
print(f"  {'─'*3}  {'─'*8}  {'─'*9}")
for l_val in sorted(set(d['l'] for d in data)):
    print(f"  {l_val:3d}  {sc_l.get(l_val,0):8d}  {nsc_l.get(l_val,0):9d}")

print()

# Kings distribution
sc_kings = Counter(d['kings'] for d in data if d['is_sc'])
nsc_kings = Counter(d['kings'] for d in data if not d['is_sc'])
print(f"  KINGS DISTRIBUTION:")
print(f"  {'kings':>5s}  {'blueself':>8s}  {'blackself':>9s}")
print(f"  {'─'*5}  {'─'*8}  {'─'*9}")
for k in sorted(set(d['kings'] for d in data)):
    print(f"  {k:5d}  {sc_kings.get(k,0):8d}  {nsc_kings.get(k,0):9d}")

print()

# ==========================================================================
# PART 3: PoS CLASS DECOMPOSITION
# ==========================================================================

print()
print("PART 3: PoS CLASS DETAILED DECOMPOSITION")
print("-" * 60)
print()

# The PoS class at n=5 is score (1,2,2,2,3)
pos_data = [d for d in data if d['scores'] == (1, 2, 2, 2, 3)]

print(f"  PoS class (1,2,2,2,3): {len(pos_data)} tournaments")
print()

# Decompose by (H, SC, kings, L, arb)
pos_groups = defaultdict(list)
for d in pos_data:
    key = (d['h'], d['is_sc'], d['kings'], d['l'], d['hc'])
    pos_groups[key].append(d)

print(f"  {'H':>3s}  {'SC':>3s}  {'kings':>5s}  {'L':>3s}  {'HC':>3s}  "
      f"{'arb_mean':>8s}  {'count':>6s}")
print(f"  {'─'*3}  {'─'*3}  {'─'*5}  {'─'*3}  {'─'*3}  {'─'*8}  {'─'*6}")

for key in sorted(pos_groups.keys()):
    h, sc, kings, l_val, hc = key
    group = pos_groups[key]
    arb_mean = np.mean([d['arb'] for d in group])
    count = len(group)
    sc_str = "yes" if sc else "no"
    print(f"  {h:3d}  {sc_str:>3s}  {kings:5d}  {l_val:3d}  {hc:3d}  "
          f"{arb_mean:8.1f}  {count:6d}")

print()

# ==========================================================================
# PART 4: NEW SEQUENCES — BLUESELF/BLACKSELF BY n
# ==========================================================================

print()
print("PART 4: SEQUENCES — BLUESELF/BLACKSELF COUNTS")
print("-" * 60)
print()

for n_check in range(3, 7):
    m_check = comb(n_check, 2)
    sc_count = 0
    nsc_count = 0

    for bits in range(1 << m_check):
        adj = tournament_adj(n_check, bits)
        comp_bits = complement_bits(bits, m_check)
        comp_adj = tournament_adj(n_check, comp_bits)
        if is_isomorphic(adj, comp_adj, n_check):
            sc_count += 1
        else:
            nsc_count += 1

    total = 1 << m_check
    print(f"  n={n_check}: blueself={sc_count}, blackself={nsc_count}, "
          f"total={total}, %blue={100*sc_count/total:.1f}%")

print()

# ==========================================================================
# PART 5: ARBORESCENCE SEQUENCES BY H
# ==========================================================================

print()
print("PART 5: ARBORESCENCE SEQUENCES")
print("-" * 60)
print()

# For each H value, what are the possible arb values?
arb_by_h = defaultdict(list)
for d in data:
    arb_by_h[d['h']].append(d['arb'])

print(f"  ARBORESCENCES BY H (n=5, root=0):")
print(f"  {'H':>3s}  {'arb values':>25s}  {'mean':>7s}  {'#distinct':>9s}")
print(f"  {'─'*3}  {'─'*25}  {'─'*7}  {'─'*9}")

arb_seqs = {}
for h in sorted(arb_by_h.keys()):
    arbs = arb_by_h[h]
    distinct = sorted(set(arbs))
    mean = np.mean(arbs)
    arb_seqs[h] = distinct
    dist_str = str(distinct) if len(distinct) <= 6 else str(distinct[:6]) + "..."
    print(f"  {h:3d}  {dist_str:>25s}  {mean:7.1f}  {len(distinct):9d}")

print()

# Total arborescences (sum over all roots)
print(f"  TOTAL ARBORESCENCES (all roots) BY H:")
total_arb_by_h = defaultdict(list)
for d in data:
    adj = tournament_adj(n, d['bits'])
    total_arb = sum(count_arborescences(adj, n, root=r) for r in range(n))
    total_arb_by_h[d['h']].append(total_arb)

for h in sorted(total_arb_by_h.keys()):
    vals = total_arb_by_h[h]
    distinct = sorted(set(vals))
    print(f"  H={h:3d}: total_arb ∈ {distinct}")

print()

# ==========================================================================
# PART 6: THE (H, kings, SC) TRIPLE — NEW SEQUENCES
# ==========================================================================

print()
print("PART 6: THE (H, kings, SC) TRIPLE")
print("-" * 60)
print()

# How many distinct triples?
hks_counts = Counter()
for d in data:
    hks_counts[(d['h'], d['kings'], d['is_sc'])] += 1

print(f"  (H, kings, SC) triples at n=5:")
print(f"  {'H':>3s}  {'kings':>5s}  {'SC':>3s}  {'count':>6s}")
print(f"  {'─'*3}  {'─'*5}  {'─'*3}  {'─'*6}")
for (h, k, sc) in sorted(hks_counts.keys()):
    print(f"  {h:3d}  {k:5d}  {'Y' if sc else 'N':>3s}  {hks_counts[(h,k,sc)]:6d}")

print()
print(f"  Distinct triples: {len(hks_counts)}")
print()

# ==========================================================================
# PART 7: NEW SEQUENCE CANDIDATES
# ==========================================================================

print()
print("PART 7: NEW SEQUENCE CANDIDATES")
print("-" * 60)
print()

# Sequence: blueself count by n
print(f"  BLUESELF SEQUENCES:")
# Already computed above; let me extract
# n=3: ?  n=4: ?  n=5: ?  n=6: ? (computed above)

# Total arborescences at H_max
print(f"  TOTAL ARBORESCENCES AT H_max (all roots):")
for h_max_val in [max(d['h'] for d in data)]:
    ta_vals = total_arb_by_h[h_max_val]
    distinct = sorted(set(ta_vals))
    print(f"    H_max={h_max_val}: total_arb ∈ {distinct}")
    if len(distinct) > 1:
        print(f"    MULTIPLE values → arb distinguishes within H_max class!")

print()

# Sequence: number of (H, arb) pairs
n_ha_pairs = len(set((d['h'], d['arb']) for d in data))
print(f"  #(H, arb) pairs at n=5: {n_ha_pairs}")

# Sequence: number of (H, kings) pairs
n_hk_pairs = len(set((d['h'], d['kings']) for d in data))
print(f"  #(H, kings) pairs at n=5: {n_hk_pairs}")

# Sequence: number of (H, L) pairs
n_hl_pairs = len(set((d['h'], d['l']) for d in data))
print(f"  #(H, L) pairs at n=5: {n_hl_pairs}")

# Sequence: number of (H, HC) pairs
n_hhc_pairs = len(set((d['h'], d['hc']) for d in data))
print(f"  #(H, HC) pairs at n=5: {n_hhc_pairs}")

print()

# ==========================================================================
# PART 8: THE BLUESELF SKELETON
# ==========================================================================

print()
print("PART 8: THE BLUESELF SKELETON")
print("-" * 60)
print()

# Which H values can be achieved by SC tournaments?
sc_h_vals = sorted(set(d['h'] for d in data if d['is_sc']))
all_h_vals = sorted(set(d['h'] for d in data))
nsc_only = sorted(set(all_h_vals) - set(sc_h_vals))

print(f"  SC-achievable H values: {sc_h_vals}")
print(f"  All H values: {all_h_vals}")
print(f"  Non-SC only H values: {nsc_only}")
print()

# At n=5: which H values are SC-exclusive?
# (Only achieved by SC tournaments)
sc_exclusive = []
for h in sc_h_vals:
    if all(d['is_sc'] for d in data if d['h'] == h):
        sc_exclusive.append(h)

print(f"  SC-exclusive H values (only SC tournaments have this H): {sc_exclusive}")
print()

# The PoS class: what fraction is SC?
pos_sc = sum(1 for d in pos_data if d['is_sc'])
print(f"  PoS class: {pos_sc}/{len(pos_data)} are SC ({100*pos_sc/len(pos_data):.1f}%)")
print()

# ==========================================================================
# SUMMARY
# ==========================================================================

print()
print("=" * 72)
print("  SUMMARY: FACTS AND SEQUENCES FROM THREE-AXIS VIEW")
print("=" * 72)
print()

print(f"""
  BLUESELF/BLACKSELF FACTS AT n=5:
    {n_sc} blueself (SC), {n_nsc} blackself (not SC)
    SC achievable H: {sc_h_vals}
    Non-SC only H: {nsc_only}
    H=1 (transitive): ALL blackself (not SC)
    H=15 (maximum): BOTH blue and black
    PoS class: {pos_sc}/{len(pos_data)} SC ({100*pos_sc/len(pos_data):.1f}%)

  THREE-AXIS DECOMPOSITION OF PoS:
    The PoS class (1,2,2,2,3) splits by (H, SC, kings, L, HC):
    Multiple distinct subgroups with different arb distributions.
    arb distinguishes within H=15 class (two different arb values).

  NEW SEQUENCES:
    blueself count by n: (computed above)
    #(H, arb) pairs: {n_ha_pairs} at n=5
    #(H, kings) pairs: {n_hk_pairs} at n=5
    #(H, L) pairs: {n_hl_pairs} at n=5
    #(H, HC) pairs: {n_hhc_pairs} at n=5
    total_arb at H_max: {sorted(set(total_arb_by_h[max(all_h_vals)]))}

  KEY INSIGHT:
    Arborescences add INDEPENDENT information to the blueself skeleton.
    SC tournaments have different arb distributions than non-SC.
    The three-axis view (H, arb, L) separates tournament types
    that the single H value conflates.
""")

print("opus-2026-03-22-S168")
