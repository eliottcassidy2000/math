#!/usr/bin/env python3
"""
hc_lens_s20v.py — kind-pasteur-2026-03-22-S20v

See the ENTIRE project through the H vs HC lens.

1. Forbidden H=7: what does HC say?
2. Score determination: is HC more or less score-determined than H?
3. Girth of Omega vs HC
4. The algebraic formula H = C_n - S_2: what about HC?
5. Paley tournaments: is Paley also the HC maximizer?
6. The source-sink embedding at ALL score classes
7. H mod n as a tournament invariant

Author: kind-pasteur-2026-03-22-S20v
"""
import sys, time
import numpy as np
from math import comb
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def count_hc(A, n):
    full = (1 << n) - 1
    dp = defaultdict(int)
    dp[(1, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[(full, v)] for v in range(n) if A[v][0])

def c3_from_scores(n, s):
    S2 = int(sum(s*s))
    return comb(n, 3) - (S2 - comb(n, 2)) // 2

def is_strongly_connected(A, n):
    # BFS from vertex 0
    visited = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for w in range(n):
            if A[v][w] and w not in visited:
                visited.add(w)
                queue.append(w)
    if len(visited) < n: return False
    # BFS on transpose
    visited = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for w in range(n):
            if A[w][v] and w not in visited:
                visited.add(w)
                queue.append(w)
    return len(visited) == n

print("="*60)
print("  SEEING EVERYTHING THROUGH THE HC LENS")
print("="*60)

# ================================================================
# EXHAUSTIVE AT n=5
# ================================================================
n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

all_data = []
for bits in range(2**m):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    s = A.sum(axis=1).astype(int)
    S2 = int(sum(s*s))
    H = count_hp(A, n)
    HC = count_hc(A, n)
    c3 = c3_from_scores(n, s)
    sc = is_strongly_connected(A, n)
    score_class = tuple(sorted(s))
    L = H - n * HC
    all_data.append({
        'bits': bits, 'H': H, 'HC': HC, 'L': L, 'c3': c3,
        'S2': S2, 'score': score_class, 'sc': sc
    })

# ================================================================
# 1. FORBIDDEN H=7 THROUGH HC LENS
# ================================================================
print(f"\n  1. FORBIDDEN H=7 THROUGH HC LENS")
print()

# H=7 would decompose as: 7 = 5*HC + L where L = 7 mod 5 = 2.
# So HC would be 1 and L = 2.
# OR: HC = 0 and L = 7 (but L can't exceed H, and 7 < 5*2=10, so HC could be 0 or 1).
# Actually: 7 = 5*1 + 2 (HC=1, L=2) OR 7 = 5*0 + 7 (HC=0, L=7).
# From our data: HC for various H values:
print("  H values and their HC distributions:")
H_to_HC = defaultdict(list)
for d in all_data:
    H_to_HC[d['H']].append(d['HC'])

for H in sorted(H_to_HC.keys()):
    hcs = sorted(set(H_to_HC[H]))
    counts = {hc: H_to_HC[H].count(hc) for hc in hcs}
    print(f"    H={H:>3d}: HC in {hcs}, L in {[H - n*hc for hc in hcs]}")

print()
print("  If H=7 existed: HC in {0, 1}. L in {7, 2}.")
print("  H=5 has HC=0 always (L=5). H=9 has HC=1 always (L=4).")
print("  H=7 (if it existed) would need HC=0 (L=7>5=H for n=5?)")
print("  Wait: L = H - n*HC. If HC=0: L=7. If HC=1: L=2.")
print("  L=7 means 7 linear (non-cyclic) HP. Possible in principle.")
print("  L=2 means 1 cycle (giving 5 HP) + 2 extra linear HP.")
print()
print("  THE HC CONSTRAINT ON H=7:")
print("  If HC=1: need 1 Ham cycle + 2 extra HP not in any cycle.")
print("  If HC=0: need 7 linear HP with no closeable cycle.")
print("  NEITHER is inherently impossible from the HC structure alone.")
print("  The impossibility must come from the CYCLE STRUCTURE (alpha_1=3 overshoot).")
print()

# ================================================================
# 2. IS HC MORE SCORE-DETERMINED THAN H?
# ================================================================
print("  2. SCORE DETERMINATION OF HC")
print()

score_to_H = defaultdict(set)
score_to_HC = defaultdict(set)
for d in all_data:
    score_to_H[d['score']].add(d['H'])
    score_to_HC[d['score']].add(d['HC'])

print("  Score class -> (H values, HC values):")
for score in sorted(score_to_H.keys()):
    Hs = sorted(score_to_H[score])
    HCs = sorted(score_to_HC[score])
    H_det = "unique" if len(Hs) == 1 else "AMBIGUOUS"
    HC_det = "unique" if len(HCs) == 1 else "AMBIGUOUS"
    if len(Hs) > 1 or len(HCs) > 1:
        print(f"    {list(score)}: H={Hs} ({H_det}), HC={HCs} ({HC_det})")

print()

# OCR for H and HC
H_arr = np.array([d['H'] for d in all_data], dtype=float)
HC_arr = np.array([d['HC'] for d in all_data], dtype=float)

# Score-conditioned means
score_H_mean = {}
score_HC_mean = {}
score_count = {}
for d in all_data:
    sc = d['score']
    if sc not in score_H_mean:
        score_H_mean[sc] = 0; score_HC_mean[sc] = 0; score_count[sc] = 0
    score_H_mean[sc] += d['H']
    score_HC_mean[sc] += d['HC']
    score_count[sc] += 1
for sc in score_H_mean:
    score_H_mean[sc] /= score_count[sc]
    score_HC_mean[sc] /= score_count[sc]

# Var(E[H|score]) / Var(H)
var_H = float(np.var(H_arr))
var_HC = float(np.var(HC_arr))

cond_H = np.array([score_H_mean[d['score']] for d in all_data])
cond_HC = np.array([score_HC_mean[d['score']] for d in all_data])

var_cond_H = float(np.var(cond_H))
var_cond_HC = float(np.var(cond_HC))

OCR_H = var_cond_H / var_H if var_H > 0 else 0
OCR_HC = var_cond_HC / var_HC if var_HC > 0 else 0

print(f"  OCR for H: {OCR_H:.6f} ({100*OCR_H:.2f}%)")
print(f"  OCR for HC: {OCR_HC:.6f} ({100*OCR_HC:.2f}%)")
print(f"  HC is {'MORE' if OCR_HC > OCR_H else 'LESS'} score-determined than H!")
print()

# ================================================================
# 3. STRONG CONNECTIVITY AND HC
# ================================================================
print("  3. STRONG CONNECTIVITY AND HC")
print()

sc_to_HC = defaultdict(list)
for d in all_data:
    sc_to_HC[d['sc']].append(d['HC'])

for sc_flag in [True, False]:
    hcs = sc_to_HC[sc_flag]
    avg_hc = sum(hcs)/len(hcs)
    max_hc = max(hcs)
    pct_zero = 100 * sum(1 for hc in hcs if hc == 0) / len(hcs)
    print(f"  SC={sc_flag}: count={len(hcs)}, avg HC={avg_hc:.2f}, max HC={max_hc}, HC=0: {pct_zero:.1f}%")

# Moon's theorem: every strongly connected tournament has a Hamiltonian cycle.
# So HC > 0 iff strongly connected?
hc_zero_sc = sum(1 for d in all_data if d['HC'] == 0 and d['sc'])
hc_pos_notsc = sum(1 for d in all_data if d['HC'] > 0 and not d['sc'])
print(f"\n  HC=0 AND strongly connected: {hc_zero_sc} (should be 0 by Moon)")
print(f"  HC>0 AND NOT strongly connected: {hc_pos_notsc} (should be 0)")
print(f"  Moon's theorem verified: HC>0 <=> strongly connected: {hc_zero_sc == 0 and hc_pos_notsc == 0}")
print()

# ================================================================
# 4. THE ALGEBRAIC FORMULA FOR HC
# ================================================================
print("  4. ALGEBRAIC FORMULA FOR HC")
print()

# H = C_n - S_2 at n<=4 (exact). What about HC?
# HC at n=3: H=1->HC=0, H=3->HC=1. So HC = (H-1)/3 = c3.
# HC at n=4: H=1->HC=0, H=3->HC=0, H=5->HC=1.
# HC(n=4) in {0, 0, 1}. Is there a score formula?

# At n=3: HC = c3 (number of 3-cycles). Since there's only 1 possible 3-cycle.
# At n=4: HC depends on... let's check.

# Build data for n=4
print("  n=4 HC analysis:")
for n_test in [4]:
    pairs_t = [(i,j) for i in range(n_test) for j in range(i+1, n_test)]
    m_t = len(pairs_t)
    score_HC_4 = defaultdict(set)
    for bits in range(2**m_t):
        A = np.zeros((n_test, n_test), dtype=np.int8)
        for k, (i,j) in enumerate(pairs_t):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        s = tuple(sorted(A.sum(axis=1).astype(int)))
        HC = count_hc(A, n_test)
        c3 = c3_from_scores(n_test, A.sum(axis=1).astype(int))
        score_HC_4[(s, c3)].add(HC)

    for key, hcs in sorted(score_HC_4.items()):
        s, c3 = key
        print(f"    scores={list(s)}, c3={c3}: HC in {sorted(hcs)}")

print()

# ================================================================
# 5. HC DISTRIBUTION AT n=5
# ================================================================
print("  5. FULL (H, HC) DISTRIBUTION AT n=5")
print()

pair_counts = defaultdict(int)
for d in all_data:
    pair_counts[(d['H'], d['HC'])] += 1

print(f"  {'H':>4s} {'HC':>4s} {'L':>4s} {'count':>6s} {'c3 range':>10s} {'SC?':>5s}")
for (H, HC) in sorted(pair_counts.keys()):
    count = pair_counts[(H, HC)]
    L = H - n * HC
    c3s = set(d['c3'] for d in all_data if d['H'] == H and d['HC'] == HC)
    scs = set(d['sc'] for d in all_data if d['H'] == H and d['HC'] == HC)
    print(f"  {H:>4d} {HC:>4d} {L:>4d} {count:>6d} {str(sorted(c3s)):>10s} {str(sorted(scs)):>5s}")

print()

# ================================================================
# 6. THE L INVARIANT
# ================================================================
print("  6. L = H - n*HC AS AN INVARIANT")
print()

L_dist = defaultdict(int)
for d in all_data:
    L_dist[d['L']] += 1

print("  L distribution at n=5:")
for L_val in sorted(L_dist.keys()):
    count = L_dist[L_val]
    pct = 100 * count / len(all_data)
    print(f"    L={L_val}: {count} tournaments ({pct:.1f}%)")

print()
# Is L always odd?
L_vals = set(d['L'] for d in all_data)
print(f"  L values: {sorted(L_vals)}")
print(f"  All L odd: {all(l % 2 == 1 for l in L_vals)}")
print(f"  L = 0 count: {sum(1 for d in all_data if d['L'] == 0)}")
print()

# L=0 means every HP is part of a cycle. These are the MOST CYCLIC tournaments.
print("  L=0 tournaments (every HP in a cycle):")
for d in all_data:
    if d['L'] == 0:
        if d == [x for x in all_data if x['L']==0][0]:  # first example
            print(f"    H={d['H']}, HC={d['HC']}, scores={d['score']}, c3={d['c3']}")
            break

# Count how many
l0_count = sum(1 for d in all_data if d['L'] == 0)
l0_H = sorted(set(d['H'] for d in all_data if d['L'] == 0))
print(f"    Total L=0: {l0_count} tournaments")
print(f"    H values with L=0: {l0_H}")
print()

# ================================================================
# 7. CORRELATION BETWEEN HC AND OTHER INVARIANTS
# ================================================================
print("  7. CORRELATIONS")
print()

c3_arr = np.array([d['c3'] for d in all_data], dtype=float)
L_arr = np.array([d['L'] for d in all_data], dtype=float)
S2_arr = np.array([d['S2'] for d in all_data], dtype=float)

corrs = {
    'H vs HC': np.corrcoef(H_arr, HC_arr)[0,1],
    'H vs c3': np.corrcoef(H_arr, c3_arr)[0,1],
    'HC vs c3': np.corrcoef(HC_arr, c3_arr)[0,1],
    'H vs L': np.corrcoef(H_arr, L_arr)[0,1],
    'HC vs L': np.corrcoef(HC_arr, L_arr)[0,1],
    'H vs S2': np.corrcoef(H_arr, S2_arr)[0,1],
    'HC vs S2': np.corrcoef(HC_arr, S2_arr)[0,1],
    'L vs S2': np.corrcoef(L_arr, S2_arr)[0,1],
}

for name, corr in sorted(corrs.items(), key=lambda x: -abs(x[1])):
    print(f"    {name:>12s}: {corr:>+.4f}")

print()
print("  KEY: HC correlates VERY strongly with H (as expected).")
print("  But L has a DIFFERENT correlation structure than H or HC.")
print("  L vs S2 tells us: how much of the NON-CYCLIC part is score-determined?")
