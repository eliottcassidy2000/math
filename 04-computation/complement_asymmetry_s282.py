#!/usr/bin/env python3
"""
complement_asymmetry_s282.py -- opus-2026-03-24-S282

MAXIMAL ASYMMETRY WITHIN PERFECT COMPLEMENT PAIRING

Every tournament T has a complement T^op (reverse all arcs).
H(T) = H(T^op), score(T) = reverse(score(T^op)), c3(T) = c3(T^op).

But the WIGGLY STRUCTURE can differ:
- Which arcs are neutral in T may differ from T^op
- The wiggly profile w(T) need not equal w(T^op)
- The wiggly neighborhood of T in Q_m may look completely different

QUESTION: How asymmetric can T and T^op be in their wiggly structure?

For SC tournaments: T ≅ T^op, so the wiggly profiles are EQUIVALENT
(related by the automorphism that maps T to T^op).

For NS tournaments: T ≇ T^op, so the wiggly profiles can differ.
The COMPLEMENT ASYMMETRY = difference between wiggly profiles of
T and T^op, projected to the meta-graph.

Author: opus-2026-03-24-S282
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def Hcount(A, n):
    dp = {}
    for v in range(n): dp[(1<<v,v)] = 1
    full = (1<<n)-1
    for S in range(1,1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            c = dp.get((S,v),0)
            if c==0: continue
            for w in range(n):
                if S&(1<<w): continue
                if A[v][w]: dp[(S|(1<<w),w)] = dp.get((S|(1<<w),w),0) + c
    return sum(dp.get((full,v),0) for v in range(n))

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

n = 5; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

banner("COMPLEMENT ASYMMETRY AT n=%d" % n)

cmap = {}; sizes = {}; tc = {}; reps = {}
for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for k,(ii,jj) in enumerate(pairs):
        if (bits>>k)&1: A[ii][jj]=1
        else: A[jj][ii]=1
    c = canon(A, n)
    if c not in cmap:
        cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
        reps[cid] = [row[:] for row in A]
    sizes[cmap[c]] += 1; tc[bits] = cmap[c]
nc = len(cmap)

info = {}
for cid in range(nc):
    A = reps[cid]
    comp_c = canon(complement(A, n), n)
    info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c],
                 'SC': cmap[comp_c] == cid, 'size': sizes[cid],
                 'aut': factorial(n)//sizes[cid]}

# For each LABELED tournament: compute wiggly profile and complement's profile
print("  For each labeled tiling T:")
print("  - w(T) = which arcs are neutral (wiggly self-loop)")
print("  - w(T^op) = which arcs are neutral in the complement")
print("  - asymmetry = Hamming distance between w(T) and w(complement(T))\n")

# Merged nodes
merged = {}; mid_map = {}; mid = 0
for cid in range(nc):
    if cid in mid_map: continue
    comp = info[cid]['comp']
    merged[mid] = (cid,) if comp == cid else (cid, comp)
    mid_map[cid] = mid
    if comp != cid: mid_map[comp] = mid
    mid += 1
nm = mid

# For each tiling: compute wiggly profile AND complement's wiggly profile
asymmetries = defaultdict(list)  # merged_id -> list of asymmetry values

comp_mask = (2**m - 1)  # all-ones = complement in bit representation

for bits in range(2**m):
    mi = mid_map[tc[bits]]
    comp_bits = bits ^ comp_mask  # complement tiling

    # Wiggly profile of T
    w_T = []
    for k in range(m):
        flipped = bits ^ (1 << k)
        w_T.append(1 if mid_map[tc[flipped]] == mi else 0)

    # Wiggly profile of T^op (complement)
    w_Top = []
    for k in range(m):
        flipped_comp = comp_bits ^ (1 << k)
        mi_comp = mid_map[tc[comp_bits]]
        w_Top.append(1 if mid_map[tc[flipped_comp]] == mi_comp else 0)

    # Asymmetry = Hamming distance between w_T and w_Top
    hamming = sum(a != b for a, b in zip(w_T, w_Top))

    # Also: do T and T^op land in the same neighbors?
    T_targets = set()
    Top_targets = set()
    for k in range(m):
        T_targets.add(mid_map[tc[bits ^ (1 << k)]])
        Top_targets.add(mid_map[tc[comp_bits ^ (1 << k)]])

    target_overlap = len(T_targets & Top_targets)
    target_union = len(T_targets | Top_targets)

    asymmetries[mi].append({
        'hamming': hamming,
        'neutral_T': sum(w_T),
        'neutral_Top': sum(w_Top),
        'target_overlap': target_overlap,
        'target_union': target_union,
        'jaccard': target_overlap / target_union if target_union > 0 else 1,
    })

print("  COMPLEMENT ASYMMETRY PER MERGED CLASS:\n")
print("  %-3s | %-3s | %-3s | %-8s | %-8s | %-8s | %-8s | %-8s" %
      ("mid", "H", "SC?", "avg_hamm", "avg_n_T", "avg_n_Top", "Jaccard", "same_n?"))
print("  " + "-"*70)

for mi in range(nm):
    data = asymmetries[mi]
    H = info[merged[mi][0]]['H']
    is_sc = info[merged[mi][0]]['SC']
    avg_hamm = np.mean([d['hamming'] for d in data])
    avg_nT = np.mean([d['neutral_T'] for d in data])
    avg_nTop = np.mean([d['neutral_Top'] for d in data])
    avg_jacc = np.mean([d['jaccard'] for d in data])
    same_neutral = avg_nT == avg_nTop

    print("  %-3d | %-3d | %-3s | %-8.2f | %-8.2f | %-8.2f | %-8.4f | %s" %
          (mi, H, "SC" if is_sc else "NS", avg_hamm, avg_nT, avg_nTop,
           avg_jacc, "YES" if same_neutral else "NO"))

banner("THE ASYMMETRY SPECTRUM")

# Distribution of Hamming distances
all_hammings = []
for mi in range(nm):
    for d in asymmetries[mi]:
        all_hammings.append(d['hamming'])

hamm_dist = Counter(all_hammings)
print("  Hamming distance distribution (T vs T^op wiggly profiles):\n")
for h in sorted(hamm_dist.keys()):
    pct = 100 * hamm_dist[h] / len(all_hammings)
    bar = "#" * int(pct / 2)
    print("    d=%d: %4d tilings (%5.1f%%) %s" % (h, hamm_dist[h], pct, bar))

avg_global = np.mean(all_hammings)
print("\n  Average Hamming distance: %.2f / %d" % (avg_global, m))
print("  = %.1f%% of arcs have different neutrality under complement" %
      (100 * avg_global / m))

banner("SC vs NS ASYMMETRY")

sc_hammings = []
ns_hammings = []
for mi in range(nm):
    is_sc = info[merged[mi][0]]['SC']
    for d in asymmetries[mi]:
        if is_sc:
            sc_hammings.append(d['hamming'])
        else:
            ns_hammings.append(d['hamming'])

print("  SC tilings: avg Hamming = %.2f" % (np.mean(sc_hammings) if sc_hammings else 0))
print("  NS tilings: avg Hamming = %.2f" % (np.mean(ns_hammings) if ns_hammings else 0))

if sc_hammings and ns_hammings:
    print("\n  SC LESS asymmetric? %s (%.2f < %.2f)" %
          (np.mean(sc_hammings) < np.mean(ns_hammings),
           np.mean(sc_hammings), np.mean(ns_hammings)))

banner("THE MEANING OF MAXIMAL ASYMMETRY")

print("""
  MAXIMAL ASYMMETRY means: T and T^op have completely DIFFERENT
  wiggly neighborhoods. Flipping the SAME tile in T vs T^op
  leads to DIFFERENT meta-graph destinations.

  PERFECT COMPLEMENT PAIRING means: T and T^op have the SAME H,
  the SAME score (reversed), the SAME c3 count. They are
  "structurally equivalent" from the complement perspective.

  THE TENSION: the complement pairing says T ≈ T^op globally,
  but the wiggly structure says T ≠ T^op locally.

  FOR SC TOURNAMENTS: T ≅ T^op, so there exists σ with σ(T) = T^op.
  This σ maps the wiggly profile of T to the wiggly profile of T^op.
  So: w(T^op) = σ(w(T)). The profiles are S_n-equivalent.
  Hamming distance can be nonzero (different LABELED positions neutral)
  but the MULTISET of neutrality patterns is preserved.

  FOR NS TOURNAMENTS: T ≇ T^op. No σ maps T to T^op.
  The wiggly profiles are GENUINELY DIFFERENT.
  The asymmetry is a REAL invariant of the complement pair.

  THE ASYMMETRY MEASURES "how much information is lost by merging":
  High asymmetry → merging loses a lot of wiggly structure.
  Low asymmetry → merging preserves most wiggly structure.
  Zero asymmetry → T and T^op are wiggly-identical (SC case).
""")

# Which class has the MOST asymmetric complement pair?
max_asym_class = max(range(nm), key=lambda mi: np.mean([d['hamming'] for d in asymmetries[mi]]))
print("  MOST ASYMMETRIC CLASS: mid=%d (H=%d, SC=%s, avg hamming=%.2f)" %
      (max_asym_class, info[merged[max_asym_class][0]]['H'],
       info[merged[max_asym_class][0]]['SC'],
       np.mean([d['hamming'] for d in asymmetries[max_asym_class]])))

min_asym_class = min(range(nm), key=lambda mi: np.mean([d['hamming'] for d in asymmetries[mi]]))
print("  LEAST ASYMMETRIC CLASS: mid=%d (H=%d, SC=%s, avg hamming=%.2f)" %
      (min_asym_class, info[merged[min_asym_class][0]]['H'],
       info[merged[min_asym_class][0]]['SC'],
       np.mean([d['hamming'] for d in asymmetries[min_asym_class]])))

print("\nDone. opus-2026-03-24-S282")
