#!/usr/bin/env python3
"""
Algebraic mechanism behind odd-n hereditary maximality.

Key question: WHY does deleting ANY vertex from the n-maximizer (odd n)
give the (n-1)-maximizer?

Three hypotheses to test:
1. VERTEX-TRANSITIVITY: Paley T_p is vertex-transitive, so all deletions
   are isomorphic. If ONE deletion is optimal, ALL are.
   But: is it true that T_p - v is the (p-1)-maximizer?

2. H-SUM IDENTITY: H(T) = (1/n) * sum_v [H(T-v) + correction].
   If the correction is minimized by the maximizer, this might force
   each H(T-v) to be maximal.

3. DESIGN STRUCTURE: The 2-(p,3,(p+1)/4) BIBD from cyclic triples
   of T_p constrains the cycle structure of deletions.

Also: at even n, the maximizer is NOT vertex-transitive (not Paley).
The deletion H-values are NOT constant. Can we characterize WHICH
deletions give the highest H?

kind-pasteur-2026-03-06-S18g
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles, conflict_graph
from itertools import permutations

MAX_H = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661}

def score_seq(T):
    return tuple(sorted(sum(T[i]) for i in range(len(T))))

def build_paley(p):
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i][j] = 1
    return T

def delete_vertex(T, v):
    n = len(T)
    verts = [i for i in range(n) if i != v]
    return [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]

def find_all_ham_paths(T):
    n = len(T)
    paths = []
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[i+1]] for i in range(n-1)):
            paths.append(perm)
    return paths

# ============================================================
# HYPOTHESIS 1: Vertex-transitivity forces constant deletion H
# ============================================================
print("=" * 70)
print("HYPOTHESIS 1: VERTEX-TRANSITIVITY AND DELETION CONSTANCY")
print("=" * 70)

for p in [3, 5, 7]:
    T = build_paley(p)
    h = hamiltonian_path_count(T)
    del_hs = [hamiltonian_path_count(delete_vertex(T, v)) for v in range(p)]
    print(f"\nT_{p}: H={h}, max H({p-1})={MAX_H[p-1]}")
    print(f"  Deletion H-values: {del_hs}")
    print(f"  All equal? {len(set(del_hs)) == 1}")
    print(f"  All = max H({p-1})? {all(h == MAX_H[p-1] for h in del_hs)}")

# At n=4,6,8 (even), check if maximizer is vertex-transitive
print(f"\n{'='*70}")
print("EVEN n: MAXIMIZER DELETION SPECTRUM")
print("=" * 70)

for n in [4, 6]:
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    # Find a maximizer
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h == max_h:
            del_hs = [hamiltonian_path_count(delete_vertex(T, v)) for v in range(n)]
            s = score_seq(T)
            print(f"\nn={n}: H={h}, score={s}")
            print(f"  Deletion H-values: {sorted(del_hs, reverse=True)}")
            print(f"  Constant? {len(set(del_hs)) == 1}")
            print(f"  Any = max H({n-1})? {MAX_H[n-1] in del_hs}")

            # Score of each deletion
            del_scores = [score_seq(delete_vertex(T, v)) for v in range(n)]
            print(f"  Deletion scores: {[str(s) for s in set(del_scores)]}")
            break

# ============================================================
# HYPOTHESIS 2: H-sum identity and path counting
# ============================================================
print(f"\n{'='*70}")
print("HYPOTHESIS 2: H-SUM IDENTITY")
print("=" * 70)

# Key identity from path theory:
# Each Hamiltonian path P of T uses n-1 arcs. When we delete vertex v,
# P either:
#   (a) v is an endpoint -> removing v gives a Ham path of T-v
#   (b) v is internal -> removing v may or may not give a Ham path of T-v
#
# For case (a): every Ham path has 2 endpoints, so each path contributes
# to exactly 2 of the n deletion counts.
#
# For case (b): removing internal v from P = ...->u->v->w->... gives
# a Ham path of T-v iff arc u->w exists in T.
#
# So: sum_v H(T-v) = 2*H(T) + (number of "internal deletions that work")

for n in [5, 7]:
    if n == 5:
        T = build_paley(5)
    else:
        T = build_paley(7)

    h = hamiltonian_path_count(T)
    del_hs = [hamiltonian_path_count(delete_vertex(T, v)) for v in range(n)]
    h_sum = sum(del_hs)

    # Count internal deletions that work
    if n <= 7:
        paths = find_all_ham_paths(T)
        internal_works = 0
        for p in paths:
            for pos in range(1, n-1):  # internal positions
                v = p[pos]
                u, w = p[pos-1], p[pos+1]
                if T[u][w]:
                    internal_works += 1

        print(f"\nn={n}: H={h}")
        print(f"  sum H(T-v) = {h_sum}")
        print(f"  2*H = {2*h}")
        print(f"  Internal deletions that work: {internal_works}")
        print(f"  Check: sum = 2*H + internal = {2*h + internal_works} {'OK' if 2*h + internal_works == h_sum else 'FAIL'}")
        print(f"  Ratio sum/H = {h_sum/h:.4f}")
        print(f"  Internal per path = {internal_works/len(paths):.4f}")

# ============================================================
# HYPOTHESIS 3: The ratio sum_v H(T-v) / H(T) across all tournaments
# Is the Paley maximizer special in having a BALANCED deletion sum?
# ============================================================
print(f"\n{'='*70}")
print("HYPOTHESIS 3: DELETION-SUM RATIO ACROSS ALL TOURNAMENTS")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    # For each tournament, compute H and sum H(T-v)
    ratio_by_h = {}  # h -> list of ratios

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h == 0:
            continue
        del_hs = [hamiltonian_path_count(delete_vertex(T, v)) for v in range(n)]
        ratio = sum(del_hs) / h

        if h not in ratio_by_h:
            ratio_by_h[h] = []
        ratio_by_h[h].append(ratio)

    print(f"\nn={n}: Deletion-sum ratio = sum H(T-v) / H(T)")
    for h in sorted(ratio_by_h.keys(), reverse=True)[:8]:
        ratios = ratio_by_h[h]
        avg_r = sum(ratios) / len(ratios)
        min_r = min(ratios)
        max_r = max(ratios)
        marker = " *** MAX" if h == max_h else ""
        print(f"  H={h}: avg ratio={avg_r:.3f}, range=[{min_r:.3f},{max_r:.3f}], count={len(ratios)}{marker}")

# ============================================================
# KEY: The "internal deletion rate" -- what fraction of
# (path, internal vertex) pairs give valid shorter paths?
# ============================================================
print(f"\n{'='*70}")
print("INTERNAL DELETION RATE vs H")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
rates = []

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    if h == 0:
        continue

    paths = find_all_ham_paths(T)
    internal_works = 0
    total_internal = 0
    for p in paths:
        for pos in range(1, n-1):
            total_internal += 1
            v = p[pos]
            u, w = p[pos-1], p[pos+1]
            if T[u][w]:
                internal_works += 1

    rate = internal_works / total_internal if total_internal > 0 else 0
    rates.append((h, rate, internal_works, total_internal))

# Sort by H descending
rates.sort(key=lambda x: -x[0])
print(f"n={n}: H vs internal deletion rate")
seen_h = set()
for h, rate, iw, ti in rates:
    if h not in seen_h:
        seen_h.add(h)
        marker = " *** MAX" if h == MAX_H[n] else ""
        print(f"  H={h}: rate={rate:.4f} ({iw}/{ti}){marker}")

# ============================================================
# DEEPER: For the Paley T_5 maximizer, which internal deletions work?
# ============================================================
print(f"\n{'='*70}")
print("T_5 INTERNAL DELETION STRUCTURE")
print("=" * 70)

T5 = build_paley(5)
paths5 = find_all_ham_paths(T5)
print(f"T_5: H={len(paths5)}, score={score_seq(T5)}")

# For each vertex, count how many paths have it as internal AND deletion works
vertex_internal = {v: 0 for v in range(5)}
vertex_internal_works = {v: 0 for v in range(5)}

for p in paths5:
    for pos in range(1, 4):  # internal positions in 5-path
        v = p[pos]
        vertex_internal[v] += 1
        u, w = p[pos-1], p[pos+1]
        if T5[u][w]:
            vertex_internal_works[v] += 1

print("Per-vertex internal deletion stats:")
for v in range(5):
    rate = vertex_internal_works[v] / vertex_internal[v] if vertex_internal[v] > 0 else 0
    print(f"  v={v}: {vertex_internal_works[v]}/{vertex_internal[v]} = {rate:.4f}")

# By vertex-transitivity these should all be equal
print(f"All rates equal? {len(set(vertex_internal_works[v] for v in range(5))) == 1}")

# ============================================================
# ALGEBRAIC: Can we express internal deletion rate in terms of
# score sequence or cycle structure?
# ============================================================
print(f"\n{'='*70}")
print("INTERNAL DELETION RATE vs SCORE/CYCLES (n=5)")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
data = []

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    if h == 0:
        continue

    s = score_seq(T)
    cycles = find_odd_cycles(T)
    c3 = len([c for c in cycles if len(c) == 3])

    paths = find_all_ham_paths(T)
    internal_works = 0
    for p in paths:
        for pos in range(1, n-1):
            u, w = p[pos-1], p[pos+1]
            if T[u][w]:
                internal_works += 1

    del_hs = [hamiltonian_path_count(delete_vertex(T, v)) for v in range(n)]
    h_sum = sum(del_hs)

    data.append((h, s, c3, internal_works, h_sum))

# Group by score
from collections import defaultdict
by_score = defaultdict(list)
for h, s, c3, iw, hs in data:
    by_score[s].append((h, c3, iw, hs))

print("Score class analysis:")
for s in sorted(by_score.keys()):
    entries = by_score[s]
    max_h = max(e[0] for e in entries)
    for h, c3, iw, hs in sorted(entries, key=lambda x: -x[0])[:2]:
        marker = " *** MAX" if h == MAX_H[n] else ""
        print(f"  score={s}: H={h}, c3={c3}, internal_works={iw}, sum={hs}{marker}")

# ============================================================
# THE CRITICAL FORMULA:
# sum_v H(T-v) = 2*H(T) + sum_{paths} sum_{internal v} [arc u->w exists]
#
# For regular tournaments (Paley), each "gap arc" u->w has the same
# probability of existing. The regularity MAXIMIZES the internal rate.
# ============================================================
print(f"\n{'='*70}")
print("GAP ARC ANALYSIS")
print("=" * 70)

# For each path P = v0->v1->...->v4, the "gap arcs" are
# (v0,v2), (v1,v3), (v2,v4) -- the arcs that would be needed
# for internal deletions to work.
# For the maximizer, how many of these gap arcs exist?

T5 = build_paley(5)
paths5 = find_all_ham_paths(T5)

gap_counts = {0: 0, 1: 0, 2: 0, 3: 0}
for p in paths5:
    gaps = sum(1 for pos in range(1, 4) if T5[p[pos-1]][p[pos+1]])
    gap_counts[gaps] += 1

print(f"T_5: Gap arc distribution across {len(paths5)} paths:")
for g in sorted(gap_counts.keys()):
    print(f"  {g} gap arcs: {gap_counts[g]} paths ({100*gap_counts[g]/len(paths5):.1f}%)")

# Compare with all n=5 tournaments
print(f"\nComparison with all n=5 tournaments (same score):")
target_score = (2, 2, 2, 2, 2)
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    if score_seq(T) != target_score:
        continue
    h = hamiltonian_path_count(T)
    paths = find_all_ham_paths(T)
    if not paths:
        continue
    gap_counts_t = {0: 0, 1: 0, 2: 0, 3: 0}
    for p in paths:
        gaps = sum(1 for pos in range(1, 4) if T[p[pos-1]][p[pos+1]])
        gap_counts_t[gaps] += 1
    avg_gaps = sum(g * c for g, c in gap_counts_t.items()) / len(paths)
    marker = " *** MAX" if h == MAX_H[n] else ""
    print(f"  bits={bits}: H={h}, avg gaps={avg_gaps:.3f}, dist={dict(gap_counts_t)}{marker}")

print("\nDone.")
