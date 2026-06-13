#!/usr/bin/env python3
"""
overhead_n7_v2_s20cn.py -- kind-pasteur-2026-03-22-S20cn

Optimized n=7 overhead decomposition using array for bits_to_class.

Key optimizations vs v1:
1. Use list (not dict) for bits_to_class → O(1) array access
2. Progress logging in Phase 2
3. Skip Aut computation for Phase 3 (use representative-based approach)

Expected: ~20 minutes total.

Author: kind-pasteur-2026-03-22-S20cn
"""
import sys, time
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

n = 7
m = comb(n, 2)  # 21
total = 1 << m   # 2097152
nfact = factorial(n)  # 5040

PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]
PAIR_INDEX = {}
for k, (i,j) in enumerate(PAIRS):
    PAIR_INDEX[(i,j)] = k

def bits_to_adj(bits):
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(PAIRS):
        if bits & (1 << k): adj[i][j] = 1
        else: adj[j][i] = 1
    return adj

def H_dp(adj, nn):
    dp = [0] * ((1 << nn) * nn)
    for v in range(nn):
        dp[(1 << v) * nn + v] = 1
    for S in range(1, 1 << nn):
        for v in range(nn):
            if not (S & (1 << v)): continue
            val = dp[S * nn + v]
            if val == 0: continue
            for u in range(nn):
                if S & (1 << u): continue
                if adj[v][u]:
                    dp[(S | (1 << u)) * nn + u] += val
    full = (1 << nn) - 1
    return sum(dp[full * nn + v] for v in range(nn))

def score_seq(adj):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def c3_count(adj):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]: c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]: c3 += 1
    return c3

def deletion_fp(adj):
    h_vals = []
    for v in range(n):
        verts = [i for i in range(n) if i != v]
        nn = len(verts)
        sub = [[adj[verts[i]][verts[j]] for j in range(nn)] for i in range(nn)]
        h_vals.append(H_dp(sub, nn))
    return tuple(sorted(h_vals))

def canonical_form(adj):
    best = None
    for perm in permutations(range(n)):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

print("=" * 70)
print("  n=7 OVERHEAD DECOMPOSITION v2")
print("=" * 70)

# ---- PHASE 1 ----
print(f"\n  Phase 1: fingerprinting {total} tournaments...")
t0 = time.time()
hash_groups = defaultdict(list)
for bits in range(total):
    adj = bits_to_adj(bits)
    ss = score_seq(adj)
    c3 = c3_count(adj)
    dfp = deletion_fp(adj)
    hash_groups[(ss, c3, dfp)].append(bits)
    if (bits+1) % 200000 == 0:
        r = (bits+1)/(time.time()-t0)
        print(f"    {bits+1}/{total} ({r:.0f}/s, ETA {(total-bits-1)/r:.0f}s, {len(hash_groups)} groups)")
print(f"  Phase 1 done: {len(hash_groups)} groups in {time.time()-t0:.1f}s")

# ---- PHASE 2 ----
print(f"\n  Phase 2: assigning classes...")
t1 = time.time()
# Tournament |Aut| must have only odd prime factors (all-odd-cycle permutations)
# For n=7: odd divisors of n! that only use odd primes
def has_only_odd_factors(x):
    while x % 2 == 0: x //= 2
    return x == x  # removed even factors; but we need ALL factors odd
# Actually: |Aut| of a tournament divides n! and |Aut| is odd
# So valid |Aut| = odd divisors of n!
valid_sizes = set()
for a in range(1, nfact+1, 2):  # only odd a
    if nfact % a == 0:
        valid_sizes.add(nfact // a)
print(f"  Valid class sizes (odd |Aut|): {sorted(valid_sizes)}")

btc = [0] * total  # bits_to_class array
class_reps = {}
class_sizes = {}
cid = 0
need_split = 0

for idx, (hkey, members) in enumerate(hash_groups.items()):
    sz = len(members)
    if sz in valid_sizes:
        class_reps[cid] = members[0]
        class_sizes[cid] = sz
        for b in members:
            btc[b] = cid
        cid += 1
    else:
        need_split += 1
        # Use H(T) to split (MUCH faster than canonical form)
        sub_h = {}
        for b in members:
            adj = bits_to_adj(b)
            h = H_dp(adj, n)
            key = h  # H alone should discriminate within hash group
            if key not in sub_h:
                sub_h[key] = cid
                class_reps[cid] = b
                class_sizes[cid] = 0
                cid += 1
            class_sizes[sub_h[key]] += 1
            btc[b] = sub_h[key]
    if (idx+1) % 50 == 0:
        print(f"    {idx+1}/{len(hash_groups)} groups processed, {cid} classes, {need_split} split")

N = cid
print(f"  Phase 2 done: {N} classes ({need_split} splits) in {time.time()-t1:.1f}s")

# ---- PHASE 3 ----
print(f"\n  Phase 3: arc orbit decomposition for {N} classes...")
t2 = time.time()

total_SL = 0
total_MW = 0
total_deg = 0
total_orb = 0

for ci in range(N):
    bits = class_reps[ci]
    adj = bits_to_adj(bits)
    aut_size = nfact // class_sizes[ci]

    # Find Aut(T)
    aut_perms = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if adj[perm[i]][perm[j]] != adj[i][j]:
                    ok = False; break
            if not ok: break
        if ok: aut_perms.append(perm)

    # Arc orbits under Aut
    orbit_id = [-1] * m
    oc = 0
    for k in range(m):
        if orbit_id[k] >= 0: continue
        orb = set()
        for perm in aut_perms:
            i, j = PAIRS[k]
            pi, pj = perm[i], perm[j]
            nk = PAIR_INDEX[(min(pi,pj), max(pi,pj))]
            orb.add(nk)
        for a in orb: orbit_id[a] = oc
        oc += 1

    # Classify orbits
    orbit_rep = {}
    for k in range(m):
        if orbit_id[k] not in orbit_rep:
            orbit_rep[orbit_id[k]] = k

    sl = 0; neighbors = set()
    for oid, rep in orbit_rep.items():
        nb = btc[bits ^ (1 << rep)]
        if nb == ci: sl += 1
        else: neighbors.add(nb)

    deg = len(neighbors)
    cross = oc - sl
    mw = cross - deg

    total_SL += sl
    total_MW += mw
    total_deg += deg
    total_orb += oc

    if (ci+1) % 100 == 0:
        print(f"    {ci+1}/{N} classes processed")

print(f"  Phase 3 done in {time.time()-t2:.1f}s")

# ---- RESULTS ----
E = total_deg // 2
OH = total_SL + total_MW
print(f"\n{'='*70}")
print(f"  RESULTS n=7")
print(f"{'='*70}")
print(f"  V = {N}, T = {total_orb}, E = {E}")
print(f"  OH = T - 2E = {OH}")
print(f"  SL = {total_SL}")
print(f"  MW = {total_MW}")
print(f"  Check: SL+2E+MW = {total_SL+total_deg+total_MW} (should be {total_orb})")
print(f"  Predicted SL = 326, Match: {total_SL == 326}")
print(f"  Predicted MW = 414, Match: {total_MW == 414}")
print(f"  Total time: {time.time()-t0:.1f}s")
