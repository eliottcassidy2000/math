#!/usr/bin/env python3
"""
overhead_n7_v3_s20cn.py -- kind-pasteur-2026-03-22-S20cn

Compute SL_7 and MW_7 using the PROVEN approach from gn_edges_n7_fast_s212.py
(which successfully computed E(G_7) = 4086 in 32 minutes).

Key: use hash-table + canonical-form-fallback for class lookup.
Add SL/MW computation per class in Phase 3.

Author: kind-pasteur-2026-03-22-S20cn
"""
import sys, time
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

n = 7; m = comb(n, 2); total = 1 << m; nfact = factorial(n)
PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]
PAIR_INDEX = {}
for k,(i,j) in enumerate(PAIRS): PAIR_INDEX[(i,j)] = k

def bits_to_adj(bits):
    adj = [[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(PAIRS):
        if bits & (1<<k): adj[i][j] = 1
        else: adj[j][i] = 1
    return adj

def H_dp(adj, nn):
    dp = [0]*((1<<nn)*nn)
    for v in range(nn): dp[(1<<v)*nn+v] = 1
    for S in range(1,1<<nn):
        for v in range(nn):
            if not (S&(1<<v)): continue
            val = dp[S*nn+v]
            if val==0: continue
            for u in range(nn):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*nn+u] += val
    return sum(dp[((1<<nn)-1)*nn+v] for v in range(nn))

def score_seq(adj):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def c3_count(adj):
    c3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if adj[i][j] and adj[j][k] and adj[k][i]: c3+=1
                if adj[i][k] and adj[k][j] and adj[j][i]: c3+=1
    return c3

def deletion_fp(adj):
    h = []
    for v in range(n):
        vts = [i for i in range(n) if i!=v]
        nn = len(vts)
        sub = [[adj[vts[i]][vts[j]] for j in range(nn)] for i in range(nn)]
        h.append(H_dp(sub, nn))
    return tuple(sorted(h))

def canonical_form(adj):
    best = None
    for perm in permutations(range(n)):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

print("="*70)
print("  n=7 OVERHEAD DECOMPOSITION v3 (proven approach)")
print("="*70)

# ---- PHASE 1 ----
print(f"\n  Phase 1: fingerprinting {total} tournaments...")
t0 = time.time()
hash_groups = defaultdict(list)
for bits in range(total):
    adj = bits_to_adj(bits)
    ss = score_seq(adj)
    c3 = c3_count(adj)
    dfp = deletion_fp(adj)
    hash_groups[(ss,c3,dfp)].append(bits)
    if (bits+1) % 200000 == 0:
        r = (bits+1)/(time.time()-t0)
        print(f"    {bits+1}/{total} ({r:.0f}/s, ETA {(total-bits-1)/r:.0f}s)")
print(f"  Phase 1 done: {len(hash_groups)} groups in {time.time()-t0:.1f}s")

# ---- PHASE 2 ----
print(f"\n  Phase 2: identifying classes (canonical form for splits)...")
t1 = time.time()
valid_sizes = {nfact//a for a in range(1,nfact+1) if nfact%a==0}

hash_to_cid = {}  # hkey -> int (single class) or dict (split)
class_reps = {}; class_sizes = {}; cid = 0
splits_done = 0; split_members = 0

for hkey, members in hash_groups.items():
    sz = len(members)
    if sz in valid_sizes:
        hash_to_cid[hkey] = cid
        class_reps[cid] = members[0]; class_sizes[cid] = sz
        cid += 1
    else:
        splits_done += 1; split_members += sz
        sub = {}
        for b in members:
            adj = bits_to_adj(b)
            cf = canonical_form(adj)
            if cf not in sub:
                sub[cf] = cid
                class_reps[cid] = b; class_sizes[cid] = 0
                cid += 1
            class_sizes[sub[cf]] += 1
        hash_to_cid[hkey] = sub
        if splits_done % 10 == 0:
            print(f"    {splits_done} groups split ({split_members} members), {cid} classes so far, {time.time()-t1:.0f}s")

N = cid
print(f"  Phase 2 done: {N} classes ({splits_done} splits, {split_members} members) in {time.time()-t1:.1f}s")

# ---- PHASE 3: SL/MW DECOMPOSITION ----
print(f"\n  Phase 3: SL/MW decomposition for {N} classes...")
t2 = time.time()

def lookup_class(bits):
    adj = bits_to_adj(bits)
    ss = score_seq(adj); c3 = c3_count(adj); dfp = deletion_fp(adj)
    entry = hash_to_cid.get((ss,c3,dfp))
    if entry is None: return -1
    if isinstance(entry, int): return entry
    cf = canonical_form(adj)
    return entry.get(cf, -1)

total_SL = 0; total_MW = 0; total_deg = 0; total_orb = 0

for ci in range(N):
    bits = class_reps[ci]
    adj = bits_to_adj(bits)
    aut_size = nfact // class_sizes[ci]

    if aut_size == 1:
        # Each arc is its own orbit
        self_arcs = 0; neighbor_map = Counter()
        for k in range(m):
            nb = lookup_class(bits ^ (1<<k))
            if nb == ci: self_arcs += 1
            else: neighbor_map[nb] += 1

        sl = self_arcs
        deg = len(neighbor_map)
        mw = sum(v-1 for v in neighbor_map.values())  # excess arcs per neighbor
        orb = m  # each arc = one orbit
    else:
        # Compute Aut(T) and arc orbits
        aut_perms = []
        for perm in permutations(range(n)):
            ok = True
            for i in range(n):
                for j in range(n):
                    if adj[perm[i]][perm[j]] != adj[i][j]: ok=False; break
                if not ok: break
            if ok: aut_perms.append(perm)

        # Arc orbits
        orbit_id = [-1]*m; oc = 0
        for k in range(m):
            if orbit_id[k]>=0: continue
            orb_set = set()
            for perm in aut_perms:
                i,j = PAIRS[k]; pi,pj = perm[i],perm[j]
                nk = PAIR_INDEX[(min(pi,pj),max(pi,pj))]
                orb_set.add(nk)
            for a in orb_set: orbit_id[a] = oc
            oc += 1

        orbit_rep = {}
        for k in range(m):
            if orbit_id[k] not in orbit_rep: orbit_rep[orbit_id[k]] = k

        sl = 0; neighbors = set()
        for oid, rep in orbit_rep.items():
            nb = lookup_class(bits ^ (1<<rep))
            if nb == ci: sl += 1
            else: neighbors.add(nb)

        deg = len(neighbors)
        cross = oc - sl
        mw = cross - deg
        orb = oc

    total_SL += sl; total_MW += mw; total_deg += deg; total_orb += orb

    if (ci+1) % 50 == 0:
        el = time.time()-t2
        print(f"    {ci+1}/{N} ({(ci+1)/el:.1f}/s, ETA {(N-ci-1)*el/(ci+1):.0f}s) SL={total_SL} MW={total_MW}")

print(f"  Phase 3 done in {time.time()-t2:.1f}s")

E = total_deg // 2
OH = total_SL + total_MW
print(f"\n{'='*70}")
print(f"  RESULTS n=7")
print(f"{'='*70}")
print(f"  V = {N} (expected 456)")
print(f"  T = {total_orb} (expected 8912)")
print(f"  E = {E} (expected 4086)")
print(f"  OH = {OH} (expected 740)")
print(f"  SL = {total_SL} (predicted 326)")
print(f"  MW = {total_MW} (predicted 414)")
print(f"  Check: SL+2E+MW = {total_SL+total_deg+total_MW} = T? {total_SL+total_deg+total_MW == total_orb}")
print(f"\n  SL match: {total_SL == 326}")
print(f"  MW match: {total_MW == 414}")
print(f"  Total time: {time.time()-t0:.1f}s ({(time.time()-t0)/60:.1f} min)")
