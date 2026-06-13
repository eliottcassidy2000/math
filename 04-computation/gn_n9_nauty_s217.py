#!/usr/bin/env python3
"""
gn_n9_nauty_s217.py — Compute |E(G_9)| using nauty
opus-2026-03-22-S217

191536 iso classes at n=9, 36 arc positions.
~6.9M lookups needed for edge computation.
"""

import sys
import subprocess
from math import comb
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)
n = 9; m = comb(n, 2)  # 36
PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]

print(f"G_{n} edge computation via nauty")
print(f"n={n}, m={m}")

# Read from nauty
t0 = time.time()
result = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
lines = [l.strip() for l in (result.stdout or '').split('\n')
         if len(l.strip()) == m and all(c in '01' for c in l.strip())]
print(f"Read {len(lines)} classes in {time.time()-t0:.2f}s")

def bits_to_adj(bits):
    adj = [[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(PAIRS):
        if bits & (1 << (m-1-k)): adj[i][j] = 1
        else: adj[j][i] = 1
    return adj

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

def H_dp(adj):
    dp = {}
    for v in range(n): dp[(1<<v,v)] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not(S&(1<<v)): continue
            val = dp.get((S,v),0)
            if val==0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]:
                    key=(S|(1<<u),u)
                    dp[key] = dp.get(key,0) + val
    return sum(dp.get(((1<<n)-1,v),0) for v in range(n))

def canonical_form(adj):
    """Degree-pruned canonical form."""
    scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    score_groups = defaultdict(list)
    for v in range(n): score_groups[scores[v]].append(v)
    groups = [score_groups[s] for s in sorted(set(scores))]
    if all(len(g)==1 for g in groups):
        perm = [g[0] for g in groups]
        return tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
    from itertools import permutations as perms
    best = None
    def gen_from_groups(gs):
        if not gs: yield []; return
        for p in perms(gs[0]):
            for rest in gen_from_groups(gs[1:]): yield list(p)+rest
    for perm in gen_from_groups(groups):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

# Phase 2: compute hash for each class
print(f"\nPhase 2: Computing (score, c3, H) for {len(lines)} classes...")
t1 = time.time()

hash_to_cids = defaultdict(list)
class_adj = []
class_bits = []

for cid, line in enumerate(lines):
    bits = int(line, 2)
    adj = bits_to_adj(bits)
    ss = score_seq(adj)
    c3 = c3_count(adj)
    h = H_dp(adj)
    hash_to_cids[(ss,c3,h)].append(cid)
    class_adj.append(adj)
    class_bits.append(bits)
    if (cid+1)%10000==0:
        print(f"  {cid+1}/{len(lines)} ({time.time()-t1:.1f}s)")

print(f"Phase 2 done: {len(hash_to_cids)} hash groups in {time.time()-t1:.1f}s")
ambiguous = sum(1 for v in hash_to_cids.values() if len(v)>1)
print(f"Ambiguous groups: {ambiguous}")

# For ambiguous groups, compute canonical forms
if ambiguous > 0:
    print(f"\nPhase 2b: Canonical forms for ambiguous groups...")
    t1b = time.time()
    canon_to_cid = {}
    ambig_cids = set()
    for cids in hash_to_cids.values():
        if len(cids) > 1:
            for c in cids: ambig_cids.add(c)

    done = 0
    for cid in ambig_cids:
        canon = canonical_form(class_adj[cid])
        canon_to_cid[canon] = cid
        done += 1
        if done % 1000 == 0:
            print(f"  {done}/{len(ambig_cids)} canons ({time.time()-t1b:.1f}s)")

    print(f"Phase 2b done: {len(ambig_cids)} canonical forms in {time.time()-t1b:.1f}s")
else:
    canon_to_cid = {}

# Phase 3: Compute edges
print(f"\nPhase 3: Computing edges...")
t2 = time.time()

def lookup_class(adj):
    ss = score_seq(adj)
    c3 = c3_count(adj)
    h = H_dp(adj)
    cids = hash_to_cids.get((ss,c3,h))
    if cids is None: return -1
    if len(cids)==1: return cids[0]
    canon = canonical_form(adj)
    return canon_to_cid.get(canon, -1)

edges = set()
self_loops = set()
deg_sum = 0

for cid in range(len(lines)):
    bits = class_bits[cid]
    nbrs = set()
    for arc in range(m):
        fb = bits ^ (1 << (m-1-arc))
        fa = bits_to_adj(fb)
        nb = lookup_class(fa)
        if nb == cid: self_loops.add(cid)
        elif nb >= 0: nbrs.add(nb)
    for nb in nbrs:
        edges.add((min(cid,nb), max(cid,nb)))
    deg_sum += len(nbrs)
    if (cid+1)%5000==0:
        elapsed = time.time()-t2
        rate = (cid+1)/elapsed
        eta = (len(lines)-cid-1)/rate
        print(f"  {cid+1}/{len(lines)} ({elapsed:.1f}s, ETA {eta:.0f}s, {len(edges)} edges)")

print(f"\nPhase 3 done in {time.time()-t2:.1f}s")

# Results
print(f"\n{'='*60}")
print(f"RESULTS: G_{n}")
print(f"{'='*60}")
print(f"  Vertices: {len(lines)}")
print(f"  EDGES: {len(edges)}")
print(f"  Self-loops: {len(self_loops)}")
print(f"  avg_deg: {deg_sum/len(lines):.2f}")
print(f"  avg_deg/m: {deg_sum/(len(lines)*m):.4f}")

T9 = 6847200
print(f"\n  T_9 = {T9}")
print(f"  T_9/E = {T9/len(edges):.4f}")
print(f"  T_9/(2E) = {T9/(2*len(edges)):.4f}")

print(f"\n*** |E(G_n)| SEQUENCE ***")
for nn,ee in [(3,1),(4,5),(5,30),(6,290),(7,4086),(8,91161),(9,len(edges))]:
    print(f"  n={nn}: {ee}")

print(f"\nTotal time: {time.time()-t0:.1f}s")
