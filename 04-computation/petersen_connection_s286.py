#!/usr/bin/env python3
"""
petersen_connection_s286.py — Is the Petersen graph hiding in the metagraph?
opus-2026-03-24-S286

The Petersen graph has:
  10 vertices, 15 edges, 3-regular, girth 5
  = Kneser graph K(5,2) = complement of Johnson graph J(5,2)

The n=5 merged metagraph has:
  10 vertices, 21 edges

21 - 15 = 6 = C(4,2) = number of tiles at n=5!

Questions:
1. Is the Petersen graph a subgraph of the n=5 merged metagraph?
2. Which 6 edges would need to be removed?
3. Do those 6 edges correspond to the 6 wiggly classes?
4. Does the Petersen structure relate to the 5 collapsed edges at n=6?
"""

import sys
import subprocess
from math import comb
from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

# Build the Petersen graph
def petersen_graph():
    """The Petersen graph on 10 vertices."""
    # Standard labeling: outer = 0-4, inner = 5-9
    # Outer pentagon: 0-1-2-3-4-0
    # Inner pentagram: 5-7-9-6-8-5
    # Spokes: 0-5, 1-6, 2-7, 3-8, 4-9
    edges = set()
    # Outer pentagon
    for i in range(5):
        edges.add((i, (i+1)%5))
    # Inner pentagram
    for i in range(5):
        edges.add((5+i, 5+(i+2)%5))
    # Spokes
    for i in range(5):
        edges.add((i, 5+i))
    # Normalize
    return {(min(a,b), max(a,b)) for (a,b) in edges}

def build_n5_metagraph():
    """Build the n=5 merged metagraph."""
    n = 5; m_total = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1, n)]
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    lines_raw = [l.strip() for l in (r.stdout or '').split('\n')
                 if len(l.strip()) == m_total and all(c in '01' for c in l.strip())]

    def b2a(bits):
        a = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if bits & (1 << (m_total-1-k)): a[i][j] = 1
            else: a[j][i] = 1
        return a
    def ss(a): return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))
    def c3(a):
        c = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if a[i][j] and a[j][k] and a[k][i]: c += 1
                    if a[i][k] and a[k][j] and a[j][i]: c += 1
        return c
    def H(a):
        dp = {}
        for v in range(n): dp[(1<<v,v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S&(1<<v)): continue
                val = dp.get((S,v), 0)
                if val == 0: continue
                for u in range(n):
                    if S&(1<<u): continue
                    if a[v][u]: dp[(S|(1<<u),u)] = dp.get((S|(1<<u),u),0) + val
        return sum(dp.get(((1<<n)-1,v),0) for v in range(n))
    def cp(a):
        sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
        sg = defaultdict(list)
        for v in range(n): sg[sc[v]].append(v)
        gs = [sg[s] for s in sorted(sg.keys())]
        best = None; cnt = 0
        def gp(gs2):
            if not gs2: yield []; return
            for p in permutations(gs2[0]):
                for rest in gp(gs2[1:]): yield list(p)+rest
        for p in gp(gs):
            f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or f < best: best = f
            cnt += 1
            if cnt > 500000: break
        return best

    cls = []; hm = defaultdict(list); cm = {}
    for i, line in enumerate(lines_raw):
        bits = int(line, 2); adj = b2a(bits)
        s = ss(adj); cc = c3(adj); h = H(adj); canon = cp(adj)
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid':i, 'score':s, 'c3':cc, 'H':h, 'sc':sc,
                    'comp_canon':cc2, 'canon':canon})
        hm[(s,cc,h)].append(i); cm[canon] = i
    for cl in cls: cl['comp_cid'] = cm.get(cl['comp_canon'], -1)

    def lk(adj):
        s = ss(adj); cc = c3(adj); h = H(adj)
        cids = hm.get((s,cc,h))
        if not cids: return -1
        if len(cids) == 1: return cids[0]
        return cm.get(cp(adj), -1)

    mn = {}; mid = 0
    for cl in cls:
        ci = cl['cid']
        if ci in mn: continue
        if cl['sc']: mn[ci] = mid; mid += 1
        else:
            mn[ci] = mid; comp = cl['comp_cid']
            if comp >= 0 and comp != ci: mn[comp] = mid
            mid += 1
    mi = {}
    for mm in range(mid):
        cids = [cl['cid'] for cl in cls if mn[cl['cid']] == mm]
        cl0 = cls[cids[0]]
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc']}

    # Build metagraph edges
    me = set()
    for cl in cls:
        b = int(''.join(str(cl['score']), ), 2) if False else 0  # placeholder
    # Actually rebuild edges from arc flips
    for cl in cls:
        bits = int(lines_raw[cl['cid']], 2)
        for arc in range(m_total):
            fb = bits ^ (1 << (m_total-1-arc))
            fa = b2a(fb)
            nb = lk(fa)
            if nb >= 0 and nb != cl['cid']:
                a2, b2 = mn[cl['cid']], mn[nb]
                if a2 != b2:
                    me.add((min(a2,b2), max(a2,b2)))

    return mi, me, cls, mn

print("=" * 72)
print("  PETERSEN GRAPH AND THE n=5 METAGRAPH")
print("  opus-2026-03-24-S286")
print("=" * 72)

# Build both graphs
pet = petersen_graph()
mi, me, cls, mn = build_n5_metagraph()

print(f"\n  Petersen: {len(pet)} edges")
print(f"  n=5 merged metagraph: {len(me)} edges")
print(f"  Difference: {len(me)} - {len(pet)} = {len(me) - len(pet)}")

# Node info
print(f"\n  MERGED METAGRAPH NODES:")
for mm in sorted(mi.keys()):
    deg = sum(1 for (a,b) in me if a == mm or b == mm)
    print(f"    Node {mm}: H={mi[mm]['H']}, SC={'SC' if mi[mm]['sc'] else 'NS'}, degree={deg}")

# Degree sequence
degs = [sum(1 for (a,b) in me if a == mm or b == mm) for mm in range(10)]
print(f"\n  Degree sequence: {sorted(degs)}")
print(f"  Petersen degree sequence: {[3]*10}")

# Check: is Petersen a subgraph?
print(f"\n  CHECKING IF PETERSEN IS A SUBGRAPH...")
# Try all 10! = 3628800 labelings... too many.
# Use degree-based pruning: Petersen is 3-regular.
# Nodes with degree < 3 in metagraph can't map to Petersen vertices.
# All nodes have degree >= 2, most have degree >= 3.

# More efficient: check all 10!/|Aut(Petersen)| = 3628800/120 = 30240 essentially distinct labelings
# Still a lot. Let me use a smarter approach.

# Build adjacency matrices
A_meta = np.zeros((10, 10), dtype=int)
for (a, b) in me:
    A_meta[a][b] = 1; A_meta[b][a] = 1

A_pet = np.zeros((10, 10), dtype=int)
for (a, b) in pet:
    A_pet[a][b] = 1; A_pet[b][a] = 1

# Check subgraph by trying permutations (brute force for 10 nodes)
# Petersen has |Aut| = 120, so we need to check 10!/120 = 30240 maps
# But let's use degree pruning
meta_degs = sorted(degs)
pet_degs = [3] * 10

# Nodes with degree 2 in metagraph can't be Petersen vertices (need deg >= 3)
low_deg = [mm for mm in range(10) if degs[mm] < 3]
high_deg = [mm for mm in range(10) if degs[mm] >= 3]
print(f"  Nodes with degree < 3: {low_deg} (can't be Petersen vertices)")
print(f"  Nodes with degree >= 3: {high_deg} ({len(high_deg)} nodes)")

if len(high_deg) >= 10:
    print(f"  All 10 nodes have degree >= 3, so Petersen subgraph is possible.")
else:
    print(f"  Only {len(high_deg)} nodes have degree >= 3, need 10.")
    print(f"  Petersen subgraph is IMPOSSIBLE (not enough high-degree nodes).")
    if len(high_deg) < 10:
        # Check with nodes of degree >= 3
        pass

# Try all permutations of the 10 nodes (3.6M — feasible in C but slow in Python)
# Let's use a smarter approach: check if the complement has the right structure
# Metagraph has 21 edges. Petersen has 15. Complement of Petersen in metagraph = 6 edges.
# These 6 edges form the "excess" — the edges in metagraph but not in Petersen.
# If Petersen is a subgraph, these 6 excess edges should form a specific pattern.

# Alternative: check if the metagraph contains a 3-regular subgraph
print(f"\n  CHECKING FOR 3-REGULAR SUBGRAPHS ON 10 VERTICES...")

# A 3-regular subgraph has 15 edges (since 10*3/2 = 15 = 21 - 6)
# We need to remove 6 edges to get a 3-regular graph
# The 6 removed edges must reduce each node's degree to exactly 3

# Current degrees: let's see what reductions are needed
target = 3
excess = [degs[i] - target for i in range(10)]
total_excess = sum(excess)
print(f"  Current degrees: {degs}")
print(f"  Excess over 3: {excess}")
print(f"  Total excess: {total_excess} (need {total_excess//2} = {total_excess/2} edges removed)")

if total_excess == 12:  # 6 edges × 2 endpoints = 12
    print(f"  ✓ Total excess = 12 = 2 × 6 → exactly 6 edges must be removed")

    # Find all sets of 6 edges whose removal makes the graph 3-regular
    found_count = 0
    all_edges_list = sorted(me)

    for removal in combinations(range(len(all_edges_list)), 6):
        removed = {all_edges_list[i] for i in removal}
        remaining = me - removed

        # Check if remaining is 3-regular
        new_degs = [0] * 10
        for (a, b) in remaining:
            new_degs[a] += 1; new_degs[b] += 1

        if all(d == 3 for d in new_degs):
            found_count += 1
            if found_count <= 5:
                removed_sorted = sorted(removed)
                # Check if this is the Petersen graph
                # Build adjacency and check girth
                A_rem = np.zeros((10, 10), dtype=int)
                for (a, b) in remaining:
                    A_rem[a][b] = 1; A_rem[b][a] = 1

                # Girth (length of shortest cycle)
                girth = 0
                for start in range(10):
                    # BFS
                    from collections import deque
                    q = deque([(start, -1, 0)])
                    visited = {start: 0}
                    while q:
                        node, parent, dist = q.popleft()
                        for nb in range(10):
                            if A_rem[node][nb] and nb != parent:
                                if nb in visited:
                                    cycle_len = dist + 1 + visited[nb]
                                    if girth == 0 or cycle_len < girth:
                                        girth = cycle_len
                                else:
                                    visited[nb] = dist + 1
                                    q.append((nb, node, dist + 1))

                # Check Petersen-ness: 3-regular, girth 5
                is_petersen_candidate = (girth == 5)

                removed_H = [(mi[a]['H'], mi[b]['H']) for (a,b) in removed_sorted]
                print(f"\n    3-regular subgraph #{found_count}:")
                print(f"    Removed edges: {removed_sorted}")
                print(f"    Removed H-pairs: {removed_H}")
                print(f"    Girth: {girth}")
                print(f"    Petersen candidate (girth=5): {is_petersen_candidate}")

                # Check spectrum
                evals = sorted(np.linalg.eigvalsh(A_rem.astype(float)))
                pet_evals = sorted(np.linalg.eigvalsh(A_pet.astype(float)))
                print(f"    Spectrum: {[f'{x:.2f}' for x in evals]}")
                print(f"    Petersen: {[f'{x:.2f}' for x in pet_evals]}")
                spec_match = np.allclose(sorted(evals), sorted(pet_evals), atol=0.01)
                print(f"    Spectrum match: {spec_match}")

    print(f"\n  Total 3-regular subgraphs found: {found_count}")
else:
    print(f"  Total excess ≠ 12, no exact 3-regular subgraph possible.")

print("\nDONE.")
