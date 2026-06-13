#!/usr/bin/env python3
"""
unmerged_asymmetry_n7_s281.py — Push the flat fiber bundle to n=7
opus-2026-03-24-S281

At n=7: 456 unmerged classes, 272 merged, 2^15 = 32768 tilings, m=15 tiles.
Verify: complement is exact symmetry of sub-flows at n=7.
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_data(n):
    m_total = comb(n, 2)
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

    cls = []
    hm = defaultdict(list); cm = {}
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

    mn = {}; mid = 0; merge_contents = {}
    for cl in cls:
        ci = cl['cid']
        if ci in mn: continue
        if cl['sc']:
            mn[ci] = mid; merge_contents[mid] = [ci]; mid += 1
        else:
            mn[ci] = mid; merge_contents[mid] = [ci]
            comp = cl['comp_cid']
            if comp >= 0 and comp != ci:
                mn[comp] = mid; merge_contents[mid].append(comp)
            mid += 1
    mi = {}
    for mm in range(mid):
        cids = merge_contents[mm]
        cl0 = cls[cids[0]]
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'], 'n_unmerged': len(cids)}

    return P, cls, lk, mn, mid, mi, b2a, m_total, merge_contents

print("=" * 72)
print("  n=7: UNMERGED ASYMMETRY — FLAT FIBER BUNDLE TEST")
print("  opus-2026-03-24-S281")
print("=" * 72)

n = 7
t0 = time.time()
print(f"\n  Building class data...", end=' ', flush=True)
P, cls, lk, mn, V, mi, b2a, m_total, merge_contents = build_data(n)
print(f"done ({time.time()-t0:.1f}s). V_merged={V}, V_unmerged={len(cls)}")

m_tiles = comb(n-1, 2)
tiles = []
for y in range(1, n-1):
    for x in range(n, y+1, -1):
        if x - y >= 2:
            tiles.append((x, y))

tile_to_P = {}; base_P = set()
for k, (i, j) in enumerate(P):
    if j == i + 1:
        base_P.add(k)
    else:
        for t_idx, (tx, ty) in enumerate(tiles):
            if (ty-1, tx-1) == (i, j):
                tile_to_P[t_idx] = k; break

def tiling_to_bits(tb):
    full_bits = 0
    for k, (i, j) in enumerate(P):
        if k not in base_P:
            for ti in range(m_tiles):
                if tile_to_P.get(ti) == k:
                    bit = (tb >> (m_tiles - 1 - ti)) & 1
                    if bit: full_bits |= (1 << (m_total - 1 - k))
                    break
    return full_bits

total = 2 ** m_tiles
print(f"  Tilings: {total}, Tiles: {m_tiles}")

# Precompute tiling → unmerged class
print(f"  Precomputing tiling classes...", end=' ', flush=True)
t1 = time.time()
tiling_unmerged = [0] * total
for tb in range(total):
    fb = tiling_to_bits(tb)
    adj = b2a(fb)
    cid = lk(adj)
    tiling_unmerged[tb] = cid
    if tb % 10000 == 0 and tb > 0:
        print(f"{tb}...", end='', flush=True)
print(f" done ({time.time()-t1:.1f}s)")

# Sample NS-NS sub-flows (checking ALL would be too slow with full W_u)
# Instead: for each tiling, flip each tile, record (unmerged_src, unmerged_dst)
# Then check complement pairing

print(f"\n  Computing sub-flows and checking complement symmetry...")
t2 = time.time()

# Count sub-flows for NS-NS edges
# For each merged pair (i,j), track flows between unmerged components
subflow = defaultdict(lambda: defaultdict(int))  # (merged_i, merged_j) → {(unmerged_a, unmerged_b): count}

processed = 0
for tb in range(total):
    cid_a = tiling_unmerged[tb]
    if cid_a < 0: continue
    ma = mn[cid_a]

    for ti in range(m_tiles):
        flipped = tb ^ (1 << (m_tiles - 1 - ti))
        cid_b = tiling_unmerged[flipped]
        if cid_b < 0: continue
        mb = mn[cid_b]

        if ma != mb:  # cross-class only
            key = (min(ma, mb), max(ma, mb))
            subflow[key][(cid_a, cid_b)] += 1

    processed += 1
    if processed % 10000 == 0:
        print(f"    {processed}/{total}...", flush=True)

dt_sub = time.time() - t2
print(f"  Done ({dt_sub:.1f}s). {len(subflow)} merged edges with sub-flows.")

# Now check complement symmetry for NS-NS edges
ns_ns_edges = [(i,j) for (i,j) in subflow.keys()
               if not mi[i]['sc'] and not mi[j]['sc']
               and mi[i]['n_unmerged'] == 2 and mi[j]['n_unmerged'] == 2]

print(f"\n  NS-NS edges to check: {len(ns_ns_edges)}")

complement_preserving = 0
not_preserving = 0
total_checked = 0

for (i, j) in ns_ns_edges:
    A, Ac = merge_contents[i]
    B, Bc = merge_contents[j]

    ab = subflow[(i,j)].get((A, B), 0) + subflow[(i,j)].get((B, A), 0)
    abc = subflow[(i,j)].get((A, Bc), 0) + subflow[(i,j)].get((Bc, A), 0)
    acb = subflow[(i,j)].get((Ac, B), 0) + subflow[(i,j)].get((B, Ac), 0)
    acbc = subflow[(i,j)].get((Ac, Bc), 0) + subflow[(i,j)].get((Bc, Ac), 0)

    # Check: w(A→B) = w(A^c→B^c) AND w(A→B^c) = w(A^c→B)
    # Need to be careful: A→B includes both directions in unmerged symmetric W
    # Actually subflow records directed (src, dst), and W_u is symmetric
    # So subflow[(i,j)][(A,B)] counts A→B transitions where i<j in merged

    # Simpler: just get the 4 sub-flow totals directly
    w_ab = subflow[(i,j)].get((A, B), 0)
    w_abc = subflow[(i,j)].get((A, Bc), 0)
    w_acb = subflow[(i,j)].get((Ac, B), 0)
    w_acbc = subflow[(i,j)].get((Ac, Bc), 0)

    # Also count reverse direction
    w_ba = subflow[(i,j)].get((B, A), 0)
    w_bca = subflow[(i,j)].get((Bc, A), 0)
    w_bac = subflow[(i,j)].get((B, Ac), 0)
    w_bcac = subflow[(i,j)].get((Bc, Ac), 0)

    # Total per sub-edge (undirected)
    s_AB = w_ab + w_ba
    s_ABc = w_abc + w_bca
    s_AcB = w_acb + w_bac
    s_AcBc = w_acbc + w_bcac

    if s_AB == s_AcBc and s_ABc == s_AcB:
        complement_preserving += 1
    else:
        not_preserving += 1

    total_checked += 1

print(f"\n  COMPLEMENT SYMMETRY CHECK (NS-NS edges):")
print(f"  Complement-preserving: {complement_preserving}/{total_checked} "
      f"({100*complement_preserving/total_checked:.1f}%)")
print(f"  Not preserving: {not_preserving}/{total_checked}")

if complement_preserving == total_checked:
    print(f"\n  *** FLAT FIBER BUNDLE VERIFIED AT n=7! ***")
    print(f"  All {total_checked} NS-NS edges have perfect complement pairing.")

# Also check SC-NS edges
sc_ns_edges = [(i,j) for (i,j) in subflow.keys()
               if mi[i]['sc'] != mi[j]['sc']]
print(f"\n  SC-NS edges: {len(sc_ns_edges)}")

balanced = 0
imbalanced = 0
for (i, j) in sc_ns_edges:
    if mi[i]['sc']:
        sc_node, ns_node = i, j
    else:
        sc_node, ns_node = j, i

    if mi[ns_node]['n_unmerged'] < 2: continue
    SC = merge_contents[sc_node][0]
    C, Cc = merge_contents[ns_node][0], merge_contents[ns_node][1]

    w1 = subflow[(i,j)].get((SC, C), 0) + subflow[(i,j)].get((C, SC), 0)
    w2 = subflow[(i,j)].get((SC, Cc), 0) + subflow[(i,j)].get((Cc, SC), 0)

    if w1 == w2:
        balanced += 1
    else:
        imbalanced += 1

print(f"  SC-NS balanced (SC→C = SC→C^c): {balanced}/{balanced+imbalanced}")
if balanced == balanced + imbalanced:
    print(f"  *** SC-NS PERFECTLY BALANCED AT n=7! ***")

# Summary stats
sc_count = sum(1 for mm in mi if mi[mm]['sc'])
ns_count = V - sc_count
print(f"\n  SUMMARY n=7:")
print(f"  SC classes: {sc_count}, NS classes: {ns_count}")
print(f"  NS-NS edges: {len(ns_ns_edges)}, all complement-preserving: "
      f"{'YES' if not_preserving == 0 else 'NO'}")
print(f"  SC-NS edges: {len(sc_ns_edges)}, all balanced: "
      f"{'YES' if imbalanced == 0 else 'NO'}")

total_time = time.time() - t0
print(f"\n  Total time: {total_time:.1f}s")
print("\nDONE.")
