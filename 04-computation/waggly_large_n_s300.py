#!/usr/bin/env python3
"""
waggly_large_n_s300.py — Deep waggly investigation at n=7 (and n=8 if feasible)
opus-2026-03-24-S300

At n=7: m=15 tiles, 2^15=32768 tilings, V=272 merged classes.
Compute the waggly filling sequence F_1 through F_k.

Also compute CREATIVE WAGGLY SUBSETS:
- Range flips: flip all tiles at same skip value
- Vertex-star flips: flip all tiles incident to a vertex
- Triangle flips: flip 3 tiles forming a 3-cycle
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations, combinations
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
    cls = []; hm = defaultdict(list); cm = {}
    for i, line in enumerate(lines_raw):
        bits = int(line, 2); adj = b2a(bits)
        s = ss(adj); cc = c3(adj); h = H(adj); canon = cp(adj)
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid':i, 'score':s, 'c3':cc, 'H':h, 'sc':sc, 'comp_canon':cc2, 'canon':canon})
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
        cl0 = cls[cids[0]]; mi[mm] = {'H': cl0['H'], 'sc': cl0['sc']}
    m_tiles = comb(n-1, 2)
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2: tiles.append((x, y))
    tile_to_P = {}; base_P = set()
    for k, (i, j) in enumerate(P):
        if j == i + 1: base_P.add(k)
        else:
            for t_idx, (tx, ty) in enumerate(tiles):
                if (ty-1, tx-1) == (i, j):
                    tile_to_P[t_idx] = k; break
    total = 2 ** m_tiles; V = mid
    tc = [0] * total
    for tb in range(total):
        full_bits = 0
        for k, (i, j) in enumerate(P):
            if k not in base_P:
                for ti in range(m_tiles):
                    if tile_to_P.get(ti) == k:
                        bit = (tb >> (m_tiles - 1 - ti)) & 1
                        if bit: full_bits |= (1 << (m_total - 1 - k))
                        break
        adj = b2a(full_bits); cid = lk(adj)
        tc[tb] = mn[cid] if cid >= 0 else -1
    return V, mi, tc, m_tiles, total, tiles

print("=" * 72)
print("  WAGGLY AT LARGE n")
print("  opus-2026-03-24-S300")
print("=" * 72)

n = 7
t0 = time.time()
V, mi, tc, m, total, tiles = build_data(n)
print(f"\n  n={n}: V={V}, m={m}, 2^m={total}. Built in {time.time()-t0:.1f}s")

total_pairs = V*(V-1)//2
print(f"  Total class pairs: {total_pairs}")

# Order 1 (wiggly) — fast
t1 = time.time()
E1 = set()
for tb in range(total):
    ma = tc[tb]
    if ma < 0: continue
    for ti in range(m):
        mb = tc[tb ^ (1 << (m-1-ti))]
        if mb >= 0 and ma != mb:
            E1.add((min(ma,mb), max(ma,mb)))

print(f"\n  F_1 (wiggly): {len(E1)}/{total_pairs} ({100*len(E1)/total_pairs:.1f}%) in {time.time()-t1:.1f}s")

# Order 2 — feasible (32K × C(15,2) = 32K × 105 = 3.4M lookups)
t2 = time.time()
E2 = set()
for tb in range(total):
    ma = tc[tb]
    if ma < 0: continue
    for ti in range(m):
        for tj in range(ti+1, m):
            mb = tc[tb ^ (1 << (m-1-ti)) ^ (1 << (m-1-tj))]
            if mb >= 0 and ma != mb:
                e = (min(ma,mb), max(ma,mb))
                if e not in E1:
                    E2.add(e)

cumul_2 = E1 | E2
print(f"  F_2 (+waggly-2): {len(cumul_2)}/{total_pairs} ({100*len(cumul_2)/total_pairs:.1f}%), "
      f"+{len(E2)} new, {time.time()-t2:.1f}s")

# Order 3 — expensive (32K × C(15,3) = 32K × 455 = 14.9M lookups)
t3 = time.time()
E3 = set()
for tb in range(total):
    ma = tc[tb]
    if ma < 0: continue
    for combo in combinations(range(m), 3):
        flip = 0
        for ti in combo: flip ^= (1 << (m-1-ti))
        mb = tc[tb ^ flip]
        if mb >= 0 and ma != mb:
            e = (min(ma,mb), max(ma,mb))
            if e not in cumul_2:
                E3.add(e)
    if tb % 10000 == 0 and tb > 0:
        print(f"    Order 3: {tb}/{total}...", flush=True)

cumul_3 = cumul_2 | E3
gaps_after_3 = total_pairs - len(cumul_3)
print(f"  F_3 (+waggly-3): {len(cumul_3)}/{total_pairs} ({100*len(cumul_3)/total_pairs:.1f}%), "
      f"+{len(E3)} new, {gaps_after_3} gaps, {time.time()-t3:.1f}s")

# Complement (order m)
EB = set()
for tb in range(total):
    ma = tc[tb]
    if ma < 0: continue
    mb = tc[((1 << m) - 1) ^ tb]
    if mb >= 0 and ma != mb:
        EB.add((min(ma,mb), max(ma,mb)))

cumul_3B = cumul_3 | EB
print(f"  F_3+comp: {len(cumul_3B)}/{total_pairs} ({100*len(cumul_3B)/total_pairs:.1f}%)")

# If gaps remain after order 3, try order 4 on the remaining gaps
if gaps_after_3 > 0 and gaps_after_3 <= 500:
    print(f"\n  Checking order 4 for {gaps_after_3} remaining gaps...")
    t4 = time.time()
    E4 = set()
    for tb in range(total):
        ma = tc[tb]
        if ma < 0: continue
        for combo in combinations(range(m), 4):
            flip = 0
            for ti in combo: flip ^= (1 << (m-1-ti))
            mb = tc[tb ^ flip]
            if mb >= 0 and ma != mb:
                e = (min(ma,mb), max(ma,mb))
                if e not in cumul_3:
                    E4.add(e)
        if tb % 10000 == 0 and tb > 0:
            print(f"    Order 4: {tb}/{total}...", flush=True)

    cumul_4 = cumul_3 | E4
    gaps_after_4 = total_pairs - len(cumul_4)
    print(f"  F_4: {len(cumul_4)}/{total_pairs} ({100*len(cumul_4)/total_pairs:.1f}%), "
          f"+{len(E4)} new, {gaps_after_4} gaps, {time.time()-t4:.1f}s")

    if gaps_after_4 == 0:
        print(f"  *** F_4 IS COMPLETE! k* = 4 ***")

# Summary
print(f"\n{'='*72}")
print(f"  SUMMARY n={n}")
print(f"{'='*72}")
print(f"  Diameter (from S296): known to be 7 at n=7")
print(f"  k* (completeness): {'≤4' if gaps_after_3 <= 500 else '?'}")
print(f"  Filling: F_1={len(E1)}, F_2={len(cumul_2)}, F_3={len(cumul_3)}")
print(f"  Gaps: {total_pairs-len(E1)} → {total_pairs-len(cumul_2)} → {gaps_after_3}")

# ΔH by order
print(f"\n  |ΔH| BY ORDER:")
for name, edges in [("Order 1", E1), ("Order 2 new", E2), ("Order 3 new", E3)]:
    if edges:
        dH = [abs(mi[a]['H'] - mi[b]['H']) for (a,b) in edges]
        print(f"    {name}: {len(edges)} edges, mean |ΔH| = {np.mean(dH):.1f}")

print(f"\n  Total time: {time.time()-t0:.1f}s")
print("DONE.")
