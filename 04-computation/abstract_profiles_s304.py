#!/usr/bin/env python3
"""
abstract_profiles_s304.py — Distance profiles applied abstractly
opus-2026-03-24-S304

Apply the distance profile concept creatively:
1. Profile from EVERY class (not just transitive)
2. Distance-regular check (constant intersection numbers?)
3. Cross-profiles: embed classes in low-dim via landmarks
4. Profile generating functions and their properties
5. Weighted profiles (by H, fiber, |Aut|)
6. The "dual profile": given H, what's the distance distribution?
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, deque, Counter
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
    m = comb(n-1, 2)
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
    total = 2 ** m; V = mid
    tc = [0] * total
    for tb in range(total):
        full_bits = 0
        for k, (i, j) in enumerate(P):
            if k not in base_P:
                for ti in range(m):
                    if tile_to_P.get(ti) == k:
                        bit = (tb >> (m - 1 - ti)) & 1
                        if bit: full_bits |= (1 << (m_total - 1 - k))
                        break
        adj = b2a(full_bits); cid = lk(adj)
        tc[tb] = mn[cid] if cid >= 0 else -1
    fiber = Counter()
    for tb in range(total):
        if tc[tb] >= 0: fiber[tc[tb]] += 1
    return V, mi, tc, m, total, dict(fiber)

def full_analysis(n):
    print(f"\n{'#'*72}")
    print(f"  n = {n}: ABSTRACT DISTANCE PROFILES")
    print(f"{'#'*72}")

    t0 = time.time()
    V, mi, tc, m, total, fiber = build_data(n)
    print(f"  V={V}, m={m}. Built in {time.time()-t0:.1f}s")

    # Build adjacency and ALL-PAIRS distances
    adj = defaultdict(set)
    for tb in range(total):
        ma = tc[tb]
        if ma < 0: continue
        for ti in range(m):
            mb = tc[tb ^ (1 << (m-1-ti))]
            if mb >= 0 and ma != mb:
                adj[ma].add(mb); adj[mb].add(ma)

    dist = np.full((V, V), -1, dtype=int)
    for start in range(V):
        d = [-1]*V; d[start] = 0
        q = deque([start])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if d[v] < 0: d[v] = d[u]+1; q.append(v)
        for j in range(V): dist[start][j] = d[j]

    diameter = int(np.max(dist))

    # ============================================================
    # 1. ALL PROFILES: one per class
    # ============================================================
    print(f"\n  1. PROFILE FROM EVERY CLASS:")
    profiles = {}
    for start in range(V):
        p = [0]*(diameter+1)
        for j in range(V):
            if dist[start][j] >= 0: p[dist[start][j]] += 1
        profiles[start] = tuple(p)

    # How many DISTINCT profiles?
    distinct = Counter(profiles.values())
    print(f"    Distinct profiles: {len(distinct)} out of {V}")
    for prof, count in distinct.most_common(5):
        classes_with = [c for c in range(V) if profiles[c] == prof]
        h_vals = [mi[c]['H'] for c in classes_with[:5]]
        print(f"      {prof}: {count} classes, H={h_vals}")

    # ============================================================
    # 2. ECCENTRICITY and CENTER
    # ============================================================
    print(f"\n  2. ECCENTRICITY:")
    eccs = [max(dist[i][j] for j in range(V) if dist[i][j] >= 0) for i in range(V)]
    ecc_dist = Counter(eccs)
    radius = min(eccs)
    center = [i for i in range(V) if eccs[i] == radius]

    print(f"    Diameter: {diameter}, Radius: {radius}")
    print(f"    Eccentricity distribution: {sorted(ecc_dist.items())}")
    print(f"    Center ({len(center)} nodes): H = {sorted([mi[c]['H'] for c in center])[:10]}")

    # Peripheral nodes (eccentricity = diameter)
    periph = [i for i in range(V) if eccs[i] == diameter]
    print(f"    Peripheral ({len(periph)} nodes): H = {sorted([mi[c]['H'] for c in periph])[:10]}")

    # ============================================================
    # 3. DISTANCE-REGULAR CHECK
    # ============================================================
    print(f"\n  3. DISTANCE-REGULAR CHECK:")
    # For distance d, check: for all u,v with dist(u,v)=d,
    # is #{w: dist(u,w)=d-1, dist(v,w)=1} constant?
    is_dist_regular = True
    for d in range(1, diameter+1):
        values = set()
        for u in range(V):
            for v in range(V):
                if dist[u][v] != d: continue
                count = sum(1 for w in range(V)
                           if dist[u][w] == d-1 and dist[v][w] == 1)
                values.add(count)
        if len(values) > 1:
            is_dist_regular = False
            print(f"    d={d}: intersection number c_d = {sorted(values)} (NOT constant)")
        else:
            print(f"    d={d}: c_d = {values.pop()} (constant ✓)")

    print(f"    Distance-regular: {is_dist_regular}")

    # ============================================================
    # 4. CROSS-PROFILE EMBEDDING
    # ============================================================
    print(f"\n  4. CROSS-PROFILE EMBEDDING (landmarks: transitive, max-H):")
    trans_id = 0
    max_H = max(mi[i]['H'] for i in range(V))
    reg_ids = [i for i in range(V) if mi[i]['H'] == max_H]
    reg_id = reg_ids[0]

    # 2D embedding: (dist from transitive, dist from regular)
    coords = [(dist[trans_id][i], dist[reg_id][i]) for i in range(V)]
    print(f"    Coordinates (d_trans, d_reg):")
    coord_groups = defaultdict(list)
    for i in range(V):
        coord_groups[coords[i]].append(i)

    for (dt, dr) in sorted(coord_groups.keys()):
        classes = coord_groups[(dt, dr)]
        h_vals = sorted([mi[c]['H'] for c in classes])
        print(f"      ({dt},{dr}): {len(classes)} classes, H={h_vals[:6]}{'...' if len(h_vals)>6 else ''}")

    # Is the embedding injective? (each class gets unique coords)
    n_coords = len(set(coords))
    print(f"\n    Injective with 2 landmarks: {n_coords == V} ({n_coords}/{V})")

    # If not injective, how many landmarks needed?
    if n_coords < V:
        # Add a third landmark: the class at the "center"
        center_id = center[0]
        coords3 = [(dist[trans_id][i], dist[reg_id][i], dist[center_id][i]) for i in range(V)]
        n_coords3 = len(set(coords3))
        print(f"    Injective with 3 landmarks: {n_coords3 == V} ({n_coords3}/{V})")

    # ============================================================
    # 5. PROFILE GENERATING FUNCTION ANALYSIS
    # ============================================================
    print(f"\n  5. PROFILE GENERATING FUNCTIONS:")

    # For each class, compute P_C(x) = Σ (# at distance d) × x^d
    # Then: Σ_C P_C(x) = V × P_avg(x) where P_avg counts average profile

    # Average profile
    avg_profile = np.zeros(diameter+1)
    for c in range(V):
        for d in range(diameter+1):
            avg_profile[d] += profiles[c][d]
    avg_profile /= V

    print(f"    Average profile: [{', '.join(f'{x:.2f}' for x in avg_profile)}]")
    print(f"    Sum: {sum(avg_profile):.2f} (should be {V})")

    # P_avg(-1)
    p_avg_m1 = sum((-1)**d * avg_profile[d] for d in range(diameter+1))
    print(f"    P_avg(-1) = {p_avg_m1:.4f}")

    # ============================================================
    # 6. H-WEIGHTED DISTANCE MATRIX
    # ============================================================
    print(f"\n  6. H-WEIGHTED ANALYSIS:")

    H_vals = np.array([mi[i]['H'] for i in range(V)], dtype=float)

    # Mean H at each distance from each class
    # Average over all source classes: what's the mean H of classes at distance d?
    H_by_d = np.zeros(diameter+1)
    count_by_d = np.zeros(diameter+1)
    for i in range(V):
        for j in range(V):
            d = dist[i][j]
            if d >= 0:
                H_by_d[d] += mi[j]['H']
                count_by_d[d] += 1

    for d in range(diameter+1):
        if count_by_d[d] > 0:
            print(f"    d={d}: mean H of neighbors = {H_by_d[d]/count_by_d[d]:.2f} "
                  f"(from {int(count_by_d[d]/V)} classes per source)")

    # H correlation matrix: corr(H_i, dist(i,j)) for each j
    dist_H_corr = np.corrcoef(H_vals, dist[0])[0,1] if V > 2 else 0
    print(f"\n    Corr(H, dist_from_transitive) = {dist_H_corr:.4f}")
    dist_reg_corr = np.corrcoef(H_vals, dist[reg_id])[0,1] if V > 2 else 0
    print(f"    Corr(H, dist_from_regular) = {dist_reg_corr:.4f}")

    # ============================================================
    # 7. THE SPECTRAL PROFILE
    # ============================================================
    print(f"\n  7. SPECTRAL PROFILE:")

    # The distance matrix D has eigenvalues related to the adjacency spectrum
    D = dist.astype(float)
    D[D < 0] = 0  # shouldn't happen for connected graph
    D_evals = sorted(np.linalg.eigvalsh(D), reverse=True)

    print(f"    Distance matrix eigenvalues (top 5): "
          f"[{', '.join(f'{e:.2f}' for e in D_evals[:5])}]")
    print(f"    Distance matrix eigenvalues (bottom 3): "
          f"[{', '.join(f'{e:.2f}' for e in D_evals[-3:])}]")

    # Ratio of top two
    if D_evals[1] != 0:
        print(f"    λ_0 / λ_1 = {D_evals[0]/D_evals[1]:.4f}")

    return diameter, profiles, dist

if __name__ == '__main__':
    print("=" * 72)
    print("  ABSTRACT DISTANCE PROFILES")
    print("  opus-2026-03-24-S304")
    print("=" * 72)

    for n in [5, 6, 7]:
        full_analysis(n)

    print("\nDONE.")
