#!/usr/bin/env python3
"""
distance_profile_s303.py — Distance profiles and palindromicity
opus-2026-03-24-S303

The DISTANCE PROFILE from the transitive class:
  P(d) = #{merged classes at metagraph distance d from transitive}

At n=5 this is [1,4,4,1] — PALINDROMIC.

Questions:
1. Is the profile palindromic at all n?
2. Is it the h-vector of some polytope?
3. Does it come from a representation of S_n?
4. What is the generating function?
5. Connection to the Molien series of S_n on Q_m?
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
    return V, mi, tc, m_tiles, total

print("=" * 72)
print("  DISTANCE PROFILES AND PALINDROMICITY")
print("  opus-2026-03-24-S303")
print("=" * 72)

all_profiles = {}

for n in [3, 4, 5, 6, 7]:
    print(f"\n{'#'*72}")
    print(f"  n = {n}")
    print(f"{'#'*72}")

    t0 = time.time()
    V, mi, tc, m, total = build_data(n)
    print(f"  V={V}, m={m}. Built in {time.time()-t0:.1f}s")

    # Build adjacency and compute BFS distances from transitive (node 0)
    adj_list = defaultdict(set)
    for tb in range(total):
        ma = tc[tb]
        if ma < 0: continue
        for ti in range(m):
            mb = tc[tb ^ (1 << (m-1-ti))]
            if mb >= 0 and ma != mb:
                adj_list[ma].add(mb); adj_list[mb].add(ma)

    # BFS from transitive (class 0)
    trans_id = 0  # transitive = all-zeros tiling
    dist = [-1] * V; dist[trans_id] = 0
    queue = deque([trans_id])
    while queue:
        u = queue.popleft()
        for v in adj_list[u]:
            if dist[v] < 0:
                dist[v] = dist[u] + 1
                queue.append(v)

    diameter = max(d for d in dist if d >= 0)

    # Distance profile
    profile = [0] * (diameter + 1)
    for d in dist:
        if d >= 0: profile[d] += 1

    # Check palindromicity
    is_palindrome = (profile == profile[::-1])

    print(f"  Diameter: {diameter}")
    print(f"  Distance profile from transitive: {profile}")
    print(f"  Sum: {sum(profile)} (should be {V})")
    print(f"  PALINDROMIC: {is_palindrome}")

    all_profiles[n] = profile

    # Generating function P(x) = Σ profile[d] × x^d
    print(f"  GF: P(x) = ", end="")
    terms = []
    for d, c in enumerate(profile):
        if c > 0:
            if d == 0: terms.append(str(c))
            elif d == 1: terms.append(f"{c}x")
            else: terms.append(f"{c}x^{d}")
    print(" + ".join(terms))

    # P(1) = V, P(-1) = ?
    p_minus1 = sum((-1)**d * profile[d] for d in range(len(profile)))
    print(f"  P(1) = {sum(profile)} = V")
    print(f"  P(-1) = {p_minus1}")

    # Weighted by fiber
    fiber = Counter()
    for tb in range(total):
        if tc[tb] >= 0: fiber[tc[tb]] += 1

    fiber_profile = [0] * (diameter + 1)
    for mm in range(V):
        if dist[mm] >= 0:
            fiber_profile[dist[mm]] += fiber[mm]

    print(f"  Fiber-weighted profile: {fiber_profile}")
    print(f"  Sum: {sum(fiber_profile)} (should be {total})")

    # Is fiber-weighted profile palindromic?
    fp_palindrome = (fiber_profile == fiber_profile[::-1])
    print(f"  Fiber-weighted PALINDROMIC: {fp_palindrome}")

    # H at each distance level
    print(f"\n  H DISTRIBUTION BY DISTANCE:")
    for d in range(diameter + 1):
        classes_at_d = [mm for mm in range(V) if dist[mm] == d]
        h_vals = sorted([mi[mm]['H'] for mm in classes_at_d])
        sc_count = sum(1 for mm in classes_at_d if mi[mm]['sc'])
        fib_sum = sum(fiber[mm] for mm in classes_at_d)
        mean_h = np.mean(h_vals) if h_vals else 0
        print(f"    d={d}: {len(classes_at_d)} classes, SC={sc_count}, "
              f"H={h_vals[:8]}{'...' if len(h_vals)>8 else ''}, "
              f"mean_H={mean_h:.1f}, fiber_sum={fib_sum}")

    # Distance profile from REGULAR tournament (H=max)
    max_H = max(mi[mm]['H'] for mm in range(V))
    reg_ids = [mm for mm in range(V) if mi[mm]['H'] == max_H]

    if len(reg_ids) == 1:
        reg_id = reg_ids[0]
        dist_reg = [-1] * V; dist_reg[reg_id] = 0
        queue = deque([reg_id])
        while queue:
            u = queue.popleft()
            for v in adj_list[u]:
                if dist_reg[v] < 0:
                    dist_reg[v] = dist_reg[u] + 1
                    queue.append(v)

        reg_profile = [0] * (max(d for d in dist_reg if d >= 0) + 1)
        for d in dist_reg:
            if d >= 0: reg_profile[d] += 1

        print(f"\n  Profile from REGULAR (H={max_H}): {reg_profile}")
        reg_palindrome = (reg_profile == reg_profile[::-1])
        print(f"  PALINDROMIC: {reg_palindrome}")

# Summary
print(f"\n{'='*72}")
print("  SUMMARY: DISTANCE PROFILES")
print(f"{'='*72}")
print(f"  {'n':>3} {'diam':>5} {'profile':>40} {'palindrome':>12}")
for n in sorted(all_profiles.keys()):
    p = all_profiles[n]
    is_pal = (p == p[::-1])
    print(f"  {n:>3} {len(p)-1:>5} {str(p):>40} {'✓' if is_pal else '✗':>12}")

print(f"\n{'='*72}")
print("  REPRESENTATION-THEORETIC INTERPRETATION")
print(f"{'='*72}")
print("""
  The distance profile P(d) from the transitive class is:
    P(x) = Σ_{d=0}^{diam} #{classes at distance d} × x^d

  PALINDROMICITY means: P(d) = P(diam - d) for all d.
  This is equivalent to: x^diam × P(1/x) = P(x)
  = the polynomial is SELF-RECIPROCAL.

  In representation theory, self-reciprocal polynomials arise from:

  1. GORENSTEIN RINGS: the h-vector of a Gorenstein simplicial complex
     is palindromic (Dehn-Sommerville relations). If the metagraph's
     independence complex is Gorenstein, the profile is palindromic.

  2. POINCARÉ DUALITY: if the quotient Q_m / S_n has a manifold-like
     structure with Poincaré duality, Betti numbers are palindromic.
     The distance profile would be the "Betti numbers" of the
     distance filtration.

  3. MOLIEN SERIES: the cycle index of S_n on Q_m gives a polynomial
     in which the coefficient of x^d counts orbits of Hamming weight d.
     This IS the distance profile (if the transitive = weight 0).
     Molien series of real representations are palindromic when the
     group contains -1 (the complement involution).

  4. SELF-DUALITY OF THE TILING SPACE: The complement involution
     (flip all tiles) maps weight-d tilings to weight-(m-d) tilings.
     If this involution preserves iso classes, it would give
     palindromicity of the CLASS profile.
""")

print("DONE.")
