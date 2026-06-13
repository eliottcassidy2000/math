#!/usr/bin/env python3
"""
burnside_wiggly_weights_s284.py — Burnside structure in wiggly weights
opus-2026-03-24-S284

KEY INSIGHT FROM WEB SEARCH: The wiggly weight W[i][j] should relate to
the DOUBLE COSET structure |Aut(T_i) \ S_n / Aut(T_j)|.

ORBITAL FORMULA: rank = (1/|G|) Σ fix(g)²
counts the number of orbits on ordered pairs.

For wiggly weights: W[i][j] counts labeled (tiling, tile) pairs.
Can we decompose W[i][j] in terms of Burnside invariants?

CHECK:
1. Is W[i][j] always divisible by some function of |Aut(i)|, |Aut(j)|?
2. Does W[i][j] / (fiber[i] × fiber[j]) follow a pattern?
3. Can we predict W from H values and |Aut| alone?
"""

import sys, time, subprocess
from math import comb, factorial, gcd as math_gcd
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_and_compute(n):
    """Full computation of W matrix."""
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
    def aut_size(a):
        return sum(1 for p in permutations(range(n))
                   if all(a[p[i]][p[j]] == a[i][j] for i in range(n) for j in range(n)))
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
        aut = aut_size(adj) if n <= 6 else 0
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid':i, 'score':s, 'c3':cc, 'H':h, 'sc':sc,
                    'comp_canon':cc2, 'canon':canon, 'aut':aut})
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
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'], 'aut': cl0['aut']}

    # Compute W from tiling enumeration
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

    total = 2 ** m_tiles
    V = mid
    W = np.zeros((V, V), dtype=int)
    fiber = Counter()

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
        if cid < 0: continue
        ma = mn[cid]; fiber[ma] += 1
        for ti in range(m_tiles):
            flipped = tb ^ (1 << (m_tiles - 1 - ti))
            full_bits_f = 0
            for k, (i, j) in enumerate(P):
                if k not in base_P:
                    for tj in range(m_tiles):
                        if tile_to_P.get(tj) == k:
                            bit = (flipped >> (m_tiles - 1 - tj)) & 1
                            if bit: full_bits_f |= (1 << (m_total - 1 - k))
                            break
            adj_f = b2a(full_bits_f); cid_f = lk(adj_f)
            if cid_f >= 0: W[ma][mn[cid_f]] += 1

    return V, mi, W, dict(fiber), m_tiles

def analyze(n):
    print(f"\n{'#'*72}")
    print(f"  n = {n}: BURNSIDE STRUCTURE IN WIGGLY WEIGHTS")
    print(f"{'#'*72}")

    t0 = time.time()
    V, mi, W, fiber, m = build_and_compute(n)
    print(f"  V={V}, m={m}, computed in {time.time()-t0:.1f}s")

    # 1. DIVISIBILITY: is W[i][j] divisible by gcd(|Aut(i)|, |Aut(j)|)?
    if n <= 6:
        print(f"\n  DIVISIBILITY CHECKS:")
        div_by_gcd_aut = 0; total_pairs = 0
        div_by_fiber_gcd = 0
        for i in range(V):
            for j in range(i, V):
                if W[i][j] > 0:
                    total_pairs += 1
                    g = math_gcd(mi[i]['aut'], mi[j]['aut'])
                    if W[i][j] % g == 0: div_by_gcd_aut += 1
                    fg = math_gcd(fiber[i], fiber[j])
                    if W[i][j] % fg == 0: div_by_fiber_gcd += 1

        print(f"  W[i][j] divisible by gcd(|Aut_i|, |Aut_j|): "
              f"{div_by_gcd_aut}/{total_pairs}")
        print(f"  W[i][j] divisible by gcd(fiber_i, fiber_j): "
              f"{div_by_fiber_gcd}/{total_pairs}")

    # 2. NORMALIZED WEIGHTS: W[i][j] / (fiber[i] × fiber[j])
    print(f"\n  NORMALIZED W[i][j] / (fiber_i × fiber_j):")
    nw_values = []
    for i in range(V):
        for j in range(i+1, V):
            if W[i][j] > 0 and fiber[i] > 0 and fiber[j] > 0:
                nw = W[i][j] / (fiber[i] * fiber[j])
                nw_values.append((i, j, nw))

    if nw_values:
        nws = [x[2] for x in nw_values]
        print(f"  Range: [{min(nws):.6f}, {max(nws):.6f}]")
        print(f"  Mean: {np.mean(nws):.6f}, Std: {np.std(nws):.6f}")

        # Are they rational with small denominators?
        print(f"\n  TOP 10 NORMALIZED WEIGHTS (looking for rational patterns):")
        for i, j, nw in sorted(nw_values, key=lambda x: -x[2])[:10]:
            # Try to identify as p/q with small q
            for q in range(1, 20):
                p = round(nw * q)
                if abs(nw - p/q) < 0.001:
                    print(f"    W[{i},{j}] = {W[i][j]:>6}, "
                          f"fiber=({fiber[i]},{fiber[j]}), "
                          f"norm={nw:.6f} ≈ {p}/{q}, "
                          f"H=({mi[i]['H']},{mi[j]['H']})")
                    break

    # 3. CAN W BE PREDICTED FROM (H, |Aut|, SC/NS) ALONE?
    if n <= 6:
        print(f"\n  PREDICTION TEST: Does (H_i, H_j, |Aut_i|, |Aut_j|) determine W[i][j]?")
        # Group edges by (H_i, H_j, aut_i, aut_j) and check if W is constant
        feature_groups = defaultdict(list)
        for i in range(V):
            for j in range(i+1, V):
                if W[i][j] > 0:
                    key = (mi[i]['H'], mi[j]['H'],
                           mi[i]['aut'], mi[j]['aut'],
                           mi[i]['sc'], mi[j]['sc'])
                    feature_groups[key].append(W[i][j])

        n_groups = len(feature_groups)
        constant_groups = sum(1 for vs in feature_groups.values()
                             if len(set(vs)) == 1)
        print(f"  Feature groups: {n_groups}")
        print(f"  Groups with constant W: {constant_groups}/{n_groups}")

        # Show variable groups
        variable = [(k, vs) for k, vs in feature_groups.items() if len(set(vs)) > 1]
        if variable:
            print(f"\n  VARIABLE GROUPS (same features, different W):")
            for key, values in variable[:5]:
                print(f"    H=({key[0]},{key[1]}), |Aut|=({key[2]},{key[3]}), "
                      f"SC=({key[4]},{key[5]}): W values = {sorted(set(values))}")

    # 4. THE RANK/ORBITAL FORMULA
    # rank = (1/|G|) Σ fix(g)² tells us the number of orbits on pairs
    # For S_n on tournaments: fix(σ) = 2^{po(σ)} for odd-cycle σ, 0 otherwise
    nf = factorial(n)

    from math import gcd as mgcd
    def partitions(nn, mx=None):
        if mx is None: mx = nn
        if nn == 0: yield []; return
        for f in range(min(nn, mx), 0, -1):
            for r in partitions(nn - f, f): yield [f] + r

    def ccs(nn, ct):
        c = Counter(ct); r = factorial(nn)
        for l, k in c.items(): r //= pow(l, k) * factorial(k)
        return r

    def po(ct):
        s = sum((c-1)//2 for c in ct)
        for i in range(len(ct)):
            for j in range(i+1, len(ct)):
                s += mgcd(ct[i], ct[j])
        return s

    rank_sum = 0
    for ct in partitions(n):
        ct_list = list(ct)
        if any(c % 2 == 0 for c in ct_list): continue
        fix = pow(2, po(ct_list))
        cc = ccs(n, ct_list)
        rank_sum += cc * fix * fix  # fix²

    rank = rank_sum // nf
    print(f"\n  ORBITAL FORMULA:")
    print(f"  rank = (1/n!) Σ fix(σ)² = {rank}")
    print(f"  Number of orbits on ORDERED pairs of tournaments = {rank}")
    print(f"  Number of orbits on UNORDERED pairs ≈ {(rank + V) // 2}")
    print(f"  Actual metagraph edges (from W): "
          f"{sum(1 for i in range(V) for j in range(i+1, V) if W[i][j] > 0)}")

    # 5. TRANSITION ORBITS from Burnside
    # T_n = (1/n!) Σ_{odd σ} ccs(σ) × Fix(σ) × C(f(σ), 2)
    # where f(σ) = number of fixed points
    T_sum = 0
    for ct in partitions(n):
        ct_list = list(ct)
        if any(c % 2 == 0 for c in ct_list): continue
        fix = pow(2, po(ct_list))
        cc = ccs(n, ct_list)
        f = ct_list.count(1)
        T_sum += cc * fix * comb(f, 2) if f >= 2 else 0

    T_n = T_sum // nf
    print(f"\n  TRANSITION ORBITS T_n = {T_n}")
    print(f"  W total (sum all entries) = {np.sum(W)}")
    print(f"  Expected: 2^m × m = {2**m * m}")
    print(f"  Match: {np.sum(W) == 2**m * m}")

    return V, mi, W, fiber, m

if __name__ == '__main__':
    print("=" * 72)
    print("  BURNSIDE STRUCTURE IN WIGGLY WEIGHTS")
    print("  opus-2026-03-24-S284")
    print("=" * 72)

    for n in [4, 5, 6]:
        analyze(n)

    print(f"\n{'='*72}")
    print("  KEY FINDINGS")
    print(f"{'='*72}")
    print("""
  The ORBITAL FORMULA rank = (1/n!) Σ fix(σ)² counts orbits
  on ORDERED pairs of tournaments. This gives an UPPER BOUND
  on the number of distinct W[i][j] values.

  The NORMALIZED WEIGHT W[i][j]/(fiber_i × fiber_j) is the
  transition probability between classes in the tiling model.
  It should relate to the DOUBLE COSET structure of S_n.

  DIVISIBILITY: W[i][j] is always divisible by gcd(fiber_i, fiber_j)
  (if verified). This would mean the transition probability is
  always a ratio with denominator dividing lcm(fiber_i, fiber_j).
""")

    print("DONE.")
