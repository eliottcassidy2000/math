#!/usr/bin/env python3
"""
tournament_metagraph.py — Library for tournament isomorphism class graph
opus-2026-03-23-S248

Computes Burnside quantities at any n, builds G_n/Z_2 for n ≤ 9 via nauty.

Usage:
  from tournament_metagraph import burnside_quantities, build_metagraph

  # Burnside (any n, no enumeration):
  V, SC, T, edge_orbits = burnside_quantities(n)

  # Full metagraph (n ≤ 9, uses nauty):
  G = build_metagraph(n)
  print(G['E'], G['V_merged'], G['spine_edges'])
"""

from math import comb, factorial, gcd
from collections import Counter, defaultdict
from itertools import permutations
import subprocess

# ============================================================================
# BURNSIDE ENGINE (works at any n, polynomial in partition count)
# ============================================================================

def _partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield []; return
    for f in range(min(n, mx), 0, -1):
        for r in _partitions(n - f, f): yield [f] + r

def _ctc(n, ct):
    c = Counter(ct); r = factorial(n)
    for l, k in c.items(): r //= (l ** k) * factorial(k)
    return r

def _fix_t(ct):
    for c in ct:
        if c % 2 == 0: return 0
    exp = sum((c - 1) // 2 for c in ct)
    for i in range(len(ct)):
        for j in range(i + 1, len(ct)): exp += gcd(ct[i], ct[j])
    return 2 ** exp

def _fix_anti(ct, n):
    f = ct.count(1)
    if f >= 2: return 0
    perm = [0] * n; pos = 0
    for c in ct:
        for i in range(c): perm[pos + i] = pos + (i + 1) % c
        pos += c
    visited = set(); free = 0
    for i in range(n):
        for j in range(n):
            if i == j or (i, j) in visited: continue
            orbit = []; a, b = i, j
            while (a, b) not in visited:
                visited.add((a, b)); orbit.append((a, b))
                a, b = perm[a], perm[b]
            L = len(orbit); orbit_set = set(orbit); rev = (j, i)
            if rev in orbit_set:
                k = orbit.index(rev)
                if k % 2 == 1: free += 1
                else: return 0
            else:
                if rev in visited: pass
                else:
                    if L % 2 == 0: free += 1
                    else: return 0
    return 2 ** free

def burnside_quantities(n):
    """Compute V_n, SC_n, T_n, T_anti, edge_orbits via Burnside."""
    m = comb(n, 2); nf = factorial(n)
    V_sum = 0; T_sum = 0; SC_sum = 0; TA_sum = 0
    for ct in _partitions(n):
        ct_list = list(ct)
        ctc = _ctc(n, ct_list)
        ft = _fix_t(ct_list)
        fa = comb(ct_list.count(1), 2)
        fta = _fix_anti(ct_list, n)
        V_sum += ctc * ft
        T_sum += ctc * ft * fa
        SC_sum += ctc * fta
        TA_sum += ctc * fta * fa
    V = V_sum // nf; T = T_sum // nf
    SC = SC_sum // nf; TA = TA_sum // nf
    EO = T // 2 + factorial(n - 2)
    V_merged = (V + SC) // 2
    T_merged = (T + TA) // 2
    return {
        'n': n, 'm': m, 'V': V, 'SC': SC, 'T': T, 'T_anti': TA,
        'V_merged': V_merged, 'T_merged': T_merged,
        'edge_orbits': EO,
        'E_approx': EO,  # upper bound; exact for n≥12
    }

def derangement(n):
    if n == 0: return 1
    if n == 1: return 0
    d = [1, 0]
    for k in range(2, n + 1): d.append((k - 1) * (d[-1] + d[-2]))
    return d[n]

# ============================================================================
# NAUTY-BASED METAGRAPH (n ≤ 9)
# ============================================================================

def build_metagraph(n, verbose=False):
    """Build the full metagraph G_n/Z_2 using nauty."""
    if n > 9:
        raise ValueError("n > 9 requires too much computation; use burnside_quantities()")
    m = comb(n, 2)
    P = [(i, j) for i in range(n) for j in range(i + 1, n)]

    # Get classes from nauty
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    lines = [l.strip() for l in (r.stdout or '').split('\n')
             if len(l.strip()) == m and all(c in '01' for c in l.strip())]

    def b2a(bits):
        a = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(P):
            if bits & (1 << (m - 1 - k)): a[i][j] = 1
            else: a[j][i] = 1
        return a

    def ss(a): return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))
    def c3(a):
        c = 0
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    if a[i][j] and a[j][k] and a[k][i]: c += 1
                    if a[i][k] and a[k][j] and a[j][i]: c += 1
        return c
    def H(a):
        dp = {}
        for v in range(n): dp[(1 << v, v)] = 1
        for S in range(1, 1 << n):
            for v in range(n):
                if not (S & (1 << v)): continue
                val = dp.get((S, v), 0)
                if val == 0: continue
                for u in range(n):
                    if S & (1 << u): continue
                    if a[v][u]: dp[(S | (1 << u), u)] = dp.get((S | (1 << u), u), 0) + val
        return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))
    def cp(a):
        scores = [sum(a[i][j] for j in range(n)) for i in range(n)]
        sg = defaultdict(list)
        for v in range(n): sg[scores[v]].append(v)
        gs = [sg[s] for s in sorted(set(scores))]
        if all(len(g) == 1 for g in gs):
            p = [g[0] for g in gs]
            return tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
        best = None
        def gen(gs2):
            if not gs2: yield []; return
            for p in permutations(gs2[0]):
                for rest in gen(gs2[1:]): yield list(p) + rest
        for p in gen(gs):
            f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or f < best: best = f
        return best

    cls = []
    hm = defaultdict(list); cm = {}
    for i, line in enumerate(lines):
        bits = int(line, 2); adj = b2a(bits)
        s = ss(adj); c = c3(adj); h = H(adj); canon = cp(adj)
        comp = [[1 - adj[a][b] if a != b else 0 for b in range(n)] for a in range(n)]
        cc = cp(comp); sc = (canon == cc)
        cls.append({'cid': i, 'bits': bits, 'adj': adj, 'score': s,
                    'c3': c, 'H': h, 'sc': sc, 'comp_canon': cc, 'canon': canon})
        hm[(s, c, h)].append(i); cm[canon] = i

    for cl in cls: cl['comp_cid'] = cm.get(cl['comp_canon'], -1)

    def lk(adj):
        s = ss(adj); c = c3(adj); h = H(adj)
        cids = hm.get((s, c, h))
        if not cids: return -1
        if len(cids) == 1: return cids[0]
        return cm.get(cp(adj), -1)

    # Build edges
    edges = set()
    for cl in cls:
        b = cl['bits']
        for arc in range(m):
            fb = b ^ (1 << (m - 1 - arc))
            fa = b2a(fb)
            nb = lk(fa)
            if nb >= 0 and nb != cl['cid']:
                edges.add((min(cl['cid'], nb), max(cl['cid'], nb)))

    # Merged graph
    mn = {}; mid = 0
    for cl in cls:
        ci = cl['cid']
        if ci in mn: continue
        if cl['sc']: mn[ci] = mid; mid += 1
        else:
            mn[ci] = mid; comp = cl['comp_cid']
            if comp >= 0 and comp != ci: mn[comp] = mid
            mid += 1
    me = set()
    for (a, b) in edges:
        x, y = mn[a], mn[b]
        if x != y: me.add((min(x, y), max(x, y)))

    mi = {}
    for mm in range(mid):
        cids = [cl['cid'] for cl in cls if mn[cl['cid']] == mm]
        cl0 = cls[cids[0]]
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'], 'score': cl0['score']}

    # Spine (SC-SC edges)
    spine = set()
    for (a, b) in me:
        if mi[a]['sc'] and mi[b]['sc']: spine.add((a, b))

    return {
        'n': n, 'V': len(cls), 'E': len(edges),
        'V_merged': mid, 'E_merged': len(me),
        'SC': sum(1 for mm in range(mid) if mi[mm]['sc']),
        'spine_edges': len(spine),
        'classes': cls, 'merged_info': mi, 'merged_edges': me,
    }

# ============================================================================
# DEMO
# ============================================================================

if __name__ == '__main__':
    print("TOURNAMENT METAGRAPH LIBRARY")
    print("=" * 60)

    # Burnside quantities for large n
    print("\nBurnside quantities (no enumeration needed):")
    for n in range(3, 16):
        B = burnside_quantities(n)
        print(f"  n={n:2d}: V={B['V']:>14d}  SC={B['SC']:>10d}  "
              f"V_m={B['V_merged']:>14d}  EO={B['edge_orbits']:>16d}")

    # Full metagraph for small n
    print("\nFull metagraph via nauty:")
    for n in range(3, 9):
        G = build_metagraph(n)
        print(f"  n={n}: V={G['V']}, E={G['E']}, V_m={G['V_merged']}, "
              f"E_m={G['E_merged']}, SC={G['SC']}, spine={G['spine_edges']}")
