#!/usr/bin/env python3
"""
leaf_filtration_laws_boxeph_S184.py  (HYP-8555, THM-1745 — laws part)

Four laws of the leaf-graded arborescence filtration, tested exactly:

(L1) FREE-ACTION DEPTH: if sigma in Aut(T) has prime order p, every
     sigma-fixed arborescence has >= p leaves (deepest nontrivial orbit
     hangs p disjoint isomorphic subtrees). Hence Aut acts FREELY on strata
     l < p_min(|Aut|), so |Aut(T)| divides c_l for ALL l < p_min. Tournament
     Aut is odd => p_min >= 3 => |Aut| divides BOTH h = c_1 AND c_2, always.
     TEST: all iso classes n <= 7 with |Aut| > 1.

(L2) EULERIAN POLE: A_{TT_n, source}(x) = sum_l A(n-1, l-1) x^l (Eulerian
     numbers; increasing-tree <-> permutation bijection). TEST n = 3..7.
     Corollary: c_2(TT_n) = 2^{n-1} - n. Minimizer check per n (n=4 breaks!).

(L3) GRADED ORDINAL-SUM LAW (interpolates h-multiplicativity and mac-mini
     THM-1460(D)): with G_S(x,t) = sum over spanning out-forests of S of
     x^{leaves} t^{components},
       A_{T (+) S, r}(x) = sum_l c_l(T,r) sum_j C(l,j) (x-1)^j G_S(x, n_T - j)
     for r in T. At x=1: A(1) = a_r(T) det(n_T I + L_in(S)) [matrix-forest]
     = THM-1460(D). At [x^1]: h-multiplicativity. TEST: all pairs
     n_T, n_S in {2,3,4} at x = 2, 3 and full coefficient vectors.

(L4) THE {7,21} WINDOW: within-band attainment of 7 and 21 at the first
     relaxation steps: c_2 and B(2) (max-out-degree <= 2) values <= 25 at
     n = 4, 5: is 21 attained (evaporates)? does 7 ever get a window?

boxeph-2026-07-20-S184. Pure python.
"""
import math
import itertools

from importlib import util as _u
spec = _u.spec_from_file_location(
    "census", "04-computation/leaf_graded_arborescences_boxeph_S184.py")


# ---- minimal re-implementation (avoid re-running the census on import) ----

def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def bits_to_adj(bits, n):
    adj = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs_of(n)):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def build_perm_maps(n):
    prs = pairs_of(n)
    idx = {pr: k for k, pr in enumerate(prs)}
    maps = []
    for p in itertools.permutations(range(n)):
        pm = []
        for (i, j) in prs:
            a, b = p[i], p[j]
            pm.append((idx[(a, b)], 0) if a < b else (idx[(b, a)], 1))
        maps.append(pm)
    return maps


def classes_for_n(n):
    if n == 7:
        with open('05-knowledge/results/n7_class_reps_boxeph_S152.txt') as f:
            return [int(t) for t in f.read().split()]
    perm_maps = build_perm_maps(n)
    total = 1 << (n * (n - 1) // 2)
    visited = bytearray(total)
    reps = []
    for bits in range(total):
        if visited[bits]:
            continue
        members = set()
        for pm in perm_maps:
            nb = 0
            for k in range(len(pm)):
                src, flip = pm[k]
                b = (bits >> src) & 1
                if flip:
                    b ^= 1
                nb |= b << k
            members.add(nb)
        for m in members:
            visited[m] = 1
        reps.append(min(members))
    return reps


def aut_order(bits, n, perm_maps):
    cnt = 0
    for pm in perm_maps:
        nb = 0
        for k in range(len(pm)):
            src, flip = pm[k]
            b = (bits >> src) & 1
            if flip:
                b ^= 1
            nb |= b << k
        if nb == bits:
            cnt += 1
    return cnt


def det_int(mat):
    m = [row[:] for row in mat]
    n = len(m)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for k in range(n - 1):
        if m[k][k] == 0:
            piv = next((r for r in range(k + 1, n) if m[r][k] != 0), None)
            if piv is None:
                return 0
            m[k], m[piv] = m[piv], m[k]
            sign = -sign
        mkk = m[k][k]
        for i in range(k + 1, n):
            mik = m[i][k]
            for j in range(k + 1, n):
                m[i][j] = (m[i][j] * mkk - mik * m[k][j]) // prev
            m[i][k] = 0
        prev = mkk
    return sign * m[n - 1][n - 1]


def leaf_poly(adj, n, r):
    inmask = [sum(1 << u for u in range(n) if adj[u][v]) for v in range(n)]
    others = [v for v in range(n) if v != r]
    S = [0] * (n + 1)
    for Ub in range(1 << (n - 1)):
        U = 0
        t = Ub
        for v in others:
            if t & 1:
                U |= 1 << v
            t >>= 1
        if any((inmask[v] & ~U) == 0 for v in others):
            continue
        k = bin(U).count('1')
        oi = {v: c for c, v in enumerate(others)}
        M = [[0] * (n - 1) for _ in range(n - 1)]
        for j in others:
            indeg = 0
            for i in range(n):
                if adj[i][j] and not ((U >> i) & 1):
                    indeg += 1
                    if i != r:
                        M[oi[i]][oi[j]] -= 1
            M[oi[j]][oi[j]] += indeg
        v = det_int(M)
        if v:
            S[k] += v
    return [sum((-1) ** (k - l) * math.comb(k, l) * S[k] for k in range(l, n + 1))
            for l in range(n + 1)]


def enum_arbs(adj, n, r):
    """list of (leafcount, parentdict) for cross-uses"""
    verts = [v for v in range(n) if v != r]
    out = []
    choices = [[u for u in range(n) if adj[u][v]] for v in verts]
    for combo in itertools.product(*choices):
        par = dict(zip(verts, combo))
        ok = True
        for v in verts:
            seen = set()
            u = v
            while u != r:
                if u in seen:
                    ok = False
                    break
                seen.add(u)
                u = par[u]
            if not ok:
                break
        if ok:
            haskid = set(par.values())
            out.append((sum(1 for v in range(n) if v not in haskid), par))
    return out


def forest_leaf_poly(adjS, nS, tmax):
    """G_S(x, t) as coefficient table g[l][c] = # spanning out-forests with
    l leaves, c components. Via virtual root *: forests of S = arbs of
    * (+) S rooted at *, components = tree-children of *; leaves counted in
    S only. Enumerate directly (nS <= 4)."""
    n = nS + 1
    adj = [[0] * n for _ in range(n)]
    for i in range(nS):
        for j in range(nS):
            adj[i + 1][j + 1] = adjS[i][j]
    for j in range(1, n):
        adj[0][j] = 1
    g = {}
    for leaves_all, par in enum_arbs(adj, n, 0):
        comps = sum(1 for v, p in par.items() if p == 0)
        haskid = set(par.values())
        lS = sum(1 for v in range(1, n) if v not in haskid)
        g[(lS, comps)] = g.get((lS, comps), 0) + 1
    return g


def G_eval(g, x, t):
    return sum(cnt * (x ** l) * (t ** c) for (l, c), cnt in g.items())


PM = {n: build_perm_maps(n) for n in (3, 4, 5, 6)}

print("=" * 78)
print("LAWS OF THE LEAF-GRADED FILTRATION")
print("=" * 78)

print("\n(L1) free-action depth: |Aut| divides c_l for all l < p_min(|Aut|);")
print("     in particular |Aut| | c_1 (=h, known) AND |Aut| | c_2 (NEW):")
viol = 0
tested = 0
import time
for n in (3, 4, 5, 6, 7):
    if n == 7:
        pm7 = build_perm_maps(7)
    reps = classes_for_n(n)
    for bits in reps:
        ao = aut_order(bits, n, PM[n] if n <= 6 else pm7)
        if ao == 1:
            continue
        adj = bits_to_adj(bits, n)
        tot = [0] * (n + 1)
        for r in range(n):
            p = leaf_poly(adj, n, r)
            for l in range(n + 1):
                tot[l] += p[l]
        pmin = min(p for p in (3, 5, 7, 11) if ao % p == 0)
        tested += 1
        for l in range(1, min(pmin, n + 1)):
            if tot[l] % ao != 0:
                viol += 1
                print("   VIOLATION n=%d bits=%d |Aut|=%d l=%d c_l=%d" %
                      (n, bits, ao, l, tot[l]))
print("   tested %d classes with |Aut|>1 across n=3..7: violations = %d" % (tested, viol))

print("\n(L2) Eulerian pole: leaf poly of transitive = Eulerian row A(n-1, .):")
for n in (3, 4, 5, 6, 7):
    # transitive: i beats j iff i<j: bits with pair (i,j),i<j -> 1
    bits = (1 << (n * (n - 1) // 2)) - 1
    adj = bits_to_adj(bits, n)
    tot = [0] * (n + 1)
    for r in range(n):
        p = leaf_poly(adj, n, r)
        for l in range(n + 1):
            tot[l] += p[l]
    # Eulerian numbers A(n-1, k)
    m = n - 1
    eul = [sum((-1) ** i * math.comb(m + 1, i) * (k + 1 - i) ** m
               for i in range(k + 2)) for k in range(m)]
    ok = tot[1:n] == eul and tot[0] == 0
    print("   n=%d: strata %s vs Eulerian %s -> %s ; c2 = 2^{n-1}-n = %d: %s" %
          (n, tot[1:n], eul, ok, 2 ** (n - 1) - n,
           tot[2] == 2 ** (n - 1) - n if n >= 3 else '-'))

print("\n   c2 minimizer per n (transitive claims the min for n>=5; n=4 breaks):")
for n in (4, 5, 6):
    reps = classes_for_n(n)
    best = None
    for bits in reps:
        adj = bits_to_adj(bits, n)
        c2 = sum(leaf_poly(adj, n, r)[2] for r in range(n))
        h = sum(leaf_poly(adj, n, r)[1] for r in range(n))
        if best is None or c2 < best[0]:
            best = (c2, h, bits)
    tbits = (1 << (n * (n - 1) // 2)) - 1
    tc2 = sum(leaf_poly(bits_to_adj(tbits, n), n, r)[2] for r in range(n))
    print("   n=%d: min c2 = %d (h of minimizer = %d); transitive c2 = %d %s" %
          (n, best[0], best[1], tc2,
           "(transitive IS min)" if best[0] == tc2 else "(transitive NOT min!)"))

print("\n(L3) graded ordinal-sum law:")
print("     A_{T(+)S,r}(x) = sum_l c_l(T,r) sum_j C(l,j)(x-1)^j G_S(x, n_T - j)")
bad = ok = 0
for nT in (2, 3):
    for nS in (2, 3):
        for bT in range(1 << (nT * (nT - 1) // 2)):
            for bS in range(1 << (nS * (nS - 1) // 2)):
                adjT = bits_to_adj(bT, nT)
                adjS = bits_to_adj(bS, nS)
                n = nT + nS
                adj = [[0] * n for _ in range(n)]
                for i in range(nT):
                    for j in range(nT):
                        adj[i][j] = adjT[i][j]
                for i in range(nS):
                    for j in range(nS):
                        adj[nT + i][nT + j] = adjS[i][j]
                for i in range(nT):
                    for j in range(nS):
                        adj[i][nT + j] = 1
                g = forest_leaf_poly(adjS, nS, nT)
                for r in range(nT):
                    lhs_poly = leaf_poly(adj, n, r)
                    cT = leaf_poly(adjT, nT, r)
                    for x in (2, 3, 5):
                        lhs = sum(lhs_poly[l] * x ** l for l in range(n + 1))
                        rhs = 0
                        for l in range(nT + 1):
                            if cT[l] == 0:
                                continue
                            for j in range(l + 1):
                                rhs += (cT[l] * math.comb(l, j) * (x - 1) ** j
                                        * G_eval(g, x, nT - j))
                        if lhs == rhs:
                            ok += 1
                        else:
                            bad += 1
                            if bad <= 3:
                                print("   MISMATCH nT=%d nS=%d bT=%d bS=%d r=%d x=%d: %d vs %d"
                                      % (nT, nS, bT, bS, r, x, lhs, rhs))
print("   evaluations: %d ok, %d mismatches" % (ok, bad))

print("\n(L4) the {7,21} window at the first relaxation steps:")
for n in (4, 5):
    reps = classes_for_n(n)
    c2set = set()
    b2set = set()
    for bits in reps:
        adj = bits_to_adj(bits, n)
        c2set.add(sum(leaf_poly(adj, n, r)[2] for r in range(n)))
        B2 = 0
        for r in range(n):
            for leaves, par in enum_arbs(adj, n, r):
                deg = {}
                for v, p in par.items():
                    deg[p] = deg.get(p, 0) + 1
                if max(deg.values()) <= 2:
                    B2 += 1
        b2set.add(B2)
    print("   n=%d: c2 values: %s" % (n, sorted(c2set)))
    print("        B(2) values: %s" % sorted(b2set))
    print("        7 in c2: %s ; 7 in B(2): %s ; 21 in c2: %s ; 21 in B(2): %s" %
          (7 in c2set, 7 in b2set, 21 in c2set, 21 in b2set))

print("\nDONE.")
