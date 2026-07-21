#!/usr/bin/env python3
"""
leaf_graded_arborescences_boxeph_S184.py  (HYP-8555, THM-1745)

THE LEAF-GRADED ARBORESCENCE FILTRATION: A_{T,r}(x) = sum_arb x^{#leaves}
(spanning out-arborescences rooted at r; leaf = tree-out-degree 0).
  [x^1] Sum_r A_{T,r} = h  (Hamiltonian paths — the {7,21} stratum,
                            multiplicative under ordinal sum)
  A_{T,r}(1) = a_r         (death-star S71's even-permitting determinant,
                            factorial under ordinal sum)

METHOD (childless inclusion-exclusion, cross-checked by enumeration):
N_r(U) = # arbs rooted at r with every u in U childless = matrix-tree
cofactor after deleting all out-edges of U. S_k = sum_{|U|=k} N_r(U) =
sum_l binom(l,k) c_{r,l}; binomial inversion gives c_{r,l}. N_r(U) = 0 if
r in U (the root needs a child for n >= 2) or if some v != r has all its
in-neighbours inside U (no available parent).

CENSUS: all iso classes n = 3..6 (orbit-marking canonicalization) + all 456
classes at n = 7 (S152 reps). Analyses:
(1) c2 (2-leaf stratum) spectrum per n and pooled — holes?
(2) (h, c2) joint spectrum near the {7,21} boundary.
(3) signed evaluation Sum_r A_r(-1) — Redei-analog hunt.
(5) max-out-degree ladder B(k) at n <= 6: B(1) = h ... B(n-1) = a: at which
    relaxation step do the holes {7,21} evaporate?

boxeph-2026-07-20-S184. Pure python, no deps.
"""
import math
import itertools

# ---------- tournament plumbing ----------

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
            if a < b:
                pm.append((idx[(a, b)], 0))
            else:
                pm.append((idx[(b, a)], 1))
        maps.append(pm)
    return maps


def classes_for_n(n):
    if n == 7:
        reps = []
        with open('05-knowledge/results/n7_class_reps_boxeph_S152.txt') as f:
            for tok in f.read().split():
                reps.append(int(tok))
        return reps
    perm_maps = build_perm_maps(n)
    total = 1 << (n * (n - 1) // 2)
    visited = bytearray(total)
    reps = []
    for bits in range(total):
        if visited[bits]:
            continue
        orbit_min = bits
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
            if nb < orbit_min:
                orbit_min = nb
        for m in members:
            visited[m] = 1
        reps.append(orbit_min)
    return reps


def ham_paths(adjmask, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1 << n):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c or not ((S >> v) & 1):
                continue
            m = adjmask[v] & ~S
            while m:
                w = (m & -m).bit_length() - 1
                dp[S | (1 << w)][w] += c
                m &= m - 1
    return sum(dp[full][v] for v in range(n))


def det_int(mat):
    m = [row[:] for row in mat]
    n = len(m)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        if m[k][k] == 0:
            piv = None
            for r in range(k + 1, n):
                if m[r][k] != 0:
                    piv = r
                    break
            if piv is None:
                return 0
            m[k], m[piv] = m[piv], m[k]
            sign = -sign
        mkk = m[k][k]
        for i in range(k + 1, n):
            mik = m[i][k]
            ri = m[i]
            rk = m[k]
            for j in range(k + 1, n):
                ri[j] = (ri[j] * mkk - mik * rk[j]) // prev
            ri[k] = 0
        prev = mkk
    return sign * m[n - 1][n - 1]


def leaf_poly(adj, n, r, inmask):
    S = [0] * (n + 1)
    others = [v for v in range(n) if v != r]
    for Ubits in range(1 << (n - 1)):
        # embed U into V \ {r}
        U = 0
        t = Ubits
        for v in others:
            if t & 1:
                U |= 1 << v
            t >>= 1
        dead = False
        for v in others:
            if inmask[v] & ~U == 0:
                dead = True
                break
        k = bin(U).count('1')
        if dead:
            continue
        # L_in of the graph with out-edges of U removed, cofactor at r
        M = []
        for i in range(n):
            if i == r:
                continue
            row = []
            for j in range(n):
                if j == r:
                    continue
                row.append(0)
            M.append(row)
        # entries: L[i][j] = indeg_j(if i==j) - A[i][j]
        oi = {v: c for c, v in enumerate(others)}
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
    c = [0] * (n + 1)
    for l in range(n + 1):
        tot = 0
        for k in range(l, n + 1):
            tot += (-1) ** (k - l) * math.comb(k, l) * S[k]
        c[l] = tot
    return c


def enum_arbs_by_leaves(adj, n, r):
    verts = [v for v in range(n) if v != r]
    c = [0] * (n + 1)
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
            leaves = sum(1 for v in range(n) if v not in haskid)
            c[leaves] += 1
    return c


def maxdeg_ladder_total(adj, n):
    Btot = [0] * n
    for r in range(n):
        verts = [v for v in range(n) if v != r]
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
                deg = {}
                for v in verts:
                    deg[par[v]] = deg.get(par[v], 0) + 1
                mx = max(deg.values())
                for k in range(mx, n):
                    Btot[k] += 1
    return Btot


print("=" * 78)
print("LEAF-GRADED ARBORESCENCE CENSUS (leaf = tree-out-degree 0)")
print("=" * 78)

for n in (4, 5):
    reps = classes_for_n(n)
    bad = 0
    for bits in reps:
        adj = bits_to_adj(bits, n)
        inmask = [sum(1 << u for u in range(n) if adj[u][v]) for v in range(n)]
        for r in range(n):
            if leaf_poly(adj, n, r, inmask) != enum_arbs_by_leaves(adj, n, r):
                bad += 1
    print("cross-check n=%d (ALL classes x all roots): mismatches = %d" % (n, bad))

DATA = {}
for n in (3, 4, 5, 6, 7):
    reps = classes_for_n(n)
    rows = []
    for bits in reps:
        adj = bits_to_adj(bits, n)
        adjmask = [sum(1 << j for j in range(n) if adj[i][j]) for i in range(n)]
        inmask = [sum(1 << u for u in range(n) if adj[u][v]) for v in range(n)]
        h = ham_paths(adjmask, n)
        polys = [leaf_poly(adj, n, r, inmask) for r in range(n)]
        tot = [sum(p[l] for p in polys) for l in range(n + 1)]
        a_tot = sum(sum(p) for p in polys)
        s_minus = sum(sum(p[l] * ((-1) ** l) for l in range(n + 1)) for p in polys)
        rows.append({'bits': bits, 'h': h, 'tot': tot, 'a': a_tot, 'am1': s_minus})
        assert tot[1] == h, "bottom!=h at n=%d bits=%d (%d vs %d)" % (n, bits, tot[1], h)
    DATA[n] = rows
    print("n=%d: %d classes; bottom-stratum == h verified on all" % (n, len(rows)))

print("\n(1) c2 spectrum (2-leaf stratum, summed over roots):")
pooled = set()
for n in (4, 5, 6, 7):
    c2s = sorted(set(row['tot'][2] for row in DATA[n]))
    pooled |= set(c2s)
    mx = max(c2s)
    missing = [v for v in range(0, min(mx, 61)) if v not in set(c2s)]
    print("  n=%d: c2 range [%d,%d], #distinct %d; missing (<=60): %s" %
          (n, min(c2s), mx, len(c2s), missing))
print("  POOLED n<=7 missing c2 <= 60: %s" % [v for v in range(61) if v not in pooled])

print("\n(2) (h, c2) joint spectrum near the holes (h <= 11):")
for n in (5, 6, 7):
    pairs = {}
    for row in DATA[n]:
        pairs.setdefault(row['h'], set()).add(row['tot'][2])
    for h in sorted(pairs):
        if h <= 11:
            s = sorted(pairs[h])
            print("  n=%d h=%-3d c2 in %s%s" % (n, h, s[:14], '...' if len(s) > 14 else ''))

print("\n(3) signed evaluation Sum_r A_r(-1):")
for n in (3, 4, 5, 6, 7):
    vals = {}
    for row in DATA[n]:
        vals[row['am1']] = vals.get(row['am1'], 0) + 1
    kv = sorted(vals.items())
    law_h = all(row['am1'] == -row['h'] for row in DATA[n])
    law_h2 = all(row['am1'] == ((-1) ** (n - 1)) * row['h'] for row in DATA[n])
    law_abs = all(abs(row['am1']) == row['h'] for row in DATA[n])
    print("  n=%d: %d distinct values, first %s%s" %
          (n, len(kv), kv[:6], '...' if len(kv) > 6 else ''))
    print("        A(-1)==-h: %s ; A(-1)==(-1)^(n-1) h: %s ; |A(-1)|==h: %s" %
          (law_h, law_h2, law_abs))

print("\n(5) max-out-degree ladder B(k), n<=6 (B(1)=h, B(n-1)=a):")
for n in (5, 6):
    specs = {k: set() for k in range(1, n)}
    for row in DATA[n]:
        adj = bits_to_adj(row['bits'], n)
        Btot = maxdeg_ladder_total(adj, n)
        for k in range(1, n):
            specs[k].add(Btot[k])
    for k in range(1, n):
        s = specs[k]
        mx = max(s)
        holes = [v for v in range(0, min(mx + 1, 41)) if v not in s]
        print("  n=%d B(%d): 7 in spec: %-5s 21 in spec: %-5s missing <=40: %s" %
              (n, k, 7 in s, 21 in s, holes[:18]))

print("\nDONE.")
