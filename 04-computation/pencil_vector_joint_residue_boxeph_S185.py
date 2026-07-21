#!/usr/bin/env python3
"""
pencil_vector_joint_residue_boxeph_S185.py  (HYP-8560, THM-1760 part 2)

TRANSVERSALITY of the two strongest poly-time fingerprints at n = 7:
  - the PENCIL M_T(t,u,v) = det(tI + uD_in - vA)   (global, graded; S185)
  - klein THM-1750's arborescence VECTOR {a_r}      (per-root, ranking)
plus H (#P) as referee. Questions:
 (1) do the pencil's 11 unresolved groups split under {a_r} (+h)?
 (2) do klein's 4 resistant pairs (same spec A, sum a, {a_r}, H) split
     under the pencil?
 (3) the JOINT residue (pencil, {a_r}): the new deepest poly-time wall;
     and with H added.
boxeph-2026-07-20-S185. Pure python.
"""
import itertools

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


def Mpencil(adj, n, t, u, v):
    indeg = [sum(adj[i][j] for i in range(n)) for j in range(n)]
    M = [[0] * n for _ in range(n)]
    for i in range(n):
        M[i][i] = t + u * indeg[i]
        for j in range(n):
            if i != j and adj[i][j]:
                M[i][j] -= v
    return det_int(M)


def arb_vector(adj, n):
    indeg = [sum(adj[i][j] for i in range(n)) for j in range(n)]
    L = [[0] * n for _ in range(n)]
    for j in range(n):
        L[j][j] = indeg[j]
        for i in range(n):
            if i != j and adj[i][j]:
                L[i][j] -= 1
    out = []
    for r in range(n):
        M = [[L[i][j] for j in range(n) if j != r] for i in range(n) if i != r]
        out.append(det_int(M))
    return out


def ham(adj, n):
    adjmask = [sum(1 << j for j in range(n) if adj[i][j]) for i in range(n)]
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c or not ((S >> v) & 1):
                continue
            m = adjmask[v] & ~S
            while m:
                w = (m & -m).bit_length() - 1
                dp[S | (1 << w)][w] += c
                m &= m - 1
    return sum(dp[full][v] for v in range(n))


n = 7
with open('05-knowledge/results/n7_class_reps_boxeph_S152.txt') as f:
    reps = [int(t) for t in f.read().split()]

pencil_key = {}
vec_key = {}
hval = {}
for bits in reps:
    adj = bits_to_adj(bits, n)
    vec = []
    for t in range(8):
        for u in range(8):
            for v in range(8):
                vec.append(Mpencil(adj, n, t, u, v))
    pencil_key[bits] = tuple(vec)
    av = arb_vector(adj, n)
    vec_key[bits] = tuple(sorted(av))
    hval[bits] = ham(adj, n)

print("=" * 78)
print("PENCIL x VECTOR TRANSVERSALITY, n = 7 (456 classes)")
print("=" * 78)

def groups_by(keyfn):
    g = {}
    for b in reps:
        g.setdefault(keyfn(b), []).append(b)
    return [v for v in g.values() if len(v) > 1]

gp = groups_by(lambda b: pencil_key[b])
gv = groups_by(lambda b: vec_key[b])
gj = groups_by(lambda b: (pencil_key[b], vec_key[b]))
gjh = groups_by(lambda b: (pencil_key[b], vec_key[b], hval[b]))

print("\npencil-only unresolved groups: %d  %s" %
      (len(gp), sorted(len(x) for x in gp)))
print("vector({a_r})-only unresolved groups: %d  %s" %
      (len(gv), sorted(len(x) for x in gv)))
print("JOINT (pencil, vector) unresolved: %d  %s" %
      (len(gj), sorted(len(x) for x in gj)))
print("JOINT + H unresolved: %d  %s" % (len(gjh), sorted(len(x) for x in gjh)))

print("\n(1) pencil's groups split by the vector?")
for g in gp:
    vks = set(vec_key[b] for b in g)
    print("   group %s: %d distinct vectors -> %s" %
          (g, len(vks), "SPLIT" if len(vks) == len(g) else
           ("partial" if len(vks) > 1 else "NOT split")))

print("\n(2) klein's 4 resistant pairs: same (specA, sum a, {a_r}, H) —")
print("    locate as vector+H collisions and test the pencil on them:")
vh = {}
for b in reps:
    vh.setdefault((vec_key[b], hval[b]), []).append(b)
res4 = [(k, v) for k, v in vh.items() if len(v) > 1]
for (vk, h), grp in sorted(res4, key=lambda x: -len(x[1]))[:12]:
    pks = set(pencil_key[b] for b in grp)
    print("   (sum a=%d, H=%d) group size %d: pencil -> %d distinct  %s" %
          (sum(vk), h, len(grp), len(pks),
           "SPLIT" if len(pks) == len(grp) else
           ("partial" if len(pks) > 1 else "NOT split")))

print("\n(3) the JOINT residue in detail (the new deepest poly-time wall):")
for g in gj:
    print("   classes %s: sum a = %d, {a_r} = %s, H = %s" %
          (g, sum(vec_key[g[0]]), vec_key[g[0]], [hval[b] for b in g]))

print("\nDONE.")
