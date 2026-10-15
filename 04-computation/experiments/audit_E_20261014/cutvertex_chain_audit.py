#!/usr/bin/env python3
"""Audit E (2026-10-14), THM-4655 statement 5 (added in commit 93231731e6): the cut-vertex transfer chain.
Own implementation of the 4-state chain exactly as stated (states pre, gap, mid, post; u = e_pre + e_gap,
f = e_gap + e_post; block entries from h, h+, h- of the strong components D_1 => ... => D_k of T - z), compared with
explicit DFS enumeration of all Hamiltonian paths, for EVERY strong tournament with n <= 8 (gentourng -z -c) and every
cut vertex z: (a) H(T); (b) #paths with v before z for every v != z; (c) the mid->mid fair-coin identity."""
import subprocess
from itertools import combinations

def d6(s):
    n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63; bits += [(v >> (5 - t)) & 1 for t in range(6)]
    return n, [[bits[i * n + j] for j in range(n)] for i in range(n)]

def dfs_paths(n, A):
    out = []
    def rec(p, used):
        if len(p) == n: out.append(p[:]); return
        for v in range(n):
            if not used >> v & 1 and A[p[-1]][v]:
                p.append(v); rec(p, used | 1 << v); p.pop()
    for s in range(n): rec([s], 1 << s)
    return out

def ham_counts(A, X, z):
    """h, h+, h- of T[X] (X a tuple of vertices) by DFS over its paths."""
    X = list(X)
    if not X: return 1, 1, 1
    h = hp = hm = 0
    def rec(p, used):
        nonlocal h, hp, hm
        if len(p) == len(X):
            h += 1
            if A[z][p[0]]: hp += 1
            if A[p[-1]][z]: hm += 1
            return
        for v in X:
            if v not in used and A[p[-1]][v]:
                p.append(v); used.add(v); rec(p, used); used.discard(v); p.pop()
    for s in X: rec([s], {s})
    return h, hp, hm

def strong_components_order(n, A, verts):
    verts = list(verts)
    R = {u: {u} for u in verts}
    changed = True
    while changed:
        changed = False
        for u in verts:
            new = set(R[u])
            for w in list(R[u]):
                for x in verts:
                    if A[w][x]: new.add(x)
            if new != R[u]: R[u] = new; changed = True
    comps = []; seen = set()
    for u in verts:
        if u in seen: continue
        C = [x for x in verts if x in R[u] and u in R[x]]; seen |= set(C); comps.append(C)
    # order: D_j => D_l (all arcs forward) ; sort by reachability size (a component reaching more comes first)
    comps.sort(key=lambda C: -len(R[C[0]]))
    return comps

def chain(A, comps, z, v=None):
    P, G, M, Q = 0, 1, 2, 3
    vec = [1, 1, 0, 0]                                        # u = e_pre + e_gap
    for D in comps:
        hD, hpD, hmD = ham_counts(A, D, z)
        restrict = (v is not None and v in D)
        sub = [X for k in range(len(D) + 1) for X in combinations(D, k)]
        def hc(X): return ham_counts(A, X, z)
        def comp(X): return tuple(x for x in D if x not in X)
        ok = lambda X: (not restrict) or (v in X)
        pre_mid = sum(hc(X)[0] * hc(comp(X))[1] for X in sub if len(X) < len(D) and ok(X))
        pre_post = sum(hc(X)[2] * hc(comp(X))[1] for X in sub if 0 < len(X) < len(D) and ok(X))
        mid_mid = sum(hc(X)[0] * hc(comp(X))[0] for X in sub if ok(X))
        mid_post = sum(hc(X)[2] * hc(comp(X))[0] for X in sub if len(X) > 0 and ok(X))
        Mj = [[0] * 4 for _ in range(4)]
        Mj[P][P] = hD; Mj[P][G] = hmD; Mj[P][M] = pre_mid; Mj[P][Q] = pre_post
        Mj[G][Q] = 0 if restrict else hpD
        Mj[M][M] = mid_mid; Mj[M][Q] = mid_post
        Mj[Q][Q] = 0 if restrict else hD
        vec = [sum(vec[i] * Mj[i][j] for i in range(4)) for j in range(4)]
    return vec[G] + vec[Q]                                    # f = e_gap + e_post

def fair_coin(A, comps, z):
    for D in comps:
        sub = [X for k in range(len(D) + 1) for X in combinations(D, k)]
        tot = sum(ham_counts(A, X, z)[0] * ham_counts(A, tuple(x for x in D if x not in X), z)[0] for X in sub)
        for v in D:
            part = sum(ham_counts(A, X, z)[0] * ham_counts(A, tuple(x for x in D if x not in X), z)[0] for X in sub if v in X)
            if 2 * part != tot: return False
    return True

for n in range(4, 9):
    lines = subprocess.run(['gentourng', '-q', '-z', '-c', str(n)], capture_output=True, text=True).stdout.split()
    ncut = nH = nV = nF = tours = 0
    for s in lines:
        _, A = d6(s)
        paths = None
        for z in range(n):
            rest = [u for u in range(n) if u != z]
            comps = strong_components_order(n, A, rest)
            if len(comps) < 2: continue
            # check the forward order D_j => D_l
            for j in range(len(comps)):
                for l in range(j + 1, len(comps)):
                    assert all(A[x][y] for x in comps[j] for y in comps[l])
            if paths is None:
                paths = dfs_paths(n, A); tours += 1
            ncut += 1
            Hc = chain(A, comps, z)
            if Hc != len(paths): nH += 1
            for v in rest:
                cnt = sum(1 for p in paths if p.index(v) < p.index(z))
                if chain(A, comps, z, v) != cnt: nV += 1
            if not fair_coin(A, comps, z): nF += 1
    print(f"n={n}: strong tournaments with a cut vertex {tours}; (tournament, cut vertex) pairs {ncut}; "
          f"H mismatches {nH}; v-before-z mismatches {nV}; fair-coin failures {nF}")
