#!/usr/bin/env python3
"""procgen_tbij_20261001_orchestrator_check.py -- the orchestrator's independent audit of the tbij lane (2026-10-01).

Written from the statements in the lane note; the lane's code was not read.
Checks:
  A. Natural-bijection obstruction (Lemma N1 + Theorem O2) at n = 5, 6: compute every labelled switching class's
     stabiliser in S_n, the untwisted ("even") Euler graphs fixed by it, the stabiliser-containment bipartite graph
     between iso types, and its maximum matching: n = 5 -> 2 classes, maximum matching 1; n = 6 -> 6 classes,
     maximum matching 5. So no natural map CLASS -> EVENE is bijective on types.
  B. Theorem O1: no labelled switching class is fixed by all of S_n (n = 3..6).
  C. The n = 5 hand proof: the class stabilisers contain elements of order 5 resp. 3; the Euler graphs on 5 vertices
     with an automorphism of order 3 or 5, and which of them are even.
  D. P2: every bipartite Euler graph on n <= 6 vertices is even (untwisted).
  E. G0 (Royle et al. 2023): #even graphs = #tournaments (A000568) for n <= 5, by brute force over all labelled graphs.
  F. P3 (Higashitani-Ueyama, arXiv:2409.10904 Ex. 4.10): s_{l,n} = t_{l,n} by brute force (union-find orbits) for
     l = 3, 4 and n <= 4, with t_{4,n} = 1, 1, 3, 8 and t_{3,n} = 1, 1, 2, 4.
  G. P1: for odd n <= 7 every switching class has exactly one member with out = in (mod 4) at every vertex.
"""
import itertools
import math
import sys
import time
from collections import defaultdict

OKS = []


def ok(cond, msg):
    OKS.append(bool(cond))
    print(('[OK] ' if cond else '[FAIL] ') + msg, flush=True)


def pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


# ---- tournaments: dict (i,j) i<j -> True if i -> j
def tiling_rep(T, n):
    inL = [0] * n
    for i in range(n - 1):
        inL[i + 1] = inL[i] ^ (1 if T[(i, i + 1)] else 0)
    return tuple(T[(i, j)] ^ bool(inL[i] ^ inL[j]) for (i, j) in pairs(n))


def act_tour(T, g):
    R = {}
    for (i, j), fwd in T.items():
        a, b = (g[i], g[j]) if fwd else (g[j], g[i])
        if a < b:
            R[(a, b)] = True
        else:
            R[(b, a)] = False
    return R


def tiling_reps(n):
    P = pairs(n)
    base = [k for k, (i, j) in enumerate(P) if j == i + 1]
    for m in range(1 << len(P)):
        if any(m >> k & 1 for k in base):
            continue
        yield tuple(bool(m >> k & 1) for k in range(len(P)))


def euler_graphs(n):
    P = pairs(n)
    for m in range(1 << len(P)):
        deg = [0] * n
        for k, (i, j) in enumerate(P):
            if m >> k & 1:
                deg[i] ^= 1
                deg[j] ^= 1
        if not any(deg):
            yield frozenset(P[k] for k in range(len(P)) if m >> k & 1)


def act_graph(F, g):
    return frozenset((min(g[i], g[j]), max(g[i], g[j])) for (i, j) in F)


def sgn(F, g):
    s = 0
    for (i, j) in F:
        if g[i] > g[j]:
            s ^= 1
    return -1 if s else 1


def graph_canon(F, perms):
    return min(tuple(sorted(act_graph(F, g))) for g in perms)


def order(g):
    n = len(g)
    seen = [False] * n
    o = 1
    for x in range(n):
        if not seen[x]:
            L = 0
            y = x
            while not seen[y]:
                seen[y] = True
                y = g[y]
                L += 1
            o = o * L // math.gcd(o, L)
    return o


def class_data(n):
    P = pairs(n)
    perms = list(itertools.permutations(range(n)))
    reps = list(tiling_reps(n))
    idx = {t: k for k, t in enumerate(reps)}
    # image of every labelled class under every permutation
    img = [[0] * len(perms) for _ in reps]
    for k, t in enumerate(reps):
        T = dict(zip(P, t))
        for gi, g in enumerate(perms):
            img[k][gi] = idx[tiling_rep(act_tour(T, g), n)]
    # orbits and stabilisers
    orbit_of = [-1] * len(reps)
    orbits = []
    for k in range(len(reps)):
        if orbit_of[k] >= 0:
            continue
        o = sorted(set(img[k]))
        for x in o:
            orbit_of[x] = len(orbits)
        orbits.append((k, o))
    stab = {k: [perms[gi] for gi in range(len(perms)) if img[k][gi] == k] for k, _ in orbits}
    return perms, reps, orbits, stab, img


def euler_data(n, perms):
    E = list(euler_graphs(n))
    canon = {}
    info = {}
    for F in E:
        c = graph_canon(F, perms)
        canon[F] = c
        if c not in info:
            aut = [g for g in perms if act_graph(F, g) == F]
            info[c] = all(sgn(F, g) == 1 for g in aut)   # untwisted?
    return E, canon, info


def max_matching(left, adj):
    match = {}
    def aug(u, seen):
        for v in adj[u]:
            if v in seen:
                continue
            seen.add(v)
            if v not in match or aug(match[v], seen):
                match[v] = u
                return True
        return False
    return sum(1 for u in left if aug(u, set()))


def check_A_B_C():
    for n in (3, 4, 5, 6):
        perms, reps, orbits, stab, img = class_data(n)
        fixed_all = [k for k in range(len(reps)) if all(img[k][gi] == k for gi in range(len(perms)))]
        ok(not fixed_all, f'B: n={n}: no labelled switching class is fixed by all of S_{n} ({len(reps)} classes, {len(orbits)} types)')
        if n < 5:
            continue
        E, canon, info = euler_data(n, perms)
        even_types = [c for c, u in info.items() if u]
        adj = {}
        for k, _ in orbits:
            G = stab[k]
            types = set()
            for F in E:
                if info[canon[F]] and all(act_graph(F, g) == F for g in G):
                    types.add(canon[F])
            adj[k] = types
        nontriv = sum(1 for k, _ in orbits if len(stab[k]) > 1)
        mm = max_matching([k for k, _ in orbits], adj)
        expect = {5: (2, 2, 1), 6: (6, 6, 5)}[n]
        ok((len(orbits), nontriv, mm) == expect and len(even_types) == len(orbits),
           f'A: n={n}: {len(orbits)} class types ({nontriv} with nontrivial Aut), {len(even_types)} even Euler types, '
           f'maximum stabiliser-containment matching {mm} < {len(orbits)}: no natural bijection')
        if n == 5:
            # C: hand proof data
            orders = [sorted(set(order(g) for g in stab[k])) for k, _ in orbits]
            ok(sorted(max(o) for o in orders) == [3, 5], f'C: n=5 class stabiliser element orders {orders} (an element of order 3, resp. 5)')
            with3 = sorted({canon[F] for F in E if any(order(g) == 3 and act_graph(F, g) == F for g in perms)})
            with5 = sorted({canon[F] for F in E if any(order(g) == 5 and act_graph(F, g) == F for g in perms)})
            ev3 = [c for c in with3 if info[c]]
            ev5 = [c for c in with5 if info[c]]
            ok(len(with3) == 4 and len(with5) == 3 and ev3 == [()] and ev5 == [()],
               f'C: Euler graphs on 5 vertices with an order-3 automorphism: {len(with3)}, order-5: {len(with5)}; the only even one is the empty graph')


def check_D():
    for n in range(2, 7):
        perms = list(itertools.permutations(range(n)))
        bad = 0
        nb = 0
        for F in euler_graphs(n):
            # bipartite?
            col = {}
            bip = True
            adjl = defaultdict(list)
            for (i, j) in F:
                adjl[i].append(j)
                adjl[j].append(i)
            for s in range(n):
                if s in col:
                    continue
                col[s] = 0
                stack = [s]
                while stack:
                    x = stack.pop()
                    for y in adjl[x]:
                        if y not in col:
                            col[y] = 1 - col[x]
                            stack.append(y)
                        elif col[y] == col[x]:
                            bip = False
            if not bip:
                continue
            nb += 1
            if any(act_graph(F, g) == F and sgn(F, g) == -1 for g in perms):
                bad += 1
        ok(bad == 0, f'D: n={n}: all {nb} labelled bipartite Euler graphs are even')


def check_E():
    A000568 = [1, 1, 2, 4, 12]
    for n in range(1, 6):
        perms = list(itertools.permutations(range(n)))
        P = pairs(n)
        seen = set()
        even = 0
        for m in range(1 << len(P)):
            F = frozenset(P[k] for k in range(len(P)) if m >> k & 1)
            c = graph_canon(F, perms)
            if c in seen:
                continue
            seen.add(c)
            aut = [g for g in perms if act_graph(F, g) == F]
            if all(sgn(F, g) == 1 for g in aut):
                even += 1
        ok(even == A000568[n - 1], f'E: n={n}: even graphs {even} = tournaments A000568 = {A000568[n - 1]}')


def check_F():
    expected_t = {4: [1, 1, 3, 8], 3: [1, 1, 2, 4]}
    for l in (3, 4):
        for n in range(1, 5):
            P = pairs(n)
            m = len(P)
            idx = {p: k for k, p in enumerate(P)}
            def enc(vec):
                x = 0
                for v in reversed(vec):
                    x = x * l + v
                return x
            def dec(x):
                vec = []
                for _ in range(m):
                    vec.append(x % l)
                    x //= l
                return vec
            N = l ** m
            # permutation action on upper-triangle entries of a skew matrix: entry (i,j), i<j, value a means
            # M[i][j] = a, M[j][i] = -a
            def act(vec, g):
                new = [0] * m
                for (i, j), k in idx.items():
                    a, b = g[i], g[j]
                    if a < b:
                        new[idx[(a, b)]] = vec[k]
                    else:
                        new[idx[(b, a)]] = (-vec[k]) % l
                return new
            def switch(vec, v):
                new = vec[:]
                for (i, j), k in idx.items():
                    if i == v:
                        new[k] = (new[k] + 1) % l
                    elif j == v:
                        new[k] = (new[k] - 1) % l
                return new
            gens = [tuple(range(n))]
            for i in range(n - 1):
                g = list(range(n)); g[i], g[i + 1] = g[i + 1], g[i]
                gens.append(tuple(g))
            parent = list(range(N))
            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x
            def union(a, b):
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb
            for x in range(N):
                vec = dec(x)
                for g in gens:
                    union(x, enc(act(vec, g)))
                for v in range(n):
                    union(x, enc(switch(vec, v)))
            s = len({find(x) for x in range(N)})
            # t: modular Eulerian matrices (row sums = 0 mod l), orbits under S_n only
            def rowsums(vec):
                r = [0] * n
                for (i, j), k in idx.items():
                    r[i] += vec[k]
                    r[j] -= vec[k]
                return [x % l for x in r]
            K = [x for x in range(N) if not any(rowsums(dec(x)))]
            Kset = set(K)
            parent2 = {x: x for x in K}
            def find2(x):
                while parent2[x] != x:
                    parent2[x] = parent2[parent2[x]]
                    x = parent2[x]
                return x
            for x in K:
                vec = dec(x)
                for g in gens:
                    y = enc(act(vec, g))
                    assert y in Kset
                    ra, rb = find2(x), find2(y)
                    if ra != rb:
                        parent2[ra] = rb
            t = len({find2(x) for x in K})
            ok(s == t == expected_t[l][n - 1], f'F: l={l}, n={n}: s = {s}, t = {t} (expected {expected_t[l][n - 1]})')


def check_G():
    for n in (3, 5, 7):
        P = pairs(n)
        reps = list(tiling_reps(n))
        bad = 0
        for t in reps:
            T = dict(zip(P, t))
            cnt = 0
            for Lm in range(1 << n):
                # switch at L
                out = [0] * n
                for (i, j), fwd in T.items():
                    f = fwd ^ bool(((Lm >> i) ^ (Lm >> j)) & 1)
                    if f:
                        out[i] += 1
                    else:
                        out[j] += 1
                if all((2 * o - (n - 1)) % 4 == 0 for o in out):
                    cnt += 1
            # L and its complement give the same tournament: each member counted twice
            if cnt != 2:
                bad += 1
        ok(bad == 0, f'G: n={n}: every one of the {len(reps)} switching classes has exactly one member with out = in (mod 4) everywhere')


def main():
    t0 = time.time()
    print('==== A/B/C ====', flush=True); check_A_B_C()
    print('==== D ====', flush=True); check_D()
    print('==== E ====', flush=True); check_E()
    print('==== F ====', flush=True); check_F()
    print('==== G ====', flush=True); check_G()
    print(f'elapsed {time.time() - t0:.0f} s')
    print('ALL CHECKS PASSED' if all(OKS) else 'SOME CHECK FAILED')


if __name__ == '__main__':
    main()
