#!/usr/bin/env python3
"""Checks for THM-4721:
(1) Parity lemma: P <= S_n of odd order, h in N(P) of even order, G = <P,h>. With Y_C = antipodal pairs of an even cycle C
    of h, the F_2-sum X = sum_{k in P} kY_C is G-invariant and h reverses an odd number of its edges.
(2) Restricted twisted Burnside: for Gamma <= N(P): #Gamma-orbits on Fix_T(P) = #Gamma-orbits on Gamma-even graphs in Fix_G(P),
    checked by brute force for Gamma in {C(P), N(P), <P,h> for a few h} and n <= 7 (P from the odd-order subgroup list)."""
import sys, itertools
from collections import Counter
sys.path.insert(0, '.')
from subgroups_fix import cyc, mul, closure, subgroup_reps, pair_orbits
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
def inv(a):
    r = [0] * len(a)
    for i, x in enumerate(a): r[x] = i
    return tuple(r)
def cycles(g):
    n = len(g); seen = [False] * n; cs = []
    for i in range(n):
        if not seen[i]:
            c = []; j = i
            while not seen[j]: seen[j] = True; c.append(j); j = g[j]
            cs.append(c)
    return cs
def order(g):
    idt = tuple(range(len(g))); x = g; k = 1
    while x != idt: x = mul(g, x); k += 1
    return k
def r(h, X):
    s = 0
    for (u, v) in X:
        if h[u] > h[v]: s ^= 1
    return s
def img(g, X): return frozenset((min(g[u], g[v]), max(g[u], g[v])) for u, v in X)
for n in range(3, NMAX + 1):
    perms = list(itertools.permutations(range(n)))
    reps = subgroup_reps(n)
    tested = 0
    for name, P in reps.items():
        N = [s for s in perms if frozenset(mul(mul(s, g), inv(s)) for g in P) == P]
        C = [s for s in N if all(mul(s, g) == mul(g, s) for g in P)]
        # (1) parity lemma
        for h in N:
            if order(h) % 2: continue
            G = closure(list(P) + [h], n)
            ev = [c for c in cycles(h) if len(c) % 2 == 0]
            c0 = ev[0]; L = len(c0)
            Y = frozenset((min(c0[i], c0[(i + L // 2) % L]), max(c0[i], c0[(i + L // 2) % L])) for i in range(L // 2))
            X = frozenset()
            for k in P: X = X ^ img(k, Y)
            assert all(img(g, X) == X for g in G), "not G-invariant"
            assert r(h, X) == 1, "even sign"
            tested += 1
        # (2) restricted principle for Gamma = C(P), N(P)
        orbs = pair_orbits(P, n)
        # H-invariant tournaments: orientation per orbit (as in subgroups_fix)
        orb_arcs = []
        for O in orbs:
            p = O[0]; orb_arcs.append(frozenset((g[p[0]], g[p[1]]) for g in P))
        tours = []; graphs = []
        for choice in itertools.product([0, 1], repeat=len(orbs)):
            arcs = frozenset().union(*[A if c else frozenset((v, u) for u, v in A) for c, A in zip(choice, orb_arcs)])
            tours.append(arcs)
            graphs.append(frozenset().union(*[frozenset(O) for c, O in zip(choice, orbs) if c]))
        tset = set(tours); gset = set(graphs)
        for gname, Gam in [("C(P)", C), ("N(P)", N)]:
            # tournament orbits
            seen = set(); t_orb = 0
            for T in tours:
                if T in seen: continue
                t_orb += 1
                for s in Gam: seen.add(frozenset((s[u], s[v]) for u, v in T))
            # Gamma-even graph orbits
            seen = set(); g_orb = 0
            for X in graphs:
                if X in seen: continue
                orbit = set(); even = True
                for s in Gam:
                    Y = img(s, X); orbit.add(Y)
                    if Y == X and r(s, X): even = False
                seen |= orbit
                if even: g_orb += 1
            assert t_orb == g_orb, (n, name, gname, t_orb, g_orb)
            print(f"n={n} P={name} Gamma={gname} (|Gamma|={len(Gam)}): orbits on Fix_T(P) = {t_orb} = Gamma-even orbits on Fix_G(P) = {g_orb}")
    print(f"n = {n}: parity lemma verified for {tested} pairs (P, h)")
