#!/usr/bin/env python3
"""Circulant tournaments versus even circulant graphs of prime order p.
Brute force (p <= PMAX_BF): isomorphism types = multiplier classes (Turner 1967 / Elspas-Turner 1970 for prime order);
evenness of circ(p,S) via Alspach's theorem Aut = Z_p x| Stab(S) for S not in {empty, all}: an automorphism x -> ax+b
with ord(a) even reverses |S|/ord(a) edge orbits, so circ(p,S) is even iff 2 ord(a) | |S| for every even-order a in Stab(S).
The criterion is checked against networkx automorphism groups for p <= 13.
Formula: t_p = (1/(p-1)) sum_{d | p-1, d odd} phi(d) 2^{(p-1)/(2d)};  e_p by Moebius over the subgroups containing -1.
Also: refined counts per odd multiplier subgroup K (types with stabiliser containing K), and tests of uniform rules
S -> D_c(S) = {x : [x in S] != [cx in S]}."""
import sys, itertools
from sympy import isprime, totient, divisors, primitive_root
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
PMAX_BF = int(sys.argv[1]) if len(sys.argv) > 1 else 31
PMAX_F = int(sys.argv[2]) if len(sys.argv) > 2 else 200

def order(a, p):
    k = 1; x = a % p
    while x != 1: x = x * a % p; k += 1
    return k

def t_formula(p):
    return sum(totient(d) * 2 ** ((p - 1) // (2 * d)) for d in divisors(p - 1) if d % 2 == 1) // (p - 1)

def e_formula(p):
    hs = [h for h in divisors(p - 1) if h % 2 == 0]  # orders of subgroups containing -1
    A = {}; g = {}
    for h in sorted(hs, reverse=True):
        ups = [h2 for h2 in hs if h2 > h and h2 % h == 0]
        A[h] = 2 ** ((p - 1) // h) - sum(A[h2] for h2 in ups)
        F = 2 ** ((p - 1) // h - 1)
        g[h] = F - sum(A[h2] if (h2 // h) % 2 == 0 else g[h2] for h2 in ups)
    # types: each S with Stab = H lies in a multiplier orbit of size (p-1)/h
    tot = 0
    for h in hs:
        assert g[h] * h % (p - 1) == 0
        tot += g[h] * h // (p - 1)
    return tot

def brute(p):
    units = list(range(1, p))
    classes = [x for x in range(1, (p + 1) // 2)]  # representatives of {x, -x}
    m = len(classes)
    def canon(S):
        return min(tuple(sorted((a * x) % p for x in S)) for a in units)
    # tournaments
    tours = set()
    for signs in itertools.product([0, 1], repeat=m):
        S = frozenset((x if s == 0 else p - x) for x, s in zip(classes, signs))
        tours.add(canon(S))
    # refined: stabiliser subgroup of each tournament type
    def stab(S): return frozenset(a for a in units if frozenset((a * x) % p for x in S) == S)
    t_by_stab = {}
    for c in tours:
        K = stab(frozenset(c)); t_by_stab[K] = t_by_stab.get(K, 0) + 1
    # even circulant graphs
    def is_even(S):
        if len(S) == p - 1: return False
        if len(S) == 0: return True
        H = stab(S)
        return all((2 * order(a, p)) % 1 == 0 and len(S) % (2 * order(a, p)) == 0 for a in H if order(a, p) % 2 == 0)
    evens = set(); e_by_stab = {}
    for sub in itertools.product([0, 1], repeat=m):
        S = frozenset()
        for x, s in zip(classes, sub):
            if s: S = S | {x, p - x}
        if is_even(S):
            c = canon(S)
            if c not in evens:
                evens.add(c); H = stab(S); e_by_stab[H] = e_by_stab.get(H, 0) + 1
    return tours, t_by_stab, evens, e_by_stab

def nx_even(p, S):
    if len(S) == 0: return True
    if len(S) == p - 1: return False
    G = nx.Graph(); G.add_nodes_from(range(p))
    for x in range(p):
        for d in S: G.add_edge(x, (x + d) % p)
    edges = [(min(u, v), max(u, v)) for u, v in G.edges()]
    for mp in GraphMatcher(G, G).isomorphisms_iter():
        s = 0
        for u, v in edges:
            if mp[u] > mp[v]: s ^= 1
        if s: return False
    return True

print("p | t_p (brute, formula) | e_p (brute, formula) | refined per odd K: types with Stab >= K (tournaments / even graphs)")
for p in [q for q in range(3, PMAX_F + 1) if isprime(q)]:
    tf, ef = t_formula(p), e_formula(p)
    line = f"{p} | {tf} | {ef} | {'EQUAL' if tf == ef else 'DIFFER'}"
    if p <= PMAX_BF:
        tours, tbs, evens, ebs = brute(p)
        line = f"{p} | t ({len(tours)}, {tf}) | e ({len(evens)}, {ef}) | {'EQUAL' if len(tours) == len(evens) == tf == ef else 'DIFFER'}"
        # refined: for each odd subgroup K of Z_p^*
        units = list(range(1, p))
        oddK = set()
        for a in units:
            if order(a, p) % 2 == 1:
                K = set(); x = 1
                while True:
                    K.add(x); x = x * a % p
                    if x == 1: break
                oddK.add(frozenset(K))
        ref = []
        for K in sorted(oddK, key=len):
            tk = sum(v for H, v in tbs.items() if K <= H)
            ek = sum(v for H, v in ebs.items() if K <= H)
            ref.append(f"|K|={len(K)}: {tk}/{ek}{'' if tk == ek else ' <-- DIFFER'}")
        line += " ; " + ", ".join(ref)
        if p <= 13:
            # verify the evenness criterion against networkx
            classes = [x for x in range(1, (p + 1) // 2)]
            bad = 0
            for sub in itertools.product([0, 1], repeat=len(classes)):
                S = frozenset()
                for x, s in zip(classes, sub):
                    if s: S = S | {x, p - x}
                # recompute criterion
                H = [a for a in units if frozenset((a * x) % p for x in S) == S]
                crit = (len(S) == 0) or (len(S) < p - 1 and all(len(S) % (2 * order(a, p)) == 0 for a in H if order(a, p) % 2 == 0))
                if crit != nx_even(p, S): bad += 1
            line += f" ; criterion vs networkx: {'OK' if bad == 0 else str(bad) + ' MISMATCHES'}"
    print(line)

print("\nUniform rules S -> D_c(S): does it induce a bijection tournament types -> even circulant types? (p <= 23)")
for c in [2, -2, 4, -4, 3, -3]:
    res = []
    for p in [q for q in range(5, 24) if isprime(q)]:
        tours, tbs, evens, ebs = brute(p)
        units = list(range(1, p))
        def canon(S): return min(tuple(sorted((a * x) % p for x in S)) for a in units)
        cc = c % p
        images = set(); ok = True
        for T in tours:
            S = frozenset(T)
            D = frozenset(x for x in units if (x in S) != ((cc * x) % p in S))
            img = canon(D)
            if img not in evens or img in images: ok = False
            images.add(img)
        res.append(f"p={p}:{'bij' if ok and len(images) == len(evens) else 'no'}")
    print(f"  c = {c}: " + " ".join(res))
