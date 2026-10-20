#!/usr/bin/env python3
"""Power-equivariance obstruction for uniform circulant rules. Gamma_p = {x -> a x^j : a in Z_p^*, gcd(j, p-1) = 1}
(the affine group of the multiplicative group Z_p^* = Z_{p-1}). Any rule definable without parameters in the structure
(Z_p^*, ., S) commutes with Gamma_p. Count Gamma_p-orbits on circulant tournament sets and on even circulant graphs:
if they differ, no Gamma_p-equivariant type bijection exists at p. Also the stabiliser-containment Hall test (Lemma N1)."""
import itertools, sys
from math import gcd
from sympy import isprime
import networkx as nx
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 31
def order(a, p):
    k = 1; x = a % p
    while x != 1: x = x * a % p; k += 1
    return k
for p in [q for q in range(5, PMAX + 1) if isprime(q)]:
    units = list(range(1, p)); classes = list(range(1, (p + 1) // 2))
    js = [j for j in range(1, p - 1) if gcd(j, p - 1) == 1]
    Gamma = [(a, j) for a in units for j in js]
    def act(S, a, j): return frozenset(a * pow(x, j, p) % p for x in S)
    def is_even(S):
        if len(S) == 0: return True
        if len(S) == p - 1: return False
        H = [a for a in units if frozenset(a * x % p for x in S) == S]
        return all(len(S) % (2 * order(a, p)) == 0 for a in H if order(a, p) % 2 == 0)
    tours = set(); evens = set()
    for signs in itertools.product([0, 1], repeat=len(classes)):
        tours.add(frozenset((x if s == 0 else p - x) for x, s in zip(classes, signs)))
    for sub in itertools.product([0, 1], repeat=len(classes)):
        S = frozenset()
        for x, s in zip(classes, sub):
            if s: S = S | {x, p - x}
        if is_even(S): evens.add(S)
    def orbits(sets):
        seen = set(); orbs = []
        for S in sets:
            if S in seen: continue
            O = {act(S, a, j) for a, j in Gamma}; seen |= O; orbs.append(O)
        return orbs
    to = orbits(tours); eo = orbits(evens)
    # multiplier types for reference
    def mtypes(sets):
        seen = set(); c = 0
        for S in sets:
            if S in seen: continue
            seen |= {act(S, a, 1) for a in units}; c += 1
        return c
    # Lemma N1 Hall test for Gamma: T ~ X iff Stab_Gamma(T) <= conjugate of Stab_Gamma(X)
    def stab(S): return frozenset((a, j) for a, j in Gamma if act(S, a, j) == S)
    def conj(g, h):  # h g h^-1 in Gamma, as maps x -> a x^j
        (a, j), (b, k) = g, h
        kinv = pow(k, -1, p - 1)
        # h g h^-1 (x) = b (a (b^-1 x)^{k^-1})^{j} ... compute by evaluating on a generator
        from sympy import primitive_root
        return None
    reps_T = [next(iter(O)) for O in to]; reps_E = [next(iter(O)) for O in eo]
    B = nx.Graph()
    for i, T in enumerate(reps_T):
        B.add_node(('T', i)); ST = stab(T)
        for k, X in enumerate(reps_E):
            B.add_node(('X', k))
            # containment up to Gamma-conjugacy: test all conjugates of X's stabiliser via all elements of the orbit of X
            ok = False
            for Y in eo[k]:
                if ST <= stab(Y): ok = True; break
            if ok: B.add_edge(('T', i), ('X', k))
    M = nx.bipartite.maximum_matching(B, top_nodes=[('T', i) for i in range(len(reps_T))])
    matched = sum(1 for i in range(len(reps_T)) if ('T', i) in M)
    print(f"p = {p}: multiplier types T {mtypes(tours)} E {mtypes(evens)} | Gamma-orbits T {len(to)} E {len(eo)} "
          f"{'EQUAL' if len(to) == len(eo) else 'DIFFER'} | Hall (Gamma): matched {matched} of {len(to)}"
          f"{'' if matched == len(to) and len(to) == len(eo) else '  <-- no Gamma-equivariant type bijection'}")
