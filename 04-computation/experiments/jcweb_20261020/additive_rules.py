#!/usr/bin/env python3
"""HYP-9342: (1) uniqueness of Gamma-equivariant type bijections at p = 7, 11, 19 (count perfect matchings of the Gamma
stabiliser-containment graph on Gamma-orbits; the Paley rule beta is one of them); (2) additive candidate rules at
p = 13, 17, 29: R1(S) = {x : |S cap (x + S)| odd}, R2(S) = {x : |S cap (x - S)| odd}, R3(S) = {x : #{(a,b) in S^2 : a + b = x} odd},
R4(S) = {x : |S cap (2x - S)| odd} (3-term progressions centred at x); test whether the rule is symmetric, even-valued, and
bijective on multiplier types."""
import itertools
from math import gcd
from sympy import isprime
def order(a, p):
    k = 1; x = a % p
    while x != 1: x = x * a % p; k += 1
    return k
def setup(p):
    units = list(range(1, p)); classes = list(range(1, (p + 1) // 2))
    tours = set()
    for signs in itertools.product([0, 1], repeat=len(classes)):
        tours.add(frozenset((x if sg == 0 else p - x) for x, sg in zip(classes, signs)))
    def is_even(S):
        if len(S) == 0: return True
        if len(S) == p - 1: return False
        H = [a for a in units if frozenset(a * x % p for x in S) == S]
        return all(len(S) % (2 * order(a, p)) == 0 for a in H if order(a, p) % 2 == 0)
    evens = set()
    for sub in itertools.product([0, 1], repeat=len(classes)):
        S = frozenset()
        for x, sg in zip(classes, sub):
            if sg: S = S | {x, p - x}
        if is_even(S): evens.add(S)
    def canon(S): return min(tuple(sorted(a * x % p for x in S)) for a in units)
    return units, tours, evens, canon
# (2) additive rules
for p in (13, 17, 29):
    units, tours, evens, canon = setup(p)
    ttypes = {canon(S) for S in tours}; etypes = {canon(S) for S in evens}
    def R1(S): return frozenset(x for x in units if sum(1 for a in S if (x + a) % p in S) % 2)
    def R2(S): return frozenset(x for x in units if sum(1 for a in S if (a - x) % p in S) % 2)
    def R3(S): return frozenset(x for x in units if sum(1 for a in S if (x - a) % p in S) % 2)
    def R4(S): return frozenset(x for x in units if sum(1 for a in S if (2 * x - a) % p in S) % 2)
    for name, R in (("R1 |S∩(x+S)| odd", R1), ("R2 |S∩(S+x)| odd", R2), ("R3 reps of x as a+b odd", R3), ("R4 3-AP centred x odd", R4)):
        imgs = {}; symmetric = True; evenval = True
        for S in tours:
            X = R(S)
            if frozenset((p - x) % p for x in X) != X: symmetric = False
            if X not in evens: evenval = False
            imgs.setdefault(canon(S), set()).add(canon(X) if symmetric else None)
        welldef = all(len(vv) == 1 for vv in imgs.values())
        injective = welldef and len({next(iter(vv)) for vv in imgs.values()}) == len(imgs)
        print(f"p = {p} {name}: symmetric {symmetric}, even-valued {evenval}, well defined on types {welldef}, injective on types {injective} (types {len(ttypes)} vs even types {len(etypes)})")
# (1) uniqueness of Gamma-equivariant type bijections at p = 7, 11, 19 (perfect matchings of the Lemma N1 graph)
import networkx as nx
for p in (7, 11, 19):
    units, tours, evens, canon = setup(p)
    js = [j for j in range(1, p - 1) if gcd(j, p - 1) == 1]
    Gam = [(a, j) for a in units for j in js]
    def act(S, a, j): return frozenset(a * pow(x, j, p) % p for x in S)
    def orbits(sets):
        seen = set(); orbs = []
        for S in sets:
            if S in seen: continue
            O = {act(S, a, j) for a, j in Gam}; seen |= O; orbs.append(O)
        return orbs
    to, eo = orbits(tours), orbits(evens)
    def stab(S): return frozenset(g for g in Gam if act(S, *g) == S)
    B = nx.Graph()
    for i, O in enumerate(to):
        T = next(iter(O)); ST = stab(T)
        for k, OX in enumerate(eo):
            if any(ST <= stab(Y) for Y in OX): B.add_edge(('T', i), ('X', k))
    # count perfect matchings by brute force (small)
    Ts = [('T', i) for i in range(len(to))]
    def count(i, used):
        if i == len(Ts): return 1
        tot = 0
        for x in B[Ts[i]]:
            if x not in used: used.add(x); tot += count(i + 1, used); used.remove(x)
        return tot
    npm = count(0, set())
    # the Paley rule image orbit for each tournament orbit
    QR = frozenset(pow(x, 2, p) for x in units)
    def beta(S):
        X = frozenset(x for x in units if (x in S) != (x in QR)); X = X | frozenset(p - x for x in X)
        return X if (len(X) // 2) % 2 == 0 else frozenset(units) - X
    ok = all(any(beta(next(iter(O))) in OX for OX in eo) for O in to)
    print(f"p = {p}: Gamma-orbits T {len(to)} E {len(eo)}; perfect matchings of the Gamma stabiliser-containment graph: {npm}; Paley rule maps orbits into even orbits: {ok}")
