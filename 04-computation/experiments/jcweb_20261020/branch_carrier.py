#!/usr/bin/env python3
"""HYP-9172 (b): decorated carriers for N_lev(n) = sum over switching-class types C of (1 - rho(Aut C)).
For n in (4, 6, 8): list the classes with even |Aut C|; for each: |Aut|, Sylow-2 order, number of Aut-conjugacy classes of
involutions and of nontrivial 2-elements, the per-class contribution 1 - rho, the formula sum_{[x] 2-element != 1} rho(C(x)),
and the decorated flip tau: (C, g) -> (C + M_g, g) on S_n-orbits of pairs (class, involution): fixed or not."""
import sys, itertools
sys.path.insert(0, '/Users/e/Documents/GitHub/math-wt-mersenne-561/04-computation/experiments')
import procgen_tbij_20261001_nat as N
import procgen_tbij_20261001_lib as L
from collections import Counter
def perm_order(p):
    idt = tuple(range(len(p))); x = p; k = 1
    while x != idt: x = N.perm_mul(x, p); k += 1
    return x is idt or k
def mul(a, b): return N.perm_mul(a, b)
def inv(a): return N.inverse(a)
def order(p):
    idt = tuple(range(len(p))); x = p; k = 1
    while x != idt: x = tuple(x[p[i]] for i in range(len(p))); k += 1
    return k
for n in (4, 6, 8):
    orbs = N.class_orbits(n)
    even = [o for o in orbs if o.order % 2 == 0]
    print(f"\n=== n = {n}: classes {len(orbs)}, even-Aut classes {len(even)} ===")
    tot = 0; Pc = []
    for o in even:
        G = N.group_elements(o.gens, n) if o.gens else [tuple(range(n))]
        G = list(G)
        ords = Counter(order(g) for g in G)
        odd = sum(c for k, c in ords.items() if k % 2 == 1)
        rho = odd / len(G); tot += 1 - rho
        twoel = [g for g in G if order(g) > 1 and (order(g) & (order(g) - 1)) == 0]
        # conjugacy classes of 2-elements and involutions
        def cls(x): return frozenset(mul(mul(h, x), inv(h)) for h in G)
        classes2 = {}
        for x in twoel: classes2.setdefault(cls(x), x)
        inv_classes = [x for c, x in classes2.items() if order(x) == 2]
        # formula sum_{[x]} rho(C_G(x))
        s = 0.0
        for c, x in classes2.items():
            Cx = [h for h in G if mul(h, x) == mul(x, h)]
            s += sum(1 for h in Cx if order(h) % 2 == 1) / len(Cx)
        syl = 1
        while len(G) % (2 * syl) == 0: syl *= 2
        # decorated flip on (C, g): is C + M_g isomorphic to C via an h commuting with g?
        T = o.rep
        fixed_flags = []
        for g in inv_classes:
            # M_g as a flip vector: reverse arcs between x and g(x)
            Tm = list(T)
            for v in range(n):
                w = g[v]
                if w > v:
                    # reverse arc between v and w
                    if (Tm[v] >> w) & 1: Tm[v] &= ~(1 << w); Tm[w] |= 1 << v
                    else: Tm[w] &= ~(1 << v); Tm[v] |= 1 << w
            Tm = tuple(Tm)
            c1 = N.class_canon_batch([T])[0]; c2 = N.class_canon_batch([Tm])[0]
            fixed_flags.append(c1 == c2)
        Pc.append((len(G), syl, len(inv_classes), len(classes2), round(1 - rho, 4), round(s, 4), fixed_flags))
    for row in sorted(Pc): print("  |Aut| %4d Syl2 %2d  inv-classes %d  2-elt-classes %d  1-rho %.4f  sum rho(C(x)) %.4f  tau-fixed(C,g) %s" % row)
    print(f"  N_lev({n}) = {tot:.4f}; decorated pairs (C, involution class): {sum(r[2] for r in Pc)}, of which tau-fixed {sum(sum(r[6]) for r in Pc)}")
