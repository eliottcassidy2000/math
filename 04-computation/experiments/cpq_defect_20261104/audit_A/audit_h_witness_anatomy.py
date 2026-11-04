#!/usr/bin/env python3
"""audit_A: witness anatomy claims of the note (section 2.4) and HYP-9463 (1), using audit_e2_orbit's enumeration.
  * (3,7) k = 4: the note says the 17 witnesses are one K_{7,7} and sixteen K_{3,3}-type configurations
    (2 cross x 2 circulant classes x 4 colours of the extra point).
  * first-defect rows: HYP-9463 (1) says the witness at k = d is unique with E_even = E_wit = E_T = 0; check which
    rows actually satisfy that ((3,5) k=2, (5,7) k=2, (7,11) k=4, (11,13) k=2 versus (3,7) k=4, (9,q) rows).
"""
import sys, itertools
from collections import defaultdict
import pynauty
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from audit_e2_orbit import dil_classes, mult_even, coloured_types, full_object

def witnesses(p, q, k):
    n = p + q + k
    eP = [C for C in dil_classes(p, False) if mult_even(p, C)]; eQ = [C for C in dil_classes(q, False) if mult_even(q, C)]
    cGe = [x for x in coloured_types(k, False) if x[2]]
    types = {}
    for CP in eP:
        for CQ in eQ:
            for cross in (0, 1):
                for (obj, col, _) in cGe:
                    adj = full_object(p, q, k, CP, CQ, cross, obj, col, False)
                    g = pynauty.Graph(n, directed=False, adjacency_dict=adj)
                    cert = pynauty.certificate(g)
                    if cert in types:
                        types[cert][0] += 1; types[cert][2].append((CP, CQ, cross, obj, col)); continue
                    edges = [(u, v) for u in adj for v in adj[u] if u < v]
                    gens = pynauty.autgrp(g)[0]
                    even = all((sum(1 for (u, v) in edges if gg[u] > gg[v]) & 1) == 0 for gg in gens)
                    types[cert] = [1, even, [(CP, CQ, cross, obj, col)]]
    return types

def describe(p, q, k, d):
    CP, CQ, cross, obj, col = d
    # K_{3,3}-type: C_P empty, exactly three fixed points of colour (x=1, y=cross) pairwise non-adjacent
    xs = [(c & 1, c >> 1) for c in col]
    triple = [t for t in range(k) if xs[t] == (1, cross)]
    k33 = (len(CP) == 0 and len(triple) >= 3)
    kqq = (len(CP) == 0 and cross == 1 and sum(1 for t in range(k) if xs[t][1] == 1) == q - p and all(xs[t] in ((0, 1), (0, 0)) for t in range(k)))
    return "C_P=%s C_Q=%s cross=%d F-edges=%s colours(x,y)=%s %s" % (sorted(CP), sorted(CQ), cross, sorted(obj), xs,
            "[K_{q,q}-type]" if kqq else ("[K_{3,3}-type]" if k33 else "[other]"))

if __name__ == "__main__":
    for (p, q, k) in [(3, 5, 2), (5, 7, 2), (3, 7, 4), (7, 11, 4), (11, 13, 2), (9, 7, 2), (9, 5, 4)]:
        types = witnesses(p, q, k)
        W = [(cert, v) for cert, v in types.items() if not v[1]]
        E_even = sum(v[0] - 1 for v in types.values() if v[1]); E_wit = sum(v[0] - 1 for _, v in W)
        print("(%d,%d) k=%d: W=%d E_even=%d E_wit=%d" % (p, q, k, len(W), E_even, E_wit))
        if (p, q, k) == (3, 7, 4) or len(W) <= 2:
            for cert, v in W:
                print("   ", describe(p, q, k, v[2][0]), "(N-even orbits in type: %d)" % v[0])
        if (p, q, k) == (3, 7, 4):
            kinds = defaultdict(int)
            for cert, v in W: kinds[describe(p, q, k, v[2][0]).split('[')[1]] += 1
            print("    kinds:", dict(kinds))
