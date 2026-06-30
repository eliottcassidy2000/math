#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TOURNAMENTS AS INTRANSITIVITY AMONG n THINGS (not graphs among n nodes). mac-mini-2026-06-30-S34.

Reframe: a tournament on n THINGS = a complete preference/dominance relation. Its ONLY content is its
INTRANSITIVITY -- how far it is from a total order (a ranking), and the PATTERN of that failure. Nodes/arcs
are scaffolding; the meaning is "which sets of things have no consistent ranking" = the (odd) cycles.

  - TRANSITIVE tournament  = a total order = ZERO intransitivity = the orderable/"rational" point (the cusp,
    H=1). It is a COBOUNDARY: comes from a ranking potential f (a->b iff f(a)>f(b)).
  - 3-CYCLE = the irreducible quantum of intransitivity = the Condorcet paradox (rock-paper-scissors among
    3 things). Moon: any intransitivity contains a 3-cycle. It is ODD (length 3).
  - intransitivity = the FAILURE of a ranking potential = the cohomology class H^1 = the cycle space.

This script grounds: (1) the order/intransitivity dimension split; (2) the intransitivity SPECTRUM
(cyclicity per iso class); (3) orderability -> 0 (intransitivity is generic, the Condorcet phenomenon);
(4) the count of distinct intransitivity PATTERNS (GF(2) cycle-space shadow = even graphs).
"""
from __future__ import annotations
import functools, itertools
from collections import defaultdict
from math import comb, factorial
print = functools.partial(print, flush=True)


def cyclicity(adj, n):
    """# of cyclic (intransitive) triples = C(n,3) - sum_v C(outdeg(v),2)."""
    out = [sum(1 for j in range(n) if adj[i][j]) for i in range(n)]
    return comb(n, 3) - sum(comb(o, 2) for o in out)


def iso_classes(n):
    prs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(prs); idx = {p: b for b, p in enumerate(prs)}
    P = list(itertools.permutations(range(n)))
    seen = {}
    for bits in range(1 << m):
        adj = [[False]*n for _ in range(n)]
        for b, (i, j) in enumerate(prs):
            if (bits >> b) & 1: adj[i][j] = True
            else: adj[j][i] = True
        canon = min(tuple(1 if adj[s[i]][s[j]] else 0 for i in range(n) for j in range(n) if i != j)
                    for s in P)
        if canon in seen: continue
        seen[canon] = cyclicity(adj, n)
    return seen


def main():
    print("=" * 82)
    print("TOURNAMENTS AS INTRANSITIVITY AMONG n THINGS (mac-mini-S34)")
    print("=" * 82)

    # ---- (1) the dimension split: order vs intransitivity ----
    print("\n[1] A preference among n things splits: ORDER (rankable part) (+) INTRANSITIVITY (the rest):")
    print(f"    {'n':>2} {'pairs C(n,2)':>12} {'ORDER dim (n-1)':>16} {'INTRANSIT dim C(n-1,2)':>23}")
    for n in range(3, 8):
        print(f"    {n:>2} {comb(n,2):>12} {n-1:>16} {comb(n-1,2):>23}")
    print("    ORDER = cut space (a ranking potential, a COBOUNDARY); INTRANSITIVITY = cycle space (H^1).")
    print("    => the n-1 'score' d.o.f. are the orderable part; the C(n-1,2) cycle d.o.f. are pure")
    print("       intransitivity -- the part NO ranking can explain. Intransitivity DOMINATES for large n.")

    # ---- (2) the intransitivity spectrum (cyclicity per iso class) ----
    print("\n[2] INTRANSITIVITY SPECTRUM -- cyclicity (# rock-paper-scissors triples) per iso class:")
    for n in (3, 4, 5, 6):
        cl = iso_classes(n)
        dist = defaultdict(int)
        for c in cl.values(): dist[c] += 1
        spec = dict(sorted(dist.items()))
        zero = spec.get(0, 0)
        cmax = max(spec)
        print(f"    n={n}: cyclicity -> #classes {spec}; TRANSITIVE (c=0): {zero} class (unique, the cusp); "
              f"max intransitivity c={cmax}")
    print("    => c=0 (transitive = the ONE total order up to relabeling) is unique; everything else is")
    print("       intransitive. The metagraph is organized by intransitivity: 0 at the center (cusp).")

    # ---- (3) orderability -> 0: the Condorcet phenomenon ----
    print("\n[3] ORDERABILITY vanishes (intransitivity is GENERIC) -- the Condorcet phenomenon:")
    print(f"    {'n':>2} {'P(transitive=fully orderable)':>30} {'P(has a Condorcet winner)':>27}")
    for n in range(3, 9):
        labeled = 2**comb(n, 2)
        p_trans = factorial(n) / labeled           # transitive labelings / all
        p_cw = n * 2**comb(n-1, 2) / labeled        # a thing beating all others
        print(f"    {n:>2} {p_trans:>30.6f} {p_cw:>27.6f}")
    print("    => both -> 0 fast: a random preference among many things is almost never rankable and")
    print("       usually has no overall winner. INTRANSITIVITY is the rule, not the exception.")

    # ---- (4) how many DISTINCT intransitivity patterns? (GF(2) cycle-space shadow = even graphs) ----
    print("\n[4] # of DISTINCT intransitivity patterns (GF(2) cycle-space shadow, up to relabeling):")
    print("    The undirected support of the intransitivity = an EVEN graph (cycle space). Iso classes:")
    print("    A002854 = 2, 3, 7, 16, 54 (n=3..7) = the # of essentially different intransitivity SHAPES.")
    print("    (Tournaments A000568 = 2,4,12,56,456 > these: a shape hosts several preferences once you")
    print("    overlay an ordering -- order x intransitivity, the labeled 2^(n-1) fibration, S32/HYP-3595.)")

    print("\n" + "=" * 82)
    print("THE ATOM: the 3-cycle (odd, length 3) = the irreducible intransitivity = the Condorcet paradox.")
    print("Tournaments ARE the space of intransitivity patterns among n things; the transitive order is the")
    print("trivial (cohomologically exact) point; the odd cycle generates H^1. The LRC danger relation among")
    print("RUNNERS is the same kind of object -- intransitivity among n things -- and its apex odd cycle C_7")
    print("is the irreducible resonance the floor 4cos^2(3pi/7)>0 (THM-590) certifies as non-degenerate.")
    print("=" * 82)


if __name__ == "__main__":
    main()
