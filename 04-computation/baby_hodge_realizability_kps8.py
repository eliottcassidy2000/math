#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Baby Hodge: the realizable invariant-vector region of tournaments, as a moment problem.
kind-pasteur-2026-06-14-S8.  Builds on THM-118 (c_k=tr A^k/k), THM-499 (H=1+2(c3+c5)+4D),
THM-502 (closed-walk census), THM-462 (c3 gap-free).

The "cycle-class map":  T  |->  (c3, c5, [c7], H)  =  the tournament's invariant vector.
"Baby Hodge": which integer vectors are in the IMAGE (realized by an actual T)?  The
HOLES (realizable-region lattice points not hit) are the combinatorial analog of a Hodge
class that is NOT algebraic.  Power-sum reframe: (3c3, 5c5) = (tr A^3, tr A^5) = MOMENTS
of the tournament spectrum; the realizable region = a truncated moment problem; the
moment/positivity boundary = "Hodge inequalities"; integrality + the conflict-graph layer
(alpha_2, the THM-499 D) = the discrete refinement that cuts the holes.

COMPUTE (exact, trace speedup; exhaustive n<=7, sampled n=8):
 (A) the (c3,c5) realizable region per n: bounding box, per-c3 c5-fibers, and ALL HOLES
     (lattice points inside [0,maxc5] of a nonempty c3-fiber that are NOT realized).
 (B) classify each hole as INTERIOR (between realized c5 values in its c3-fiber = a
     'non-algebraic Hodge class', the continuous moment region permits it) vs a boundary
     artifact.
 (C) the spectral/conflict-graph test: holes should persist in (c3,c5) [spectral coords]
     but the FULL vector (c3,c5,H) realizable set is finer (H adds the alpha_2 layer).
"""

import sys, itertools
import numpy as np
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None


def all_tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        A = np.zeros((n, n), dtype=np.int64)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
        yield A


def cc(A):
    A2 = A @ A; A3 = A2 @ A; A4 = A2 @ A2; A5 = A4 @ A
    return int(np.trace(A3)) // 3, int(np.trace(A5)) // 5


def region(n):
    from collections import defaultdict
    fibers = defaultdict(set)   # c3 -> set of c5
    for A in all_tournaments(n):
        c3, c5 = cc(A)
        fibers[c3].add(c5)
    return fibers


def holes_of(fibers):
    """interior holes: per c3-fiber, lattice points strictly between min and max realized c5
    that are NOT realized = the 'non-algebraic Hodge classes'."""
    H = []
    for c3 in sorted(fibers):
        cs = sorted(fibers[c3])
        for v in range(cs[0], cs[-1] + 1):
            if v not in fibers[c3]:
                H.append((c3, v))
    return H


def main():
    print("=== Baby Hodge: the realizable (c3,c5) region and its HOLES (non-algebraic classes) ===\n", flush=True)
    for n in (5, 6, 7):
        fib = region(n)
        allc3 = sorted(fib)
        allc5 = sorted(set().union(*fib.values()))
        holes = holes_of(fib)
        print(f"n={n}: c3 in [{allc3[0]},{allc3[-1]}], c5 in [{allc5[0]},{allc5[-1]}]", flush=True)
        print(f"   per-c3 c5-fibers:", flush=True)
        for c3 in allc3:
            cs = sorted(fib[c3])
            gaps = [v for v in range(cs[0], cs[-1]+1) if v not in fib[c3]]
            print(f"      c3={c3:2d}: c5 in {cs}" + (f"   HOLES {gaps}" if gaps else ""), flush=True)
        print(f"   INTERIOR HOLES (c3,c5) [non-algebraic Hodge classes]: {holes}", flush=True)
        # global c5 holes (a c5 value realized by NO tournament but inside [0,max])
        gh = [v for v in range(allc5[-1]+1) if v not in allc5]
        print(f"   GLOBAL c5 holes (no tournament at all): {gh}", flush=True)
        print(flush=True)


if __name__ == "__main__":
    main()
