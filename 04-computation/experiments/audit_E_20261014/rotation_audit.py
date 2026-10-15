#!/usr/bin/env python3
"""Audit E (2026-10-14), THM-4655 statement 2 (rotation lemma) and its corollary.
(1) For every cyclic arrangement C of [n] (n <= 9) and sigma = identity, cut C before each vertex, and count DIRECTLY
    (by building the n paths) the pairs i<j with i before j.  Check: the number of cuts with i before j equals
    (c_i - c_j) mod n; the total is <= n(n-1)(2n-1)/6 with equality only for the identity cycle (in sigma-order);
    the average over cuts and pairs is <= (2n-1)/(3n) < 2/3.
(2) Symbolic check of the corollary: 2/3 < f (2n-1)/(3n) + (1-f)  <=>  f < n/(n+1)."""
from itertools import permutations
from fractions import Fraction as F
import sympy as sp

for n in range(2, 10):
    best = -1; argbest = []; formula_ok = True
    for perm in permutations(range(1, n)):
        cyc = (0,) + perm                       # cyclic order; vertex labels = sigma ranks
        c = {v: k for k, v in enumerate(cyc)}
        conc = 0
        per_pair = {}
        for s in range(n):                      # cut before position s
            p = [cyc[(s + t) % n] for t in range(n)]
            ps = {v: t for t, v in enumerate(p)}
            for i in range(n):
                for j in range(i + 1, n):
                    if ps[i] < ps[j]:
                        conc += 1
                        per_pair[(i, j)] = per_pair.get((i, j), 0) + 1
        for i in range(n):
            for j in range(i + 1, n):
                if per_pair.get((i, j), 0) != (c[i] - c[j]) % n:
                    formula_ok = False
        if conc > best:
            best, argbest = conc, [cyc]
        elif conc == best:
            argbest.append(cyc)
    bound = n * (n - 1) * (2 * n - 1) // 6
    avg = F(best, n * n * (n - 1) // 2)
    ident = tuple(range(n))
    print(f"n={n}: cut-count formula (c_i-c_j) mod n {'holds' if formula_ok else 'FAILS'}; max concordance {best} "
          f"(bound n(n-1)(2n-1)/6 = {bound}); average {avg} vs (2n-1)/(3n) = {F(2*n-1, 3*n)}; "
          f"maximisers {argbest} {'= identity only' if argbest == [ident] else ''}")

f, n = sp.symbols('f n', positive=True)
expr = sp.simplify(f * (2 * n - 1) / (3 * n) + (1 - f) - sp.Rational(2, 3))
print("f(2n-1)/(3n) + (1-f) - 2/3 =", sp.factor(expr), "; positive iff f <", sp.solve(sp.Eq(expr, 0), f))
