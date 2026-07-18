#!/usr/bin/env python3
"""
Content law, uniform-in-n attempt  (mac-mini-2026-07-18-S111) — THM-1006 sec.I
==============================================================================
(I1) CLEAN FORM (proved): for a tight n-set, val=1 <=> no element divisible by n+1.
     => CONTENT LAW: a primitive tight n-set contains no multiple of n+1.
(I2) UNIFORM FAMILY EXCLUDED (proved): M({1..n-1} u {k n(n+1)}) = k(n+1)/(k n(n+1)+1)
     > 1/(n+1) for all n>=3, k>=1 -- never tight.
(I4) Verdict: full uniform law NOT proved; residue is a metric/equidistribution
     input, provably beyond counting (THM-1006 sec.H).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

def M_full(S):
    S = sorted(set(S)); best = F(0); bq = None; bval = None
    dens = {2 * a for a in S}
    for a, b in combinations(S, 2):
        dens.add(a + b); dens.add(b - a)
    dens.discard(0)
    for q in dens:
        for m in range(1, q):
            num = q
            for v in S:
                r = (v * m) % q; d = min(r, q - r)
                if d < num: num = d
                if num == 0: break
            if num > 0:
                c = F(num, q)
                if c > best: best = c; bq = q; bval = num
    return best, bval, bq

print("(I1) clean form: tight n-sets containing (n+1) are all IMPRIMITIVE?")
for n in range(3, 8):
    N = 3 * n + 2; L = F(1, n + 1); res = []
    for A in combinations(range(1, N + 1), n):
        if (n + 1) not in A: continue
        if M_full(A)[0] == L: res.append((list(A), reduce(gcd, A)))
    prim = [a for a, g in res if g == 1]
    print(f"  n={n}: tight sets containing {n+1}: {len(res)}, primitive: {len(prim)}  {res[:3]}")

print()
print("(I2) uniform family: M({1..n-1} u {k n(n+1)}) == k(n+1)/(k n(n+1)+1)?")
ok = True
for n in range(3, 10):
    for k in (1, 2, 3, 4):
        A = list(range(1, n)) + [k * n * (n + 1)]
        if len(set(A)) != n: continue
        M = M_full(A)[0]; pred = F(k * (n + 1), k * n * (n + 1) + 1)
        if M != pred: ok = False; print(f"   MISMATCH n={n} k={k}: {M} vs {pred}")
print(f"  closed form verified on all tested (n,k): {ok}")
print("  > 1/(n+1) always since k(n+1)^2 > k n(n+1)+1 <=> k(n+1) > 1.  => never tight.")

print()
print("(I4) VERDICT: content law NOT proved uniformly. Counting is provably insufficient")
print("     (THM-1006 sec.H: capacity+primitivity satisfiable for every val=2..13);")
print("     the residue is a metric/equidistribution input (klein CRUX(C), codex-S64).")
