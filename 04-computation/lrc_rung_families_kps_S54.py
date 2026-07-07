#!/usr/bin/env python3
"""
kps-2026-07-07-S54 (part 2) -- DEEP CENSUS of the rung families of the single-scale
13-family M-spectrum.  The S54 census found the bottom is a FAREY LADDER
1/14 < 2/27 < 1/13 < 2/25 < 1/12 < ... with the AP uniquely at 1/14.  Now characterize the
finite/infinite parametric families at each rung and hunt closed-form M (the n=13 analogue
of the n=12 residue-liar formula M({1..11,13,12k}) = k/(12k+5)).

  (1) SPECTRUM QUANTIZATION: collect ALL distinct M-values <= 0.11 for single-scale
      13-families; confirm they are Farey rungs (three-gap quantization, mac-mini HYP-4412).
  (2) ONE-SWAP families {1..13}\{j} u {X}: M as a function of (j, X); closed forms; which
      are near-tight; finite (bounded X) vs infinite (all X) lonely.
  (3) The 2/27 first-excited family: identify + closed form.
"""
from fractions import Fraction
import numpy as np
from math import gcd
from functools import reduce
from itertools import combinations
from collections import defaultdict

def Mw(v):
    v = [x for x in v if x]
    S = sum(abs(x) for x in v)
    Q = min(4 * S, 2 * max(abs(x) for x in v) + 2)
    va = np.array(v, dtype=np.int64)
    bn, bd = 0, 1
    for q in range(2, Q + 1):
        a = np.arange(1, q)
        r = np.outer(va, a) % q
        d = np.minimum(r, q - r)
        col = d.min(axis=0)
        qb = int(col.max())
        if qb * bd > bn * q:
            bn, bd = qb, q
    return Fraction(bn, bd)

print("=" * 78)
print("(1) SPECTRUM QUANTIZATION: distinct M-values <= 0.11 for single-scale 13-families")
print("=" * 78)
seen = defaultdict(int)
for mx in range(13, 20):
    for combo in combinations(range(1, mx), 12):
        v = combo + (mx,)
        if reduce(gcd, v) != 1:
            continue
        M = Mw(v)
        if M <= Fraction(11, 100):
            seen[M] += 1
print("  distinct M-values in (1/14, 0.11]:")
prev = None
allfarey = True
for M in sorted(seen):
    # is M a mediant-neighbour of the Farey ladder? check small denominator
    tag = ""
    if M.denominator <= 30:
        tag = "(small-denom Farey rung)"
    else:
        tag = "(LARGER denom -- check)"; allfarey = False
    med = ""
    if prev is not None:
        m = Fraction(prev.numerator + M.numerator, prev.denominator + M.denominator)
        med = "  mediant(prev,this)=%s" % m
    print("    %-8s = %.6f  x%-4d %s%s" % (M, float(M), seen[M], tag, med))
    prev = M
print("  => all small-denominator Farey rungs (three-gap quantized)?", allfarey)

print()
print("=" * 78)
print("(2) ONE-SWAP families {1..13}\\{j} u {X}: M(j, X)")
print("=" * 78)
base = list(range(1, 14))
for j in range(1, 14):
    rest = [x for x in base if x != j]
    # scan X from 14 upward, record M; find low-M X and any closed-form period
    lowX = []
    Mvals = {}
    for X in range(14, 400):
        if X in rest:
            continue
        v = tuple(sorted(rest + [X]))
        if reduce(gcd, v) != 1:
            continue
        M = Mw(list(v))
        Mvals[X] = M
        if M <= Fraction(1, 13):
            lowX.append((X, M))
    # summarize
    minM = min(Mvals.values())
    argmin = [X for X, M in Mvals.items() if M == minM]
    lowstr = ", ".join("X=%d:%s" % (X, M) for X, M in lowX[:6])
    print("  remove %2d: min M=%-7s at X in %s ... low(<=1/13): %s" %
          (j, str(minM), str(argmin[:4]), lowstr if lowstr else "(none)"))

print()
print("=" * 78)
print("(3) The 2/27 first-excited family (remove 10, vary X = c*10?) + closed forms")
print("=" * 78)
rest10 = [x for x in base if x != 10]
print("  {1..9,11,12,13} u {X}, X = 20,30,...,200 (multiples of 10) and neighbours:")
for X in [20, 30, 40, 50, 100, 200, 17, 18, 19, 21, 23]:
    v = tuple(sorted(rest10 + [X]))
    if reduce(gcd, v) != 1:
        print("    X=%-4d  (gcd>1, skip)" % X); continue
    print("    X=%-4d  M=%s = %.6f" % (X, Mw(list(v)), float(Mw(list(v)))))

# hunt: for each j, is there an INFINITE arithmetic family {1..13}\{j} u {X: X in c+q*Z} with const M?
print()
print("  INFINITE-FAMILY hunt: remove j, X = X0 + q*t (t=0..6); is M constant along t?")
for j in [10, 13, 1, 7]:
    rest = [x for x in base if x != j]
    # try to find a period q where M is eventually constant
    for q in [13, 14, 27, 25, 12]:
        Ms = []
        for t in range(0, 7):
            X = 14 + q * t if j != 13 else 14 + q * t
            v = tuple(sorted(rest + [X]))
            if reduce(gcd, v) != 1:
                Ms.append(None); continue
            Ms.append(Mw(list(v)))
        if len(set(m for m in Ms if m is not None)) == 1 and None not in Ms:
            print("    remove %2d, X=14+%d*t: M CONSTANT = %s for t=0..6" % (j, q, Ms[0]))
print()
print("DONE -- structural readout in the session reflection.")
