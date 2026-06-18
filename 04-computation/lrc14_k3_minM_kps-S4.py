#!/usr/bin/env python3
"""kind-pasteur-S4: min M over k>=3 covering primitive S3 13-sets (the genuinely OPEN case)
and the structure of the minimizers (binders, mod-7/14, near k/7?). Sampled, feasible. Exact."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random, sys
random.seed(13)
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def MvalAt(S):
    b = F(0); at = None
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v; at = t
    return b, at
def is_cov(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def klarge(S): return sum(1 for v in S if v > 13)

CAP = 130
best = (F(1), None, None); cnt = 0; allmins = []
smalls = list(range(1, 14))
for drop in [3, 4, 5, 6]:
    Dlist = list(combinations(smalls, drop)); random.shuffle(Dlist)
    for D in Dlist[:18]:
        P = [s for s in smalls if s not in D]
        pool = list(range(14, CAP+1))
        for _ in range(450):
            L = random.sample(pool, drop)
            S = sorted(P + L)
            if len(S) != 13 or reduce(gcd, S) != 1 or not is_cov(S): continue
            if klarge(S) < 3 or max(S) < 13*min(S): continue
            cnt += 1
            M, at = MvalAt(S)
            if M < best[0]: best = (M, S, at)
            allmins.append((M, S, at))
print("k>=3 covering primitive S3 sets sampled (Vmax<=%d): %d" % (CAP, cnt)); sys.stdout.flush()
if best[1]:
    M, S, at = best
    b = [v for v in S if nrm(v*at) == M]
    print("MIN M(k>=3) = %s = %.5f   M*14 = %.4f   (1/14 = %.5f)" % (M, float(M), float(M*14), float(F(1,14))))
    print("  minimizer S =", S)
    print("  witness tau =", at, "  binders =", b)
    print("  binders: mod7 =", [v%7 for v in b], " mod14 =", [v%14 for v in b], " large? =", [v>13 for v in b])
    n7 = min((abs(at-F(k,7)), k) for k in range(7))
    print("  dist(tau, k/7) = %.5f (k=%d)" % (float(n7[0]), n7[1]))
    allmins.sort(key=lambda x: x[0])
    print("  12 tightest k>=3 sets (M, binders, any-large-binder):")
    for M, S, at in allmins[:12]:
        b = [v for v in S if nrm(v*at) == M]
        print("    M=%.5f=%s  binders=%s  large-binder=%s  mod7(b)=%s"
              % (float(M), M, b, any(v > 13 for v in b), [v%7 for v in b]))
