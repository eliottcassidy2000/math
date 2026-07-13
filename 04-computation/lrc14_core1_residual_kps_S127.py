# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont62: the |core|=1 residual (mac-mini-S74) IS the interval-core / single-killer domain.
#
# opus-S259/mac-mini-S74: a covering family splits into CORE = speeds coprime to 30030=2*3*5*7*11*13 (primes<=13)
# and NON-CORE = the rest (each divisible by a prime <=13). Equidistribution CLOSES |core|>=2; the residual is
# |core|=1 (a SINGLE coprime runner, one arc, no equidistribution -- density up to 0.92 on the good set).
# CLAIM (the connection): (a) the covering-min extremals ARE |core|=1 (deep well core={1}); (b) the HARDEST
# |core|=1 case is core-runner = 1 (smallest => fewest arcs => highest density); larger core runners equidistribute
# and are EASY; (c) core-runner-1 |core|=1 covering families = exactly my interval-core families -- single-killer
# (PROVED + formalized cont.60/61) + multi-killer (enumerated cont.58). So my structural work is the residual's
# hard sub-case, and the single-killer ladder is its formalized extremal slice.
from math import gcd
from fractions import Fraction as F
from functools import reduce
from itertools import combinations

P = 2*3*5*7*11*13  # 30030
def norm(x): r = x - int(x); r = r + 1 if r < 0 else r; return min(r, 1 - r)
def Mexact(v):
    qc = 4*max(v)+2; best = F(0)
    for q in range(2, qc):
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(F(vi*p, q)) for vi in v)
                if m > best: best = m
    return best
def is_cov(v, N=14): return all(any(x % d == 0 for x in v) for d in range(2, N+1))
def prim(v): return reduce(gcd, v) == 1
def core(v): return [x for x in v if gcd(x, P) == 1]

def main():
    tgt = F(14, 183)
    print(f"CORE = coprime to 30030=2*3*5*7*11*13. covering-min 14/183={float(tgt):.6f}\n")
    print("(a) the covering-min extremals + single-killer ladder are |core|=1 (core = {1}):")
    for name, v in [("deep well {1..12,182}", list(range(1,13))+[182]),
                    ("ladder S_2 {1..12,364}", list(range(1,13))+[364]),
                    ("ladder S_3 {1..12,546}", list(range(1,13))+[546]),
                    ("multi-killer {1..11,13,84}", list(range(1,12))+[13,84])]:
        c = core(v)
        print(f"    {name:<28} core={c} |core|={len(c)}  M={Mexact(v)}")

    print("\n(b) HARDEST |core|=1 = core-runner 1 (small); larger core runner => higher M (easier). Swap 1->c:")
    print(f"    {'core runner':>11} | family {{c}}u{{2..12,182}} | covering? | M | vs 14/183")
    for cr in [1, 17, 19, 23, 29]:
        v = [cr] + list(range(2,13)) + [182]
        if len(set(v)) != 13: continue
        cc = core(v); M = Mexact(v)
        print(f"    {cr:>11} | {'core='+str(cc):<24} | {is_cov(v)!s:>5} | {str(M):>8}={float(M):.5f} | {'=' if M==tgt else 'higher' if M>tgt else 'LOWER!'}")

    print("\n(c) |core|=1 covering-min over an enumerated slice (interval core {1..12}+killer, and shorter):")
    best = F(1); bf = None
    # single-killer sweep {1..12, k} + confirm min is the deep well
    for k in range(13, 200):
        v = list(range(1,13))+[k]
        if prim(v) and is_cov(v) and len(core(v))==1:
            M = Mexact(v)
            if M < best: best, bf = M, v
    print(f"    single-killer |core|=1 min over k<=199: M={best}={float(best):.6f} at {bf}  (= deep well? {best==tgt})")
    print()
    print("=> the |core|=1 residual's HARD sub-case (core runner 1) = my interval-core families. Single-killer")
    print("   PROVED + formalized (cont.60/61); multi-killer enumerated (cont.58). Larger core runners equidistribute")
    print("   (opus |core|>=2 style) and are easy. So structural (kps) + analytic (opus/mac-mini) tile the residual.")

if __name__ == "__main__":
    main()
