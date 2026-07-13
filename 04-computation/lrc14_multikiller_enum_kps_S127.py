# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont58: the MULTI-KILLER primitive covering families -- the open frontier of the
# covering-min lower bound (opus-S253 Picture 3; single-killer closed by opus + kps cont.57 + mac-mini S68).
#
# A single-killer family = interval core {1..12} + ONE outlier (the deep well {1..12,182}). MULTI-KILLER =
# a shorter interval core {1..k}, k<12, plus 2+ outliers supplying the missing divisors {k+1,..,12,13,14}.
# QUESTION: does any multi-killer primitive DC family beat the single-killer floor 14/183? Enumerate exactly.
from math import gcd
from fractions import Fraction as F
from functools import reduce
from itertools import combinations

def norm(x): r = x - int(x); r = r + 1 if r < 0 else r; return min(r, 1 - r)
def Mfloat(v):
    qcap = 3 * max(v) + 2; best = 0.0
    for q in range(2, qcap):
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(vi * p / q) for vi in v)
                if m > best: best = m
    return best
def Mexact(v):
    qcap = 4 * max(v) + 2; best = F(0); bt = None
    for q in range(2, qcap):
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(F(vi * p, q)) for vi in v)
                if m > best: best, bt = m, F(p, q)
    return best, bt
def is_cov(v, N=14): return all(any(x % d == 0 for x in v) for d in range(2, N + 1))
def prim(v): return reduce(gcd, v) == 1

def enum(kcore, nout, Bout):
    """interval core {1..kcore} + nout outliers in (kcore, Bout], primitive DC. kcore+nout=13."""
    assert kcore + nout == 13
    core = list(range(1, kcore + 1))
    tgt = F(14, 183)
    best = 1.0; cands = []
    cnt = 0
    for outs in combinations(range(kcore + 1, Bout + 1), nout):
        v = core + list(outs)
        if not prim(v): continue
        if not is_cov(v): continue
        cnt += 1
        mf = Mfloat(v)
        if mf < best + 1e-9:
            if mf < best - 1e-9: cands = []
            best = mf; cands.append(v)
    # exact-verify minimizers
    bestE = F(1); bf = None; bt = None
    for v in cands[:40]:
        me, t = Mexact(v)
        if me < bestE: bestE, bf, bt = me, v, t
    return cnt, bestE, bf, bt, best

def main():
    tgt = F(14, 183)
    print(f"single-killer floor (deep well {{1..12,182}}) = 14/183 = {float(tgt):.6f}\n")
    print(f"{'structure':<28} | #prim-DC | {'min M exact':>12} = float   | at t* | minimizer | vs 14/183")
    for (k, no, B) in [(12,1,200),(11,2,200),(10,3,120),(9,4,80)]:
        cnt, mE, bf, bt, mf = enum(k, no, B)
        cmp = "BELOW!" if mE < tgt else ("= (tie)" if mE == tgt else "above")
        label = f"{{1..{k}}}+{no} out (<= {B})"
        print(f"{label:<28} | {cnt:>7} | {str(mE):>12} = {float(mE):.6f} | {bt} | {bf} | {cmp}")
    print()
    print("=> if all multi-killer minima are >= 14/183 (single-killer strictly the global min), the deep well")
    print("   is the UNIQUE covering-min and multi-killer never binds; that is the structural claim to prove.")

if __name__ == "__main__":
    main()
