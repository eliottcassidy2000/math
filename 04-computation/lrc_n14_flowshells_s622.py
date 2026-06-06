#!/usr/bin/env python3
"""
S622 — LRC n=14, C'(14): twisted involutions on flow shells.
Convention (THM-398): n=14 runners, 13 speeds, gap 1/n=1/14, M(S)=max_t min ||v t||.
C'(14): if 14|v for some speed, then M(S) > 1/14 (loose).  We probe the ALL-SHORT
residual (THM-398 §4) — the configs where local perturbation at a 14-clock witness
cannot clear v's arc — and analyse WHERE the escape point t* lives, in the
CRT grid Z/14 = Z/2 x Z/7 and the 2n-1 = 27 = 3^3 shells.
"""
from fractions import Fraction as Fr
from math import gcd
import itertools, sys
sys.path.insert(0, '04-computation')
from lrc_tight_enum_s621 import norm           # ||x|| for Fraction

def gap_and_argmax(speeds):
    """exact M(S)=max_t min_i||v_i t|| and the set of argmax times t* in (0,1/2]."""
    V = [abs(v) for v in speeds]; cands = set()
    for i in range(len(V)):
        vi = V[i]
        for k in range(0, 2*vi+1):
            t = Fr(2*k+1, 2*vi)
            if 0 < t <= Fr(1, 2): cands.add(t)
        for j in range(i):
            vj = V[j]
            for d in (vi+vj, abs(vi-vj)):
                if d == 0: continue
                kk = 1
                while Fr(kk, d) <= Fr(1, 2):
                    cands.add(Fr(kk, d)); kk += 1
    best = Fr(0); arg = []
    for t in cands:
        m = min(norm(v*t) for v in V)
        if m > best: best, arg = m, [t]
        elif m == best: arg.append(t)
    return best, sorted(arg)

def safe_components(Sp, level):
    """connected components (as exact intervals) of G = {t in [0,1): ||v t||>level for all v in Sp}.
       returns list of (a,b) lengths via the arc-endpoint arrangement."""
    bps = {Fr(0), Fr(1)}
    for v in Sp:
        for k in range(v):
            for s in (level, -level):
                t = (Fr(k)+s)/v; t -= t.numerator//t.denominator
                if 0 <= t < 1: bps.add(t)
    bps = sorted(bps); comps = []; cur = None
    for a, b in zip(bps, bps[1:]):
        mid = (a+b)/2
        safe = all(norm(v*mid) > level for v in Sp)
        if safe:
            cur = (cur[0], b) if cur else (a, b)
        else:
            if cur: comps.append(cur); cur = None
    if cur: comps.append(cur)
    return comps

def analyse(S):
    n = 14; lvl = Fr(1, n)
    S = sorted(S)
    mults = [v for v in S if v % n == 0]
    M, arg = gap_and_argmax(S)
    loose = M > lvl
    # residual test: for the (first) multiple v=nw, are all components of G(S\{v}) short (<=2/(n^2 w))?
    info = ""
    if mults:
        v = mults[0]; w = v//n; Sp = [x for x in S if x != v]
        comps = safe_components(Sp, lvl)
        rho2 = Fr(2, n*v)                       # v's arc diameter 2/(n v)
        longest = max((b-a for a,b in comps), default=Fr(0))
        allshort = longest <= rho2
        info = f" v={v}(w={w}) Bprime_long={longest} vs 2/nv={rho2} all_short={allshort}"
    return M, loose, arg, info, mults

if __name__ == "__main__":
    n = 14
    print(f"C'(14): multiple-of-14 configs must be loose (M>1/14={float(Fr(1,14)):.5f}).")
    print("Probing near-tight worry configs (AP with one speed -> a multiple of 14).\n")
    base = list(range(1, 14))               # AP {1..13}, the tight set at n=14
    # replace each AP element by 14, 28, 42 (w=1,2,3): natural worry configs
    seen = set()
    for w in (1, 2, 3):
        v = n*w
        for i in range(13):
            S = base[:]; S[i] = v
            if len(set(S)) < 13: continue
            g = 0
            for x in S: g = gcd(g, x)
            if g != 1: continue
            key = tuple(sorted(S))
            if key in seen: continue
            seen.add(key)
            M, loose, arg, info, mults = analyse(S)
            tag = "LOOSE" if loose else "**TIGHT**"
            print(f"AP[{i}={base[i]}]->{v}: M={M}={float(M):.5f} {tag} t*={[str(a) for a in arg]}{info}")
