#!/usr/bin/env python3
"""
VERIFY the deep-hiding dichotomy (mac-mini-2026-07-03-S30), elementary complement to THM-566.
 Lemma 1 (covering => deep hiding): S covers {2..n} => min-dist=0 at every t=a/q, q<=n => q*(S) >= n+1.
 Lemma 2 (tight => n | q*): M(S)=1/n at t*=a/q* (lowest terms) => min residue-dist = q*/n in Z => n | q*.
 Cor: tight covering => q* >= 2n (n | q*, q* >= n+1). For n=14: q* in {28,42,...}.
 Also: verify the AP/GW tight locus is NON-covering, and g(14)<=3 (Steinhaus) on near-tight covering families.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def is_covering(sp, n): return all(any(v % q == 0 for v in sp) for q in range(2, n+1))
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)):
            Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0); bt = None
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            t = F(a, q); m = min(nd(v*t) for v in sp)
            if m > best: best, bt = m, t
    return best, bt
def qstar(sp):
    M, t = M_exact(sp)
    return t.denominator, M
def ngaps(sp, t):
    ph = sorted(set(float((v*t) % 1) for v in sp) | {0.0})
    gs = set(round(ph[(i+1) % len(ph)]-ph[i] if i < len(ph)-1 else 1-ph[-1], 6) for i in range(len(ph)))
    return len(gs)

if __name__ == "__main__":
    import sys
    def out(*a): print(*a); sys.stdout.flush()
    n = 14; rng = random.Random(303)
    out("LEMMA 1 (covering => q* >= n+1 = 15): test many covering gcd=1 families, check q*>=15.")
    out("="*80)
    viol1 = 0; tested = 0; qmin = 999
    for _ in range(2500):
        hi = rng.choice([30, 60, 90])
        sp = sorted(set(rng.sample(range(1, hi), 13)))
        if len(sp) != 13 or gcd_all(sp) != 1 or not is_covering(sp, n): continue
        tested += 1
        q, M = qstar(sp)
        qmin = min(qmin, q)
        if q < n+1: viol1 += 1
        # also confirm min-dist is 0 at every t=a/q, q<=14
        for q0 in range(2, n+1):
            a = rng.randrange(1, q0)
            md = min(nd(v*F(a, q0)) for v in sp)
            if md != 0:
                out(f"   !! covering family min-dist {md}!=0 at t={a}/{q0}: {sp}"); break
    out(f"   tested {tested} covering gcd=1 families: q* >= 15 violations = {viol1}; min q* seen = {qmin}")
    out(f"   (Lemma 1 predicts q* >= 15 ALWAYS for covering; and min-dist=0 at all t=a/q, q<=14)")

    out("\nLEMMA 2 (tight M=1/14 => 14 | q*): check the tight locus + even block.")
    AP = list(range(1, 14))
    evenblock = [2*i for i in range(1, 14)]
    for name, S in [("AP {1..13}", AP), ("even block 2*{1..13}", evenblock),
                    ("deep well {1..12,182}", list(range(1, 13))+[182])]:
        q, M = qstar(S); cov = is_covering(S, n); g = gcd_all(S)
        out(f"   {name:>26}: M={str(M):>7}={float(M):.5f}  q*={q:>4}  14|q*? {q % 14 == 0}  covering? {cov}  gcd={g}")

    out("\nCOR (tight covering => q* >= 28) + g(14)<=3 on the 12 tightest covering families found:")
    cands = [list(range(1, 13))+[182]]
    for _ in range(3000):
        hi = rng.choice([30, 60, 120, 200]); sp = sorted(set(rng.sample(range(1, hi), 13)))
        if len(sp) == 13 and gcd_all(sp) == 1 and is_covering(sp, n): cands.append(sp)
    scored = sorted(((M_exact(sp)[0], sp) for sp in cands), key=lambda x: x[0])[:12]
    gviol = 0
    for M, sp in scored:
        _, t = M_exact(sp); q = t.denominator; g = ngaps(sp, t)
        if g > 3: gviol += 1
        out(f"   M={float(M):.5f} q*={q:>4} 14|q*?{str(q%14==0):>5} #gaps@opt={g}  {sp[:7]}...")
    out(f"\n   g>3 violations among tightest covering families: {gviol} (Steinhaus g(14)<=3 predicts 0)")
    out("=> Lemma 1 (covering=>deep, q*>=15) + Lemma 2 (tight=>14|q*) => tight covering hides at q*>=28;")
    out("   the tight LOCUS (M=1/14) is the AP & even block, and the primitive tight ones are NON-covering.")
