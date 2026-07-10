#!/usr/bin/env python3
"""
THE RELATION-LATTICE RESUMMATION (klein-2026-07-10-S234, HYP-5895, THM-685).

THE TRANSFER THEOREM (THM-685, elementary):
    LM(q) = #{c in [1,q) : c*v_l mod q in [ceil(q/14), floor(13q/14)] for all l}
          = #{c : c/q in P(S)},   P(S) = {alpha : frac(v_l alpha) in [1/14,13/14] all l}
(EXACT: integer points cannot separate the real band from the rounded band), and
P(S) is a union of at most Sum(v_l) closed intervals (coordinate l enters the
band exactly v_l times per period), each contributing sampling error <= 1:

    |LM(q) - q*mu(S)| <= N(S) <= Sum_l v_l      for EVERY integer q >= 1.

So q > Sum(v)/mu(S) forces LM(q) > 0: any continuum measure floor becomes an
explicit live certificate at every large modulus (prime or not) -- and dually a
dead ruler at q certifies mu(S) <= (N + something)/q... precisely mu <= (LM+N)/q.

EXACT LATTICE CONSTANTS: for reduced direction w (gcd 1), the line measure
A_inf(w) = meas{alpha : frac(w_l alpha) in band, all l} is computed EXACTLY as
a Fraction by breakpoint sweep (breakpoints (14k+1)/(14w_l), (14k+13)/(14w_l)).
Connected constants: c2 = A2 - (6/7)^2;  c3 = A3 - (6/7)*Sum(3 pairs A2) + 2(6/7)^3.
Mobius reconstruction: mu(S) = (6/7)^13 + (6/7)^11*Sum(78 c2) + (6/7)^10*Sum(286 c3)
+ R>=4 (higher connected terms, reported exactly).

Instances: GEN (S233 adversary), DIL (quarantined near-dilation), and THE DEEP
WELL {1..12,182} (the covering-min extremizer M = 14/183 -- maximal relation
stacking; the well's depth vs its lattice is question (iii) of the lead).
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd

GEN = [12, 33, 46, 47, 68, 73, 79, 81, 85, 87, 91, 112, 120]
DIL = [20, 41, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260]
WELL = list(range(1, 13)) + [182]
R = F(6, 7)


def in_band(v, a):
    # frac(v*a) in [1/14, 13/14], exact: 14*r in [den, 13*den] with r = (v*num) % den
    r = (v * a.numerator) % a.denominator
    return a.denominator <= 14 * r <= 13 * a.denominator


def line_sweep(w):
    """Exact line measure + CLOSED-component count (isolated points included).

    P(S) is closed (closed band); an isolated point is a breakpoint whose two
    neighbor gaps are bad but which itself satisfies every closed band
    condition -- a component of measure zero that LM(q) can still hit when
    denominators align (the tight sub-configurations).  Both interval
    components and isolated points obey the +-1 sampling error, so the
    transfer bound uses K = comps + isolated."""
    pts = {F(0), F(1)}
    for v in w:
        for k in range(v):
            pts.add(F(14 * k + 1, 14 * v))
            pts.add(F(14 * k + 13, 14 * v))
    pts = sorted(pts)
    mu = F(0)
    comps = 0
    isolated = 0
    prev_good = False
    for x, y in zip(pts, pts[1:]):
        mid = (x + y) / 2
        good = all(in_band(v, mid) for v in w)
        if good:
            mu += y - x
            if not prev_good:
                comps += 1
        elif not prev_good and x > 0 and all(in_band(v, x) for v in w):
            isolated += 1  # bad-good(point)-bad: an isolated tight point
        prev_good = good
    return mu, comps + isolated


def LM(S, q):
    lo, hi = (q + 13) // 14, (13 * q) // 14
    n = 0
    for c in range(1, q):
        if all(lo <= c * v % q <= hi for v in S):
            n += 1
    return n


def classify(a, c, d):
    if a + c == d:
        return "SCHUR"
    if 2 * c == a + d:
        return "AP"
    ok = all(max(x // gcd(x, y), y // gcd(x, y)) <= 6
             for x, y in ((a, c), (a, d), (c, d)))
    return "RATIO" if ok else "-"


def study(S, name, qs, dead_scan=None):
    print(f"\n================ {name}: {S}")
    sv = sum(S)
    mu, ncomp = line_sweep(S)
    print(f"mu(S) = {mu} = {float(mu):.6f}   components N(S) = {ncomp} "
          f"(bound Sum v = {sv})")
    # ---- transfer verification
    worst = 0
    for q in qs:
        lm = LM(S, q)
        err = lm - q * mu
        worst = max(worst, abs(err))
        tag = "LIVE" if lm else "DEAD"
    print(f"transfer check over {len(qs)} moduli (primes AND composites, "
          f"q = {qs[0]}..{qs[-1]}): max|LM - q*mu| = {float(worst):.1f} "
          f"vs bound N = {ncomp} (slack {ncomp / max(float(worst), 1e-9):.1f}x); "
          f"bound holds: {worst <= ncomp}")
    qstar = None
    if mu > 0:
        qstar = sv // mu + 1  # ceil(Sum v / mu): LM > 0 guaranteed beyond
        qstar = int(qstar)
        print(f"q* = ceil(Sum v / mu) = {qstar}: LM(q) > 0 GUARANTEED for q > q*")
    if dead_scan:
        dead = [q for q in range(15, dead_scan) if LM(S, q) == 0]
        print(f"dead rulers q < {dead_scan}: {len(dead)}"
              + (f", largest = {max(dead)} (vs q* = {qstar})" if dead else ""))
    # ---- exact connected constants
    pairsA = {}
    for i, j in combinations(range(13), 2):
        g = gcd(S[i], S[j])
        pairsA[i, j] = line_sweep([S[i] // g, S[j] // g])[0]
    c2 = {k: a - R ** 2 for k, a in pairsA.items()}
    c3 = {}
    for i, j, k in combinations(range(13), 3):
        g = gcd(gcd(S[i], S[j]), S[k])
        a3 = line_sweep([S[i] // g, S[j] // g, S[k] // g])[0]
        c3[i, j, k] = a3 - R * (pairsA[i, j] + pairsA[i, k] + pairsA[j, k]) \
            + 2 * R ** 3
    top2 = sorted(c2.items(), key=lambda kv: -abs(kv[1]))[:4]
    top3 = sorted(c3.items(), key=lambda kv: -abs(kv[1]))[:6]
    print("top |c2| pairs:  " + "; ".join(
        f"({S[i]},{S[j]}) {float(v):+.5f}" for (i, j), v in top2))
    print("top |c3| triples:")
    for (i, j, k), v in top3:
        print(f"   ({S[i]},{S[j]},{S[k]}) c3 = {float(v):+.6f}  "
              f"{classify(S[i], S[j], S[k])}")
    nrel = sum(1 for (i, j, k) in c3 if classify(S[i], S[j], S[k]) != "-")
    # ---- the Mobius reconstruction
    l2 = R ** 11 * sum(c2.values())
    l3 = R ** 10 * sum(c3.values())
    rest = mu - R ** 13 - l2 - l3
    print(f"RECONSTRUCTION: mu = (6/7)^13 {float(R**13):+.6f}  "
          f"layer2 {float(l2):+.6f}  layer3 {float(l3):+.6f}  "
          f"R>=4 {float(rest):+.6f}   [relation triples: {nrel}/286]")
    print(f"   deficit mu - (6/7)^13 = {float(mu - R**13):+.6f}; "
          f"share explained by t<=3: "
          f"{float((l2 + l3) / (mu - R**13)) * 100 if mu != R**13 else 0:.1f}%")
    return dict(mu=mu, l2=l2, l3=l3, rest=rest, c3=c3, S=S)


qs_test = [139, 197, 239, 353, 383, 500, 547, 1000, 1009, 1600, 2003,
           2048, 3001, 4096, 5003, 7000, 9973, 12000, 20011]
g = study(GEN, "GEN (generic adversary)", qs_test)
d = study(DIL, "DIL (quarantined near-dilation)", qs_test)
w = study(WELL, "DEEP WELL {1..12,182} (covering-min extremizer, M = 14/183)",
          qs_test, dead_scan=2000)

print("\n================ CROSS-CHECKS")
key = tuple(sorted([GEN.index(12), GEN.index(73), GEN.index(85)]))
print(f"GEN (12,73,85): exact line c3 = {float(g['c3'][key]):+.6f} "
      f"vs S233 measured dev3/q at q=5003 = -0.012679 "
      f"vs subtorus -17/1372 = {float(F(-17,1372)):+.6f}")
print(f"GEN continuum layers (l2, l3) = ({float(g['l2']):+.4f}, {float(g['l3']):+.4f}) "
      f"vs S233 finite-q signed (L2sig, L3sig) at q=5003 = (-0.0023, -0.0647)")
print(f"WELL sanity: M(WELL) = 14/183 > 1/14 strictly => mu > 0: mu = {w['mu']} "
      f"= {float(w['mu']):.8f}")
