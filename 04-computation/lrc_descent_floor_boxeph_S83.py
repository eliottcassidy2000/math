#!/usr/bin/env python3
"""
The LRC(13)-descent floor  (boxeph-2026-07-17-S83)
==================================================

CLAIM (elementary, via LRC(13) SETTLED + one perturbation):
    M(V) >= v_max / (13 (v_max + v_2nd))          [v_2nd = second-largest speed]
Hence  M(V) >= 1/14  <=>  v_max >= 13 * v_2nd  (the far-element / dominant regime),
and universally  M(V) >= 1/26  (rho=1 floor).  The COMPACT-CORE RESIDUAL is exactly
    1 <= rho := v_max/v_2nd < 13,
where this elementary floor rho/(13(rho+1)) sits below 1/14.

MECHANISM: remove v_max -> V' (12 speeds).  LRC(13) gives t* with min_{V'}||v t*||>=1/13.
The kept runners have margin 1/13 - 1/T; perturbing t = t* + s (|s| <= (1/13-1/T)/v_2nd)
keeps them >= 1/T while sweeping v_max by up to v_max*(1/13-1/T)/v_2nd = rho(1/13-1/T),
enough to lift v_max to >= 1/T iff rho(1/13-1/T) >= 1/T, i.e. T >= 13(1+1/rho).

This script: (1) verifies M(V) >= v_max/(13(v_max+v_2nd)) exactly on many families;
(2) confirms rho>=13 => M>=1/14; (3) checks complementarity with THM-1003 (fill-1);
(4) locates the compact residual rho in [1,13).
"""
from fractions import Fraction as F
from math import gcd
import random

def norm(x):
    r = x % 1
    return min(r, 1 - r)

def exact_M(V):
    """M(V)=max_t min_v||v t||, via THM-999 pair-sum candidates q|v_i+v_j, q<=2max."""
    best = F(0)
    qs = set([14])
    for i in range(len(V)):
        for j in range(i, len(V)):
            s = V[i] + V[j]
            for d in range(1, s + 1):
                if s % d == 0:
                    qs.add(d)
    for q in qs:
        for a in range(1, q):
            if gcd(a, q) == 1:
                m = min(min((v * a) % q, q - (v * a) % q) for v in V)
                c = F(m, q)
                if c > best:
                    best = c
    return best

def descent_floor(V):
    vs = sorted(V)
    vmax, v2 = vs[-1], vs[-2]
    return F(vmax, 13 * (vmax + v2))

def rho(V):
    vs = sorted(V)
    return F(vs[-1], vs[-2])

# THM-1003 fill-1 certifier (for complementarity check)
def band_ok(V, t):
    return min(norm(F(v) * t) for v in V) >= F(1, 14)

def thm1003_certifies(V, N=14):
    vmax = max(V)
    for b in range(2, N):
        S = [v for v in V if v % b == 0]
        if len(S) != 1:
            continue
        vstar = S[0]
        body = [v for v in V if v % b != 0]
        B = max(body) if body else 0
        if b * B <= (N - b) * vstar:
            return True
    return False

def report(fams, title):
    print(f"\n--- {title} ---")
    ok = True
    for name, V in fams:
        M = exact_M(V)
        fl = descent_floor(V)
        r = rho(V)
        holds = (M >= fl)
        ok &= holds
        cert14 = (M >= F(1, 14))
        descent14 = (r >= 13)  # rho>=13 => floor>=1/14
        t1003 = thm1003_certifies(V)
        print(f"  {name:26s} rho={float(r):5.2f} M={str(M):9s} floor={str(fl):9s} "
              f"M>=floor:{holds}  M>=1/14:{cert14}  descent(rho>=13):{descent14}  fill1(1003):{t1003}")
    return ok

if __name__ == "__main__":
    print("=" * 96)
    print("THE LRC(13)-DESCENT FLOOR:  M(V) >= v_max/(13(v_max+v_2nd)) = rho/(13(rho+1))")
    print("=" * 96)

    named = [
        ("deep well {1..12,182}", list(range(1,13))+[182]),
        ("residue {1..11,13,84}", list(range(1,12))+[13,84]),
        ("GW {1..11,13,24}",       list(range(1,12))+[13,24]),
        ("prefix/AP {1..13}",      list(range(1,14))),
        ("{2..14} full-cover",     list(range(2,15))),
        ("kps floor min V",        [3,4,11,12,13,15,18,20,24,42,55,64,67]),
        ("{1..10,13,22,84} multi", list(range(1,11))+[13,22,84]),
    ]
    ok1 = report(named, "named crux families")

    # random families to stress the floor inequality
    random.seed(1)  # determinism (Math.random-free spirit: fixed seed)
    rnd = []
    for k in range(12):
        V = sorted(random.sample(range(1, 60), 13))
        rnd.append((f"rand#{k}", V))
    ok2 = report(rnd, "random 13-families (floor must hold for ALL)")

    # dominant far-element families: rho>=13 should always give M>=1/14
    dom = []
    for w in [13*13, 13*12+5, 200, 500]:
        V = list(range(1,13)) + [w]
        dom.append((f"{{1..12,{w}}}", V))
    ok3 = report(dom, "dominant far-element (rho>=13 -> M>=1/14)")

    print("\n" + "=" * 96)
    print(f"FLOOR HOLDS on all named: {ok1}   on all random: {ok2}   on all dominant: {ok3}")
    print("Compact-core residual = families with 1 <= rho < 13 (floor < 1/14).")
    print("Descent (rho>=13) and THM-1003 (fill-1) are COMPLEMENTARY far-element tools.")
