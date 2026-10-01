#!/usr/bin/env python3
"""procgen_mlr_20260930_lrc.py -- Part A (LRC side, existence) of the multiplicative lonely runner lane.

A1  kappa of the boxes {2^j 3^k : j<J, k<K} (continuous time): exact values, the 1/5 rigidity set.
A2  the x2x3 lonely spectrum I(t) = inf_{j,k>=0} ||2^j 3^k t||: computer-assisted certificate of its
    top six values (1/5, 1/7, 1/10, 1/11, 1/13, 1/14) + finite-exact census over small denominators.
A3  discrete time t = h/D: maximal loneliness L(D,J,K) for primes, composite D and Collatz gates;
    comparison with 1/(n+1), the Poisson-apex model, and the spectral corollary (small prime factors).
A4  discrete multiplicative LRC via the 6-free counting criterion (PROVED bound F_{J,K}) + the exact
    exception list for small boxes.
Prints to stdout; raises on any failed check.
"""
import math
import sys
import time
from fractions import Fraction as F

import numpy as np

import procgen_mlr_20260930_core as core
from procgen_mlr_20260930_core import check

T0 = time.time()


def hdr(s):
    print("\n" + "=" * 100 + "\n" + s + "\n" + "=" * 100)
    sys.stdout.flush()


# ----------------------------------------------------------------------------------------------
def a1_kappa_boxes():
    hdr("A1. kappa(J,K) = sup_t min_{j<J,k<K} ||2^j 3^k t||  (continuous time, exact)")
    pred = lambda J, K: F(1, 2) if J == 1 else (F(1, 3) if K == 1 else (F(1, 4) if J == 2 else F(1, 5)))
    rows = []
    for J in range(1, 7):
        for K in range(1, 5):
            V = [2 ** j * 3 ** k for j in range(J) for k in range(K)]
            if max(V) > 2000:
                continue
            kap, arg = core.kappa_exact(V)
            check(kap == pred(J, K), f"kappa({J},{K}) = {kap} != predicted {pred(J, K)}")
            rows.append((J, K, len(V), kap, arg))
    print("  J K  n   kappa   argmax     kappa*(n+1)   [all equal the PROVED table 1/2,1/3,1/4,1/5]")
    for J, K, n, kap, arg in rows:
        print(f"  {J} {K} {n:2d}   {str(kap):5s}   {str(arg):6s}   {float(kap * (n + 1)):6.2f}")
    # rigidity: for J>=3, K>=2 the set {t : min ||vt|| >= 1/5} is exactly {1,2,3,4}/5
    four = [(F(r, 5), F(r, 5)) for r in range(1, 5)]
    check(core.X_set([1, 2, 3, 4], F(1, 5)) == four, "X_{1/5}({1,2,3,4}) != {r/5}")
    nbox = 0
    for J in range(3, 17):
        for K in range(2, 11):
            X = core.X_set([2 ** j * 3 ** k for j in range(J) for k in range(K)], F(1, 5))
            check(X == four, f"X_1/5 box({J},{K}) != {{r/5}}")
            nbox += 1
    print(f"  X_(1/5)(box) = {{1/5,2/5,3/5,4/5}} exactly for all {nbox} boxes 3<=J<=16, 2<=K<=10 "
          f"(exact rational intervals); also for the AP {{1,2,3,4}} itself.")
    print("  => PROVED: kappa(V) = 1/5 for every finite 3-smooth V containing a dilate of {1,2,3,4};")
    print("     LRC margin kappa*(n+1) = (n+1)/5 -> infinity.")


# ----------------------------------------------------------------------------------------------
def a2_spectrum():
    hdr("A2. the x2x3 lonely spectrum I(t) = inf_{j,k>=0} ||2^j 3^k t||  (Furstenberg: I(t)>0 only for rationals)")
    # finite-exact census over rationals with denominator <= NMAX
    NMAX = 1200
    vals = {}
    for N in range(2, NMAX + 1):
        for r in range(1, N // 2 + 1):
            if math.gcd(r, N) == 1:
                v = core.I_rat(F(r, N))
                if v >= F(1, 25):
                    vals.setdefault(v, []).append(F(r, N))
    top = sorted(vals, reverse=True)
    print(f"  FINITE-EXACT census, all r/N with N <= {NMAX}: largest values of I (with smallest attaining t):")
    for v in top[:16]:
        ex = sorted(vals[v], key=lambda t: (t.denominator, t))[:4]
        print(f"    I = {str(v):7s} = {float(v):.5f}   e.g. t = {', '.join(map(str, ex))}   (#t in [0,1/2]: {len(vals[v])})")
    # certificates
    print("\n  Computer-assisted certificates (box (16,10) + local covers, exact Fractions):")
    levels = [F(1, 7), F(1, 10), F(1, 11), F(1, 13), F(1, 14)]
    worst_all = 0
    for d in levels:
        t = time.time()
        res = core.certify_spectrum(d, 16, 10)
        check(res["ok"], f"spectrum certificate failed at level {d}: {res['fails'][:3]} {res['unassigned'][:3]}")
        exc_vals = sorted({v for _, v in res["exceptional"]}, reverse=True)
        census_vals = [v for v in top if v >= d]
        check(exc_vals == census_vals, f"certified exceptional values {exc_vals} != census {census_vals}")
        worst_all = max(worst_all, res["worst_multiplier"])
        print(f"    level {str(d):4s}: X has {res['ncomp']:3d} components -> {res['ncentre']:3d} centres; "
              f"exceptional t (I >= level): {len(res['exceptional']):3d}; values {', '.join(map(str, exc_vals))}; "
              f"max multiplier factor {float(res['worst_multiplier']):.3g}; {time.time() - t:.1f}s")
        last = res
    dens = sorted({t.denominator for t, _ in last["exceptional"]})
    odd6 = [N for N in dens if N % 2 and N % 3]
    check(odd6 == [5, 7, 11, 13], f"odd 6-free exceptional denominators {odd6}")
    print(f"  denominators of the 88 exceptional t at level 1/14: {dens}; coprime to 6: {odd6}")
    nd = sorted({t.denominator for t, _ in last["nonexceptional"]})
    check(all(N % 2 == 0 or N % 3 == 0 for N in nd), "a non-exceptional centre has a denominator coprime to 6")
    print(f"  non-exceptional centres (in X but I < 1/14): {len(last['nonexceptional'])}, denominators {nd} "
          f"(all divisible by 2 or 3: none is of the form h/D with D coprime to 6)")
    print("  => PROVED: {I(t) : t in R} intersected with [1/14, 1/2] = {1/5, 1/7, 1/10, 1/11, 1/13, 1/14};")
    print("     every t outside the 88 certified rationals has I(t) < 1/14.  Next census values: 5/73, 7/104, 1/15.")
    return worst_all


# ----------------------------------------------------------------------------------------------
def poisson_apex_pred(D, J, K):
    A = (J - 1) * (K - 1) / 3 + 2 * (K - 1) / 3 + (J - 1) / 2 + 1
    return math.log(D) / (2 * A)


def a3_discrete(worst):
    hdr("A3. discrete time h/D: maximal loneliness L(D,J,K) = max_h min_box ||2^j 3^k h/D||")
    print("  (i) primes D, boxes J = c ln D/ln 2, K = c ln D/ln 3")
    print("      D          c    J   K     n     L          L*(n+1)  L*n/lnD  indep lnD/2n  apex-model")
    for D in (10007, 100003, 1000003):
        for c in (0.5, 1.0, 2.0, 3.0):
            J = max(1, round(c * math.log(D) / core.LOG2))
            K = max(1, round(c * math.log(D) / core.LOG3))
            n = J * K
            L, h = core.loneliness(D, J, K)
            check(L <= F(1, 5), "L > 1/5 contradicts kappa(box) <= 1/5")
            print(f"      {D:8d}  {c:3.1f}  {J:3d} {K:3d}  {n:5d}   {float(L):.5f}   {float(L * (n + 1)):6.2f}"
                  f"   {float(L) * n / math.log(D):6.3f}   {math.log(D) / (2 * n):.5f}     {poisson_apex_pred(D, J, K):.5f}")
            sys.stdout.flush()
    # (ii) spectral corollary: D with a small 'lonely' prime factor
    print("\n  (ii) spectral corollary check: L(D,J,K) for D = q * P (P prime) and boxes containing all 3-smooth")
    print("       numbers <= (multiplier factor)*D; prediction L = I(1/q) for the smallest q in {5,7,11,13} dividing D,")
    print("       and L < 1/14 when D has no such factor")
    fac = float(worst)
    for q, I in ((5, F(1, 5)), (7, F(1, 7)), (11, F(1, 11)), (13, F(1, 13)), (1, None)):
        for P in (1009, 20011):
            D = q * P if q > 1 else core.next_prime(P * 17)
            J = math.ceil(math.log2(fac * D)) + 1
            K = math.ceil(math.log(fac * D, 3)) + 1
            L, h = core.loneliness(D, J, K)
            if I is not None:
                check(L == I, f"corollary fails: D={D} L={L} expected {I}")
            else:
                check(L < F(1, 14), f"corollary fails: D={D} prime L={L} >= 1/14")
            print(f"       D = {D:7d} (q={q:2d})  box ({J},{K})  L = {str(L):12s} = {float(L):.5f}  "
                  f"{'= I(1/q)' if I is not None else '< 1/14'}")
    # (iii) Collatz gates with the transport box (J,K) = (p,a)
    print("\n  (iii) Collatz gates D = |2^p - 3^a|, box (J,K) = (p,a), 12 <= p <= 21, 2 <= a < p, D <= 2.2e6")
    print("       p  a        D   spf     n     L          L*(n+1)   note")
    rows = []
    for p in range(12, 22):
        for a in range(2, p):
            D = abs(2 ** p - 3 ** a)
            if D < 50 or D > 2_200_000:
                continue
            n = p * a
            L, h = core.loneliness(D, p, a)
            spf = next(f for f in range(5, D + 1) if D % f == 0)
            note = ""
            if D % 5 == 0:
                check(L == F(1, 5), f"gate {(p, a)}: 5|D but L={L}")
                note = "5|D => L = 1/5 (rigid)"
            rows.append((p, a, D, spf, n, L))
            if a % 4 == 0 or p == 21:
                if True:
                    print(f"      {p:2d} {a:2d} {D:9d} {spf:5d}  {n:4d}   {float(L):.5f}   {float(L * (n + 1)):7.2f}   {note}")
    nf = sum(1 for r in rows if r[5] * (r[4] + 1) < 1)
    print(f"       {len(rows)} gate boxes computed; L*(n+1) < 1 (discrete LRC fails) on {nf} of them; "
          f"5 | D on {sum(1 for r in rows if r[2] % 5 == 0)} (all with L = 1/5).")
    print(f"       min L*(n+1) over all {len(rows)} gate boxes: {min(float(r[5] * (r[4] + 1)) for r in rows):.2f}")
    for q, v in ((7, F(1, 7)), (11, F(1, 11)), (13, F(1, 13)), (73, F(5, 73))):
        hits = [r for r in rows if r[2] % q == 0 and all(r[2] % q2 for q2 in (5, 7, 11, 13) if q2 < q or q == 73)]
        eq = sum(1 for r in hits if r[5] == v)
        print(f"       gates with q = {q:2d} | D (and no smaller lonely prime): {len(hits):2d}, of which L = {v}: {eq}")
    med = np.median([float(r[5] * (r[4] + 1)) for r in rows if r[2] % 5 and r[2] % 7])
    print(f"       median L*(n+1) over gates with 5,7 not dividing D: {med:.2f}")


# ----------------------------------------------------------------------------------------------
def H_bound(J, K, R):
    """smooth majorant of the lower-order part of the PROVED bound for F_{J,K}(R):
    F_{J,K}(R) <= (2R/3)(JK + K + J/2 + 1/2) + 4 (JK + K log2 R + J log3 R + sigma3(R)),
    sigma3(R) <= (log2 R + 1)(log3 R + 1)."""
    l2, l3 = math.log2(R), math.log(R, 3)
    return 4 * (J * K + K * l2 + J * l3 + (l2 + 1) * (l3 + 1))


def H_deriv(J, K, R):
    l2, l3 = math.log2(R), math.log(R, 3)
    return 4 * (K / (R * core.LOG2) + J / (R * core.LOG3) + (l3 + 1) / (R * core.LOG2) + (l2 + 1) / (R * core.LOG3))


def a4_counting():
    hdr("A4. discrete multiplicative LRC: the 6-free counting criterion")
    # (i) exact check of the PROVED bound #bad <= F_{J,K}(R) on instances
    nchk = 0
    worst = 0.0
    for D in (1009, 4001, 10007, 30011):
        for (J, K) in ((3, 2), (4, 3), (6, 4), (9, 6), (12, 8)):
            for R in (2, 5, 17, 60, 200):
                if R > D // 6:
                    continue
                exact = core.bad_set_exact(D, J, K, R)
                Fb = core.bad_set_bound(J, K, R)
                check(exact <= Fb, f"bad-set bound violated D={D} J={J} K={K} R={R}: {exact} > {Fb}")
                worst = max(worst, exact / Fb)
                nchk += 1
    print(f"  (i) #bad(D,J,K,R) <= F_JK(R) verified exactly on {nchk} instances (max ratio exact/F = {worst:.3f})")
    # check the explicit majorant of F on a range of R
    for (J, K) in ((3, 2), (5, 3), (8, 5)):
        for R in (2, 7, 30, 100, 400, 1500):
            Fb = core.bad_set_bound(J, K, R)
            maj = 2 * R / 3 * (J * K + K + J / 2 + 0.5) + 4 * (J * K + K * math.log2(R) + J * math.log(R, 3) + core.sigma3(R))
            check(Fb <= maj, f"majorant of F fails ({J},{K},{R}): {Fb} > {maj}")
            check(core.sigma3(R) <= (math.log2(R) + 1) * (math.log(R, 3) + 1), "sigma3 majorant fails")
    print("      explicit majorant F_JK(R) <= (2R/3)(JK+K+J/2+1/2) + 4(JK + K log2R + J log3R + sigma3(R)) checked")
    # (ii) asymptotics and explicit D0 from the analytic bound
    print("  (ii) certified loneliness from F_JK(ceil(D/(n+1))) < D-1:  L >= 1/(n+1) for all D >= D0 (explicit);")
    print("       asymptotically L >= (3/2)/(n+K+J/2+1/2) (union bound alone: ~1/(2n))")
    print("       J  K    n   (3/2)(n+1)/(n+K+J/2+1/2)   R0     D0")
    D0s = {}
    for J in range(3, 9):
        for K in range(2, 6):
            n = J * K
            c = (n + 1) - 2 / 3 * (n + K + J / 2 + 0.5)
            check(c > 0, "slope condition (J-2)(K-1)>0")
            R0 = 3
            while not (c * R0 - (n + 1) - H_bound(J, K, R0) > 1e-6 * R0 and H_deriv(J, K, R0) < c):
                R0 += max(1, R0 // 50)
            D0 = (n + 1) * (R0 - 1) + 1
            D0s[(J, K)] = D0
            if K <= 3 or J in (3, 8):
                print(f"       {J}  {K}  {n:3d}        {1.5 * (n + 1) / (n + K + J / 2 + 0.5):.3f}            {R0:6d}  {D0}")
    # (iii) below D0: exact criterion, then exact L where it fails; list exceptions to L >= 1/(n+1)
    print("  (iii) all D coprime to 6 with 5 <= D < min(D0, 20000): exact criterion F_JK(ceil(D/(n+1))) < D-1,")
    print("        else exact L(D,J,K); exceptions = D with L(D,J,K) < 1/(n+1)")
    for (J, K) in ((3, 2), (4, 2), (3, 3), (4, 3), (5, 3), (6, 4)):
        n = J * K
        D0 = D0s[(J, K)]
        cap = min(D0, 20000)
        exc = []
        nexact = 0
        for D in range(5, cap):
            if D % 2 == 0 or D % 3 == 0:
                continue
            R = -(-D // (n + 1))
            if Fcache(J, K, R) < D - 1:
                continue
            nexact += 1
            L, _ = core.loneliness(D, J, K)
            if L * (n + 1) < 1:
                exc.append(D)
        print(f"        box ({J},{K}) n={n:2d}: D0={D0}, checked D < {cap}: exact L needed for {nexact} D; "
              f"exceptions ({len(exc)}): {exc[:20]}{' ...' if len(exc) > 20 else ''}")
        sys.stdout.flush()


_FC = {}


def Fcache(J, K, R):
    if (J, K, R) not in _FC:
        _FC[(J, K, R)] = core.bad_set_bound(J, K, R)
    return _FC[(J, K, R)]


def main():
    core.lean_malloc()
    a1_kappa_boxes()
    worst = a2_spectrum()
    a3_discrete(worst)
    a4_counting()
    print(f"\n[Part A done in {time.time() - T0:.0f}s, peak RSS {core.mem_mb():.0f} MB]")
    check(core.mem_mb() < 500, "memory budget exceeded")
    print("PART A CHECKS PASSED")


if __name__ == "__main__":
    main()
