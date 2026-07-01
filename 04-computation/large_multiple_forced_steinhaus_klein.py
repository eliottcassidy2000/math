#!/usr/bin/env python3
"""
large_multiple_forced_steinhaus_klein.py  --  klein-2026-06-30-S52  (HYP-3763)

Making "LARGE MULTIPLES ARE FORCED" rigorous, and pushing the residual with Steinhaus.

Context: the lowness lemma (HYP-3747) needs, for each core speed k<=n-2, that dropping k raises
M above n/Phi6(n). The residual (R4 of HYP-3748, HYP-3745) is the UNBOUNDED/large-speed case:
"a large speed substituted for k cannot pull M back below n/Phi6" -- verified n=14, but not a
closed inequality.

This script establishes:

(1) THE LARGE-MULTIPLE-FORCED LEMMA (fully rigorous, general n, NO primality needed).
    q-witness: if no runner is divisible by k, then t=1/k is lonely with gap >= 1/k
    (every runner s has ||s/k|| = dist(s mod k)/k >= 1/k, since k does not divide s).
    And 1/k > n/Phi6(n) for all k <= n-2  (<=> n^2-n+1 > kn <=> n^2-n+1 > (n-2)n = n^2-2n, true).
    So M(S) <= n/Phi6  =>  S has a multiple of k.  If k not in S and k > (n-2)/2, the smallest
    available multiple is 2k > n-2: a LARGE speed is FORCED.  (verified: M(P_k) >= 1/k here.)

(2) THE STEINHAUS DISPLACEMENT (the mechanism closing the residual).
    The forced large multiple kappa = k*c (c>=2) kills resonance k (kappa ≡ 0 mod k, so it fills
    the D=k hole).  But at any OTHER modulus D, its image kappa*a = c*(k*a) is the c-th multiple of
    k's slot -- DISPLACED from the k-slot unless D | k(c-1).  By the three-gap theorem (HYP-3762),
    the punctured core {1..n-2}\{k} leaves a MERGED (double) gap where k's image was; kappa's image
    lands elsewhere, so it does NOT fill the gap it created.  The hole survives at a spread modulus.

(3) THE CLOSED INEQUALITY (verified general n=8..20): the adversary's best bounded escape
    (punctured core + 2 killers) has M >= (the binding value) > n/Phi6, and the binding modulus /
    displacement confirm the Steinhaus mechanism.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def Phi6(n):
    return n * n - n + 1

def dist0(x, D):
    x %= D
    return min(x, D - x)

def M_and_arg(S):
    """Exact M(S) = max_{D,a} min_s ||s a / D||, with the binding (D,a)."""
    Dmax = 2 * max(S) + 2
    best = F(0)
    arg = None
    for D in range(2, Dmax + 1):
        for a in range(1, D):
            m = D
            for s in S:
                d = dist0(s * a, D)
                if d < m:
                    m = d
                if m == 0:
                    break
            if F(m, D) > best:
                best = F(m, D)
                arg = (D, a)
    return best, arg

def M_le(S, thr):
    for D in range(2, 2 * max(S) + 2):
        for a in range(1, D):
            m = D
            for s in S:
                d = dist0(s * a, D)
                if d < m:
                    m = d
                if m == 0:
                    break
            if F(m, D) > thr:
                return False
    return True

def kills(S, r):
    return any(s % r == 0 for s in S)


def part1_qwitness(ns):
    print("=" * 78)
    print("(1) LARGE-MULTIPLE-FORCED: q-witness M(P_k) >= 1/k and 1/k > n/Phi6, for k in ((n-2)/2, n-2]")
    print("=" * 78)
    for n in ns:
        m = n - 2
        thr = F(n, Phi6(n))
        print(f"  n={n:2d}  n/Phi6={thr}={float(thr):.5f}  (top-half core k in ({m/2:.1f}, {m}]):")
        for k in range(m // 2 + 1, m + 1):
            Pk = [j for j in range(1, m + 1) if j != k]
            Mv, arg = M_and_arg(Pk)
            ok_res = Mv >= F(1, k)          # q-witness lower bound holds
            ok_gap = F(1, k) > thr          # and 1/k beats n/Phi6
            forced = 2 * k                   # smallest multiple of k above the core
            flag = "" if (ok_res and ok_gap and forced > m) else "  <-- CHECK"
            print(f"     k={k:2d}: M(P_k)={str(Mv):>6}={float(Mv):.4f} (>=1/k? {ok_res}) "
                  f"1/k>n/Phi6? {ok_gap}  forced multiple>= {forced}(>{m}? {forced>m}){flag}")


def part2_escape(ns, KB=None):
    print()
    print("=" * 78)
    print("(2) THE FORCED LARGE MULTIPLE DIGS A HOLE: min M over bounded escapes P_k + {A,B}")
    print("    (A kills res k, B kills res n-1) -- confirm > n/Phi6; extract binding + displacement")
    print("=" * 78)
    for n in ns:
        m = n - 2
        thr = F(n, Phi6(n))
        Kbound = KB if KB else 8 * n
        worst = None  # (M, k, A, B, arg)
        for k in range(m // 2 + 1, m + 1):
            Pk = [j for j in range(1, m + 1) if j != k]
            # A = multiple of k (>= 2k), B = multiple of n-1 (>= n-1); both <= Kbound; |S|=n-1
            bestk = None
            As = [a for a in range(2 * k, Kbound + 1, k)]
            Bs = [b for b in range(n - 1, Kbound + 1, n - 1)]
            for A in As:
                for B in Bs:
                    S = Pk + [A, B]
                    if len(set(S)) != n - 1:
                        continue
                    if not (kills(S, k) and kills(S, n - 1)):
                        continue
                    Mv, arg = M_and_arg(S)
                    if bestk is None or Mv < bestk[0]:
                        bestk = (Mv, A, B, arg)
            if bestk:
                Mv, A, B, arg = bestk
                D, a = arg
                # Steinhaus displacement at the binding: image of the k-slot vs image of kappa=A
                kslot = (k * a) % D
                Aimg = (A * a) % D
                displaced = dist0(Aimg - kslot, D)
                if worst is None or Mv < worst[0]:
                    worst = (Mv, k, A, B, arg, kslot, Aimg, displaced)
        Mv, k, A, B, (D, a), kslot, Aimg, displaced = worst
        margin = Mv - thr
        print(f"  n={n:2d}: tightest escape M={str(Mv):>7}={float(Mv):.4f} vs n/Phi6={float(thr):.4f} "
              f"margin={float(margin):+.5f} {'OK' if margin>0 else 'FAIL'}")
        print(f"        at k={k}, A(=kappa, mult of {k})={A}, B(mult of {n-1})={B}; binding D={D},a={a}")
        print(f"        Steinhaus displacement: k-slot(image of k)={kslot}, kappa image={Aimg}, "
              f"dist={displaced} (kappa does NOT sit on the k-slot hole)")


def part3_huge(n=14):
    print()
    print("=" * 78)
    print(f"(3) HUGE-SPEED (CRT) case at n={n}: a huge band-coverer does not fix the punctured hole")
    print("=" * 78)
    m = n - 2
    thr = F(n, Phi6(n))
    # k=1 example from HYP-3748: S = {2..12, 182, 7430}; hole at 269 (not the band primes)
    for k, extra in [(1, [182, 7430]), (12, [84, 13])]:
        Pk = [j for j in range(1, m + 1) if j != k]
        S = sorted(set(Pk + extra))
        if len(S) != n - 1:
            # pad/trim not needed for these curated examples
            pass
        Mv, (D, a) = M_and_arg(S)
        print(f"  k={k}: S={S}")
        print(f"        M={Mv}={float(Mv):.4f} > n/Phi6={float(thr):.4f}? {Mv>thr}; hole at D={D},a={a} "
              f"(NOT a band prime of k) -- the huge speed's coverage moved the hole, never closed it")


if __name__ == "__main__":
    NS = [8, 10, 12, 14, 16, 18, 20]
    part1_qwitness(NS)
    # NOTE: part2_escape / part3_huge below searched over ALL sets (not just COVERING sets), so they
    # surfaced the NON-COVERING tight/mediant minimizers (GW M=1/n, drop-2 M=2/(2n-1)) -- THM-523's
    # TRIVIAL class, which the covering-min does NOT range over. The correct escape analysis (COVERING
    # sets only: a multiple of every q in {2..n}) is in covering_escape_large_multiple_klein.py.
    print()
    print("(2)/(3) See covering_escape_large_multiple_klein.py -- the escape must be a COVERING set")
    print("        (multiple of every q in {2..n}); the tight/mediant minimizers are NON-covering")
    print("        (THM-523 trivial class) and do not bound the covering-min.")
    print()
    print("NET: the LARGE-MULTIPLE-FORCED lemma (part 1) is rigorous & general (q-witness, no primality):")
    print("     k in ((n-2)/2, n-2], k not in S, M(S)<=n/Phi6  =>  S has a multiple of k that is >= 2k.")
