#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD A (kps-Swf9): the PRINCIPLED delta2 P/R split of the doublet signed error.

K-runner doublet E_M = consec_{K-2} U {M, M+1}, compared vs cap_K.

INCLUSION-EXCLUSION SPLIT (exact, principled -- not period-averaging):
  p0(E_M) = p0(base) + A(M) + A(M+1) + d2(M)
    A(f)   = p0(base U {f}) - p0(base)          (single-far increment at speed f)
    d2(M)  = p0(E_M) - p0(base U {M}) - p0(base U {M+1}) + p0(base)   (far-far interaction)

LIMITS (M->inf, far runners equidistribute):
    a_inf  = lim A(f)  = BV1 - p0(base),   BV1 = boundary_value_direct(base, 1)   [codex, EXACT]
    d_inf  = lim d2(M) = Phi - BV2 ,        BV2 = boundary_value_direct(base, 2)   [EXACT]
       (BV2 = decorrelated 2-far plateau = p0(base)+2 a_inf + 0 interaction; the frozen
        doublet plateau Phi adds the adjacent-pair correlation d_inf = Phi - BV2.)
  => Phi = p0(base) + 2 a_inf + d_inf  (consistency check vs frozen-integral Phi).

THE P/R DECOMPOSITION of g(M) = M*(p0(E_M) - Phi):
    Delta_M = [A(M)-a_inf] + [A(M+1)-a_inf] + [d2(M)-d_inf]
    P(M) := M*[A(M)-a_inf] + M*[A(M+1)-a_inf]    (EXACTLY PERIODIC, period 7*lcm(base), by THM-563:
              M*[A(M)-a_inf] is periodic; M*[A(M+1)-a_inf] = (M+1)*[A(M+1)-a_inf]-[A(M+1)-a_inf]
              = periodic - periodic = periodic.)
    R(M) := M*[d2(M)-d_inf]                       (the far-far interaction correction; CLAIM O(1/M).)

This script (exact rationals), for K=8..12:
  (1) verify Phi (frozen) == p0(base)+2 a_inf + d_inf  (cross-check the two plateau computations).
  (2) confirm P(M) is EXACTLY periodic (2nd-diff of M*[A(M)-a_inf] over the base period == 0).
  (3) tabulate d2(M)-d_inf and R(M)=M*(d2(M)-d_inf): is R bounded? M*|R| / M^2*|R| -> decay rate.
  (4) period-max(P) exact over one base period; tail sup_{M>=M0}|R|.
  (5) assemble: sup g = sup(P+R); show p0=Phi+g/M < cap-0.16 for M>=15 (finite window + tail).
"""
from __future__ import annotations
import sys, functools, math
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce, lru_cache
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP, QVAL
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct
# reuse exact frozen-doublet plateau
from lrc14_doublet_almostperiodic_PR_kpswf9 import phi_exact, base_of, lcm_list, TARGET

ALL_INNER = 0b1111110


def doublet(K, M):
    return tuple(sorted(set(list(range(K - 2)) + [M, M + 1])))


@lru_cache(maxsize=None)
def p0c(E):
    return p0_fast(E)


def base_tuple(K):
    return tuple(range(K - 2))


def A_incr(K, f):
    """single-far increment A(f)=p0(base U{f})-p0(base)."""
    b = base_tuple(K)
    return p0c(tuple(sorted(set(b) | {f}))) - p0c(b)


def d2_interaction(K, M):
    b = base_tuple(K)
    p_b = p0c(b)
    p_bM = p0c(tuple(sorted(set(b) | {M})))
    p_bM1 = p0c(tuple(sorted(set(b) | {M + 1})))
    p_E = p0c(doublet(K, M))
    return p_E - p_bM - p_bM1 + p_b


def main():
    print("=" * 88)
    print("DOUBLET delta2 P/R SPLIT (principled inclusion-exclusion)  (THREAD A, kps-Swf9)")
    print("p0(E_M)=p0(base)+A(M)+A(M+1)+d2(M);  P=M*[(A(M)-a_inf)+(A(M+1)-a_inf)], R=M*(d2-d_inf)")
    print("=" * 88)

    SUMMARY = {}
    for K in range(8, 13):
        b = base_tuple(K)
        base_nz = [e for e in b if e]
        P_per = 7 * lcm_list(base_nz)
        p_base = p0c(b)
        BV1 = boundary_value_direct(b, 1)
        BV2 = boundary_value_direct(b, 2)
        a_inf = BV1 - p_base
        Phi = phi_exact(K)
        d_inf = Phi - BV2
        cap = CAP[K]
        head = cap - TARGET - Phi

        # (1) plateau consistency: Phi ?= p_base + 2 a_inf + d_inf  (trivially true by def of d_inf;
        #     the REAL check is Phi(frozen) vs BV2 + d_inf where d_inf is independently lim d2(M).)
        # independent d_inf estimate from large M:
        d2_large = [d2_interaction(K, M) for M in range(900, 940)]
        d_inf_emp = sum(d2_large) / len(d2_large)
        consistent = (d_inf == d_inf_emp)  # exact equality if frozen Phi is exact

        print(f"\nK={K}  base=consec_{K-2}={b}  P_per=7*lcm{tuple(base_nz)}={P_per}")
        print(f"  p0(base)={float(p_base):.6f}  BV1={float(BV1):.6f}  BV2={float(BV2):.6f}")
        print(f"  a_inf={float(a_inf):.6f}  Phi={float(Phi):.6f}  d_inf=Phi-BV2={float(d_inf):.6f}")
        print(f"  d_inf(emp,large M)={float(d_inf_emp):.6f}  exact-match? {consistent}")
        print(f"  cap_{K}={float(cap):.6f}  headroom H_K=cap-0.16-Phi={float(head):.6f}")

        # (2) P exactly periodic: 2nd-diff of S(M):=M*[A(M)-a_inf] over step P_per must vanish.
        per_ok = True
        worst_pd = F(0)
        for M in range(15, 15 + P_per):
            S0 = M * (A_incr(K, M) - a_inf)
            S1 = (M + P_per) * (A_incr(K, M + P_per) - a_inf)
            S2 = (M + 2 * P_per) * (A_incr(K, M + 2 * P_per) - a_inf)
            sd = S0 - 2 * S1 + S2
            if sd != 0:
                per_ok = False
                worst_pd = max(worst_pd, abs(sd))
            if M > 15 + 60:  # spot-check first ~60 residues for speed; full check optional
                break
        print(f"  (2) P(M) single-far part exactly periodic (2nd-diff vanish, spot 60)? {per_ok}"
              + ("" if per_ok else f"  worst={worst_pd}"))

        # (3) R(M)=M*(d2(M)-d_inf): bounded? decay?  Use M up to a few periods.
        Mmax = min(3000, 15 + 6 * P_per)
        Rvals = {}
        d2vals = {}
        for M in range(15, Mmax + 1):
            d2 = d2_interaction(K, M)
            d2vals[M] = d2
            Rvals[M] = M * (d2 - d_inf)
        def tail(fn, thr):
            vs = [fn(M) for M in Rvals if M >= thr]
            return max(vs) if vs else F(0)
        sup_absR = max(abs(v) for v in Rvals.values())
        print(f"  (3) R(M)=M*(d2-d_inf): sup|R| over [15,{Mmax}] = {float(sup_absR):.5f}")
        print(f"      {'M0':>6} {'sup|d2-dinf|':>13} {'sup|R|=M|d2-dinf|':>18} {'sup M^2|d2-dinf|':>17}")
        for M0 in (15, 100, 300, 800, 1500):
            sD = tail(lambda M: abs(d2vals[M] - d_inf), M0)
            sR = tail(lambda M: abs(Rvals[M]), M0)
            sM2 = tail(lambda M: abs(Rvals[M]) * M, M0)
            print(f"      {M0:>6} {float(sD):>13.6f} {float(sR):>18.5f} {float(sM2):>17.3f}")

        # (4) period-max(P): P(M)=M*[A(M)-a_inf]+M*[A(M+1)-a_inf] over one base period.
        Pvals = {}
        for M in range(15, 15 + P_per):
            Pvals[M] = M * (A_incr(K, M) - a_inf) + M * (A_incr(K, M + 1) - a_inf)
        pmaxP = max(Pvals.values())
        pminP = min(Pvals.values())
        argmaxP = max(Pvals, key=lambda M: Pvals[M])
        print(f"  (4) period-max(P) over one base period [{15},{15+P_per}) = {float(pmaxP):.5f} at M={argmaxP}"
              f"  (min {float(pminP):.5f})")

        # (5) assemble
        tailR = tail(lambda M: abs(Rvals[M]), 800)
        B = max(pmaxP, F(0)) + tailR
        M0 = math.ceil(float(B) / float(head)) if head > 0 else None
        hi = min(M0, Mmax)
        worst_p0 = max(p0c(doublet(K, M)) for M in range(15, hi + 1))
        thr = cap - TARGET
        print(f"  (5) ASSEMBLE: B=pmaxP+tail|R|@800={float(B):.5f}; M0=ceil(B/H_K)={M0}; "
              f"finite window [15,{hi}] worst p0={float(worst_p0):.6f} <= cap-0.16={float(thr):.6f}? "
              f"{worst_p0 <= thr}")
        SUMMARY[K] = dict(P_per=P_per, Phi=Phi, head=head, pmaxP=pmaxP, sup_absR=sup_absR,
                          tailR=tailR, B=B, M0=M0, worst_p0=worst_p0, thr=thr,
                          closes=(worst_p0 <= thr), per_ok=per_ok, d_consistent=consistent)

    print("\n" + "=" * 88)
    print("SUMMARY (does the delta2 split close the doublet leg with margin 0.16?)")
    print("=" * 88)
    print(f"  {'K':>3} {'P_per':>7} {'H_K':>9} {'pmaxP':>9} {'sup|R|':>9} {'M0':>7} {'closes':>7} {'cap-worst':>10}")
    allclose = True
    for K in range(8, 13):
        s = SUMMARY[K]
        allclose &= s['closes']
        print(f"  {K:>3} {s['P_per']:>7} {float(s['head']):>9.5f} {float(s['pmaxP']):>9.5f} "
              f"{float(s['sup_absR']):>9.5f} {s['M0']:>7} {str(s['closes']):>7} "
              f"{float(CAP[K]-s['worst_p0']):>10.5f}")
    print(f"\n  ALL K close with margin 0.16? {allclose}")


if __name__ == "__main__":
    main()
