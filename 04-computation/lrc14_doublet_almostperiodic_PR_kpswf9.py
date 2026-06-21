#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD A (kps-2026-06-21-Swf9): the ALMOST-PERIODIC P/R DECOMPOSITION of the binding
genuine-wide DOUBLET signed error.

K-RUNNER convention (CORRECTED indexing): the binding genuine-wide maximizer for K runners is
    E_M = consec_{K-2} U {M, M+1}     ( base = {0,1,...,K-3} (K-2 elts) + far doublet {M,M+1} )
=> |E| = K runners, compared against cap_K.  (The opus signed-bound script used
consec_{k-1}+{M,M+1} which is the (k+1)-runner family mislabeled k; that off-by-one made the
crude test 'fail at k=9'. Here everything is K-runner-correct.)

THM-563 closes SINGLE-far by exact periodicity of w*Delta_w. The DOUBLET is NOT exactly periodic
(2nd-diff ~0.01-0.03 nonzero, decaying) -> ALMOST-periodic. Decompose
    g(M) := M * (p0(E_M) - Phi) = P(M) + R(M)
  - Phi  = EXACT plateau lim_{M->inf} p0(E_M), computed via the FROZEN-DOUBLET integral.
  - P(M) = the exactly-periodic part (period P_per = 7*lcm(base)).
  - R(M) = g(M) - P(M): the w-dependent breakpoint correction (claimed O(1/M) decaying).

ASSEMBLY (margin-to-cap target 0.16):  p0(E_M) = Phi + g(M)/M.
  headroom H_K = cap_K - 0.16 - Phi   (>0 needed). For M>M0:
     p0(E_M) - (cap_K-0.16) = g(M)/M - H_K <= [period-max(P)+tail|R|]/M - H_K < 0
  once M > [period-max(P)+tail|R|]/H_K =: M0.  Window [15,M0] checked EXACTLY.
"""
from __future__ import annotations
import sys, functools, math
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd, comb
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP, QVAL, MARGIN

ALL_INNER = 0b1111110
TARGET = F(16, 100)  # required margin cap_K - sup p0 >= 0.16


def lcm_list(xs):
    return reduce(lambda a, b: a * b // gcd(a, b), xs, 1)


def doublet(K, M):
    """K-runner binding doublet: base consec_{K-2}={0..K-3}, far={M,M+1}. |E|=K."""
    return tuple(sorted(set(list(range(K - 2)) + [M, M + 1])))


def base_of(K):
    return list(range(K - 2))


# ===================================================================================
# EXACT PLATEAU Phi = lim_{M->inf} p0(E_M)  via the FROZEN-DOUBLET model.
#
# p0(E) = meas{x in[0,1): {floor(7 e x) mod 7 : e in E} covers inner sectors 1..6}.
# As M->inf with far={M,M+1}: at base coordinate y=x the two far runners have phases
# theta=M x mod1 and (M+1)x = theta + x mod1. theta equidistributes uniformly & indep of
# the base cell; the OFFSET between the two far runners is exactly x=y (slow). So the frozen
# far law at base-coord y: theta~Unif[0,1), far sectors a=floor(7 theta), b=floor(7(theta+y)).
#   Phi = INT_{y in[0,1)} Pr_theta[ base_sec(y) U {a,b} >= inner 1..6 ] dy   (EXACT rational).
# ===================================================================================

def base_sectors_at_mid(base_nz, ymid):
    mask = 0
    for e in base_nz:
        sec = int(e * ymid * 7) % 7  # floor(7 e y) mod 7 (y>=0)
        mask |= 1 << sec
    return mask


def frozen_far_cover_prob(missing_mask, y):
    """Pr_{theta~U[0,1)}[ {floor(7 theta), floor(7(theta+y))} covers missing_mask ]. Exact."""
    miss = [s for s in range(7) if (missing_mask >> s) & 1]
    if len(miss) > 2:
        return F(0)
    miss_set = set(miss)
    fy = 7 * y
    total = F(0)
    for a in range(7):
        lo_u, hi_u = F(a), F(a + 1)
        nlo = math.floor(a + fy) + 1
        nhi = math.ceil(a + 1 + fy) - 1
        bps = [lo_u]
        for n in range(nlo, nhi + 1):
            ub = F(n) - fy
            if lo_u < ub < hi_u:
                bps.append(ub)
        bps.append(hi_u)
        bps = sorted(set(bps))
        for ulo, uhi in zip(bps, bps[1:]):
            if uhi <= ulo:
                continue
            umid = (ulo + uhi) / 2
            b = int(umid + fy) % 7
            if miss_set <= {a % 7, b}:
                total += (uhi - ulo) / 7
    return total


def phi_exact(K):
    base_nz = [e for e in base_of(K) if e]  # 1..K-3
    l = lcm_list(base_nz)
    bps = {F(0), F(1)}
    for e in base_nz:
        for j in range(7 * e + 1):
            bps.add(F(j, 7 * e))
    for n in range(8):
        bps.add(F(n, 7))
    bps = sorted(b for b in bps if F(0) <= b <= F(1))
    total = F(0)
    for ylo, yhi in zip(bps, bps[1:]):
        if yhi <= ylo:
            continue
        ymid = (ylo + yhi) / 2
        bmask = base_sectors_at_mid(base_nz, ymid) & ALL_INNER
        if bmask.bit_count() == 6:
            total += (yhi - ylo)
            continue
        missing = (~bmask) & ALL_INNER
        total += (yhi - ylo) * frozen_far_cover_prob(missing, ymid)
    return total


# ===================================================================================
# P/R SPLIT.  g(M) = M*(p0(E_M)-Phi). Build P(M) as the exactly-periodic component by
# PERIOD-AVERAGING over the base period P_per=7*lcm(base): for residue r mod P_per,
#   P(r) := mean over a finite set of M ≡ r (mod P_per), large M, of g(M).
# Then R(M)=g(M)-P(M mod P_per). If g is almost-periodic with period P_per and decaying
# correction, P captures the periodic part and R->0. We measure the decay of R.
# (We use several large blocks of M to define P stably; see verify in main.)
# ===================================================================================

def g_of(K, M, Phi):
    return M * (p0_fast(doublet(K, M)) - Phi)


def step2_split_and_decay(K, Phi, Mmax=3000, verbose=True):
    """Build P(M) by period-averaging g over P_per=7*lcm(base); R=g-P; measure decay.
    Returns dict with period-max(P), R-tail bounds, and the global sup of g(M)."""
    base_nz = [e for e in base_of(K) if e]
    P_per = 7 * lcm_list(base_nz)
    # compute g over [15, Mmax]
    g = {}
    for M in range(15, Mmax + 1):
        g[M] = g_of(K, M, Phi)
    # P(r): period-average using the LATE blocks (largest M available) for each residue.
    # Use the last full period present in [15,Mmax] plus the one before to stabilize.
    # P(r) := average of g over the available M ≡ r (mod P_per) within a late window.
    late_lo = max(15, Mmax - 6 * P_per + 1)
    by_res = {}
    for M in range(late_lo, Mmax + 1):
        by_res.setdefault(M % P_per, []).append(g[M])
    Pfun = {r: (sum(v) / len(v)) for r, v in by_res.items()}
    # period-max of P over residues we have
    pmaxP = max(Pfun.values()) if Pfun else None
    pminP = min(Pfun.values()) if Pfun else None
    # R(M) = g(M) - P(M mod P_per)  (only where residue present)
    R = {M: g[M] - Pfun[M % P_per] for M in g if (M % P_per) in Pfun}
    # decay tables: bucket by M-decade, report max |R|, max M*|R|, max M^2*|R|
    sup_g = max(g.values())
    arg_sup_g = max(g, key=lambda M: g[M])
    # tail sup |R| for M >= threshold
    def tail_supR(thr):
        vs = [abs(R[M]) for M in R if M >= thr]
        return max(vs) if vs else F(0)
    def tail_supMR(thr):
        vs = [abs(R[M]) * M for M in R if M >= thr]
        return max(vs) if vs else F(0)
    def tail_supM2R(thr):
        vs = [abs(R[M]) * M * M for M in R if M >= thr]
        return max(vs) if vs else F(0)
    return dict(K=K, P_per=P_per, pmaxP=pmaxP, pminP=pminP, sup_g=sup_g, arg_sup_g=arg_sup_g,
                R=R, g=g, Pfun=Pfun, tail_supR=tail_supR, tail_supMR=tail_supMR,
                tail_supM2R=tail_supM2R, Mmax=Mmax)


def main():
    print("=" * 86)
    print("DOUBLET ALMOST-PERIODIC P/R DECOMPOSITION  (THREAD A, kps-Swf9, K-runner-correct)")
    print("E_M = consec_{K-2} U {M,M+1}    g(M)=M*(p0(E_M)-Phi)=P(M)+R(M)   target margin 0.16")
    print("=" * 86)

    # STEP 1: exact Phi, margins, headroom
    print("\n[STEP 1] EXACT plateau Phi (frozen-doublet integral); margins & headroom")
    print("-" * 86)
    PHI = {}
    for K in range(8, 13):
        ph = phi_exact(K)
        PHI[K] = ph
        cap = CAP[K]
        head = cap - TARGET - ph
        print(f"  K={K}: Phi={ph}={float(ph):.6f}  cap={float(cap):.6f}  "
              f"cap-Phi={float(cap-ph):.6f}")
        print(f"        headroom H_K = cap-0.16-Phi = {float(head):.6f}  "
              f"({'OK >0' if head > 0 else 'NEG'})  15*H_K={float(15*head):.4f}")
    print("\n  (Need sup_{M>=15} g(M) <= 15*H_K via finite window + tail; tightest K has smallest H_K.)")

    # STEP 2/3: P/R split + decay, and global sup of g(M). Mmax per K (period-aware).
    print("\n[STEP 2/3] P/R split (period-avg) + decay of R; global sup g(M)=M*Delta_M")
    print("-" * 86)
    MMAX = {8: 4200, 9: 6000, 10: 6000, 11: 5400, 12: 5400}
    res_by_K = {}
    for K in range(8, 13):
        Phi = PHI[K]
        Mmax = MMAX[K]
        r = step2_split_and_decay(K, Phi, Mmax=Mmax)
        res_by_K[K] = r
        head = CAP[K] - TARGET - Phi
        print(f"\n  K={K}  P_per=7*lcm(base)={r['P_per']}  Mmax={Mmax}")
        print(f"    global sup_M g(M) over [15,{Mmax}] = {float(r['sup_g']):.5f} at M={r['arg_sup_g']}  "
              f"(15*H_K={float(15*head):.4f})  sup_g<15*H_K? {r['sup_g'] < 15*head}")
        print(f"    period-max P(r): max={float(r['pmaxP']):.5f}  min={float(r['pminP']):.5f}")
        print(f"    R decay (tail from threshold M0):")
        print(f"      {'M0':>6} {'sup|R|':>10} {'M0*sup|R|':>12} {'sup M|R|(tail)':>15} {'sup M^2|R|':>14}")
        for M0 in (15, 100, 300, 1000, 2000):
            sR = r['tail_supR'](M0)
            sMR = r['tail_supMR'](M0)
            sM2R = r['tail_supM2R'](M0)
            print(f"      {M0:>6} {float(sR):>10.5f} {float(M0*sR):>12.4f} {float(sMR):>15.4f} {float(sM2R):>14.2f}")

    # STEP 4: ASSEMBLY -- finite window M0 + tail bound -> p0 < cap-0.16
    print("\n[STEP 4] ASSEMBLY: p0(E_M)=Phi+g(M)/M < cap-0.16 for all M>=15")
    print("-" * 86)
    print(f"  {'K':>3} {'H_K':>9} {'pmaxP':>9} {'tail|R|@1000':>13} {'B=pmax+tail':>12} {'M0=ceil(B/H_K)':>14} {'window<cap-.16':>15}")
    for K in range(8, 13):
        r = res_by_K[K]
        Phi = PHI[K]
        head = CAP[K] - TARGET - Phi
        # B = bound on |g(M)| for M >= M0: use period-max(P) (an upper env of P) + tail sup|R|
        tail = r['tail_supR'](1000)
        B = max(r['pmaxP'], F(0)) + tail
        M0 = math.ceil(float(B) / float(head)) if head > 0 else None
        # finite window [15, M0]: exact check p0 <= cap-0.16
        cap = CAP[K]
        thr = cap - TARGET
        hi = min(M0, r['Mmax'])
        worst_p0 = max(p0_fast(doublet(K, M)) for M in range(15, hi + 1))
        win_ok = worst_p0 <= thr
        print(f"  {K:>3} {float(head):>9.5f} {float(r['pmaxP']):>9.5f} {float(tail):>13.5f} "
              f"{float(B):>12.5f} {M0:>14} {str(win_ok):>15}")
        print(f"        finite window [15,{hi}]: worst p0={float(worst_p0):.6f} <= cap-0.16={float(thr):.6f}? {win_ok}  "
              f"(at boundary cap-worst={float(cap-worst_p0):.5f})")
    print("\n" + "=" * 86)
    print("CLOSURE: if (a) sup_g<15*H_K (tail bounded), (b) R decays (M|R| bounded -> R=O(1/M)),")
    print("(c) finite window [15,M0] has p0<=cap-0.16, then the doublet leg CLOSES with margin 0.16.")
    return PHI, res_by_K


if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()
