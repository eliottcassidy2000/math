#!/usr/bin/env python3
"""
lrc14_covering_autocorr_xintegral_klein_S287.py
================================================
klein-2026-07-13-S287 (owner: build the middle-order x-integral for the covering side).

THE CONSTRUCTION (covering analogue of THM-729). Good set G' = {t in [0,1): ||w t|| >= 1/14 for all
speeds w in S} (the set of 1/14-lonely times). Bad arc D_w = {||w t|| < 1/14}, measure 1/7. opus's
per-core correlation is eps_v = Cov(1_{D_v}, 1_{G'}) = meas(D_v ^ G') - |D_v||G'|, where |D_v|=1/7.

Because 1_{D_v}(t)=h(v t) with h the measure-1/7 arc indicator, its Fourier coeffs live on the
v-grid: hat(1_{D_v})(l) = hbar(l/v) if v|l else 0. So
    eps_v = sum_{m != 0} hbar(m) chat_{m v},   chat_l = hat(1_{G'})(l).
Cauchy-Schwarz + Wiener-Khinchin + the Poisson v-grid sampling formula give the POSITIVE-DEFINITE bound
    |eps_v|^2 <= (6/49) * [ (1/v) sum_{j=0}^{v-1} A(j/v)  -  |G'|^2 ],           (*)
where A(tau) = meas(G' ^ (G'-tau)) is the good-set AUTOCORRELATION (>=0, positive-definite), and the
bracket = the v-grid DISCREPANCY of A (its grid-mean minus true mean |G'|^2), manifestly >=0. This is
the middle-order x-integral: a positive-definite SPATIAL discrepancy of the good-set self-overlap on
the core-v-grid -- NO Fourier expansion of the product (avoids the S266 divergence).

This script (1) verifies (*) holds for every core v, on the deep well and the near-AP residue; (2)
checks the ACID TEST (mac-mini-S83): does the metric object rank {1..11,13,84} (L=.005) as MORE stuck
than the deep well (L=.024), the ordering that every STRUCTURAL deficit gets WRONG?
"""
import numpy as np

NG = 1 << 21  # 2,097,152 grid points (fine; FFT-friendly). aliasing not delicate for a single set's autocorr.
THR = 1.0/14.0

def good_indicator(S):
    """1_{G'} on the grid: 1 where ||w t||>=1/14 for all w in S."""
    t = np.arange(NG, dtype=np.float64)/NG
    g = np.ones(NG, dtype=np.float64)
    for w in S:
        frac = (w*t) % 1.0
        dist = np.minimum(frac, 1.0-frac)          # ||w t||
        g *= (dist >= THR)
    return g

def analyze(S, name):
    g = good_indicator(S)
    L = g.mean()                                    # |G'|
    # autocorrelation A(tau) = mean_t g(t) g(t+tau), via FFT.  A on the grid, A[k]=A(k/NG).
    G = np.fft.rfft(g)
    A = np.fft.irfft(G*np.conj(G), n=NG)/NG          # circular autocorr, A[0]=|G'|... actually mean of g^2 = L
    # A[k] = (1/NG) sum_t g[t] g[t+k]  -> A[0] = L (since g^2=g). mean_k A[k] = L^2.
    Abar = A.mean()                                  # = L^2
    results=[]
    for v in sorted(set(S)):
        # eps_v exact on grid: meas(D_v ^ G') - (1/7)|G'|
        t = np.arange(NG)/NG
        frac=(v*t)%1.0; dist=np.minimum(frac,1.0-frac)
        Dv = (dist < THR).astype(np.float64)
        eps = (Dv*g).mean() - (1.0/7.0)*L
        # RHS of (*): (6/49)*[ (1/v) sum_{j=0}^{v-1} A(j/v) - L^2 ]
        # A(j/v) via nearest grid index round(j/v * NG)
        idx = np.round((np.arange(v)/v)*NG).astype(np.int64) % NG
        grid_mean = A[idx].mean()                    # (1/v) sum_j A(j/v)
        disc = grid_mean - Abar                      # v-grid discrepancy of A (>=0)
        rhs = (6.0/49.0)*disc
        ok = eps*eps <= rhs + 1e-12
        results.append((v, eps, eps*eps, rhs, disc, ok))
    return L, Abar, results

print("="*94)
print("COVERING MIDDLE-ORDER x-INTEGRAL (THM-729 analogue): |eps_v|^2 <= (6/49)*[v-grid disc of A_G']")
print("A(tau)=meas(G' ^ (G'-tau)) good-set autocorrelation (positive-definite). NG=%d, thr=1/14"%NG)
print("="*94)

FAMILIES = [
    ([1,2,3,4,5,6,7,8,9,10,11,12,182], "deep well {1..12,182}   L~.024"),
    ([1,2,3,4,5,6,7,8,9,10,11,13,84],  "near-AP residue {1..11,13,84} L~.005"),
]
summary=[]
for S,name in FAMILIES:
    L,Abar,res = analyze(S,name)
    print("\n%s   |G'|=%.5f   mean-autocorr L^2=%.6f"%(name,L,Abar))
    print("   %4s %12s %12s %12s %10s %4s"%("v","eps_v","eps_v^2","(6/49)disc","disc_v","(*)"))
    allok=True; maxdisc=0.0; sumeps2=0.0
    for v,eps,eps2,rhs,disc,ok in res:
        allok &= ok; maxdisc=max(maxdisc,disc); sumeps2+=eps2
        flag = "OK" if ok else "**FAIL**"
        # print the tightest / a few
        if v in (min(s for s in [r[0] for r in res]),) or v>=13 or not ok or v==max(r[0] for r in res):
            print("   %4d %12.6f %12.3e %12.3e %10.3e %4s"%(v,eps,eps2,rhs,disc,flag))
    print("   -> bound (*) holds for ALL cores: %s ;  max disc_v=%.3e ; sum eps_v^2=%.3e"%(allok,maxdisc,sumeps2))
    summary.append((name,L,maxdisc,sumeps2))

print("\n"+"="*94)
print("ACID TEST (mac-mini-S83): a METRIC object must rank the near-AP residue as MORE stuck (smaller L).")
print("   family                              |G'|      max_v disc_v   sum_v eps_v^2")
for name,L,md,se in summary:
    print("   %-36s %.5f   %.4e   %.4e"%(name,L,md,se))
print("\n   Structural deficits get this WRONG (residue dominates deep well on Schur/Q2 yet has smaller L).")
print("   If disc_v / sum eps_v^2 is LARGER for the residue, the autocorrelation discrepancy is a")
print("   faithful METRIC surrogate (tracks L, not structure). If not, it inherits the same decoupling.")
print("done.")
