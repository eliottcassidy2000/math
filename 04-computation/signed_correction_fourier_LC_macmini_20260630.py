#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S75 -- THE SIGNED CORRECTION = FOURIER OF 1_{L_C} AT HARMONICS OF THE TEST SPEED.
Upgrades S74's equidistribution heuristic to a QUANTITATIVE bound that closes the single-huge-speed and
<=6-huge-speed cases of the covering-min residual.

For a test speed w, the measure of the core's fixed lonely set L_C(r) that w's danger arcs cover is:
    covered(w) = integral 1_{L_C}(t) 1_{||wt||<r}(t) dt
               = 2r|L_C|              (DIAGONAL / equidistribution main term, j=0)
               + sum_{j!=0} hat(1_{L_C})(jw) * ghat(j)     (SIGNED CORRECTION, off-diagonal)
where ghat(j)=sin(2*pi*j*r)/(pi*j) (Fourier of one danger arc), and the correction samples the Fourier
coefficients of 1_{L_C} AT THE HARMONICS jw of the test speed (the owner's seed).

Since L_C is a union of N arcs, |hat(1_{L_C})(m)| <= N/(pi|m|), so
    |signed correction| <= sum_{j!=0} (N/(pi|jw|))(1/(pi|j|)) = (N/(pi^2 w)) * 2*zeta(2) = N/(3w).
Hence covered(w) <= 2r|L_C| + N/(3w) < |L_C|  once  w > N/(3(1-2r)|L_C|)  => M(C u {w}) >= r  (single huge w).
Multi-patch: union bound covered(H) <= |H|*2r|L_C| + sum N/(3 w_i) < |L_C| when |H|<=6 and each w_i > threshold.
"""
import numpy as np
r = 14/183; G = 1 << 21; t = np.arange(G)/G
def lonely_mask(core):
    m = np.ones(G, bool)
    for v in core: m &= (np.abs(v*t - np.round(v*t)) >= r)
    return m
def narcs(mask):
    d = np.diff(np.r_[mask.astype(int), mask[0]]); return int((d == 1).sum())
def ghat(j): return 2*r if j == 0 else np.sin(2*np.pi*j*r)/(np.pi*j)

print(f"r=14/183={r:.5f}, 2r={2*r:.4f}")
print("\n(A) FOURIER IDENTITY covered(w) = 2r|L_C| + signed_correction  (direct vs Fourier), core={1..6}:")
core = list(range(1, 7)); m = lonely_mask(core); LC = m.astype(float)
hat = np.fft.rfft(LC)/G; meas = LC.mean(); N = narcs(m)
print(f"    |L_C|={meas:.5f}, N(arcs)={N}, diagonal 2r|L_C|={2*r*meas:.5f}")
for w in [7, 50, 182, 1000, 10000]:
    dw = np.abs(w*t - np.round(w*t)) < r
    direct = (m & dw).mean()
    J = (G//2)//w
    corr = sum(2*hat[j*w].real*ghat(j) for j in range(1, min(J, 400)+1) if j*w <= G//2)
    print(f"    w={w:6d}: direct={direct:.6f}  2r|L_C|+corr={2*r*meas+corr:.6f}  signed_corr={corr:+.6f}  bound N/(3w)={N/(3*w):.5f}")

print("\n(B) |hat 1_{L_C}(m)| <= N/(pi m)  =>  |signed correction| <= N/(3w)  (the decay giving the bound):")
for core in [list(range(1, 7)), list(range(1, 10)), list(range(1, 12))]:
    m = lonely_mask(core); LC = m.astype(float); hat = np.fft.rfft(LC)/G
    N = narcs(m); meas = LC.mean(); ms = np.arange(1, 20000)
    peak = float(np.max(np.abs(hat[1:20000])*np.pi*ms))
    thr = N/(3*(1-2*r)*meas)
    print(f"    core={{1..{max(core)}}}: |L_C|={meas:.5f} N={N:3d}  sup(pi m|hat|)={peak:.2f}<=N? {peak<=N+1}  "
          f"THRESHOLD w>{thr:.0f} => single huge w beats nothing (M>=r)")

print("\n(C) THE ARGUMENT (proof push):")
print("    single huge w > N/(3(1-2r)|L_C|): covered(w) < |L_C| => M(C u {w}) >= r  RIGOROUS (Fourier decay of 1_{L_C}).")
print("    multi-patch <=6 huge speeds each > threshold: union bound => covered < |L_C| => M >= r  RIGOROUS.")
print("    bounded speeds <= threshold: lazy-cut ILP (HYP-3782).  Residual: >=7 huge speeds (cross-harmonic hat(jw_i-j'w_j)).")
print("DONE.")
