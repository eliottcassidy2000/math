#!/usr/bin/env python3
"""
signed_resonance_correction_klein.py  --  klein-2026-07-01-S66

BOUNDING THE SIGNED RESONANCE CORRECTION on the lonely set L_C, via the RIESZ / THETA form.

Setup (S65 / THM-515). Construction C = {1..n-2, n(n-1)} (n=14 -> {1..12, 182}). Lonely set at level r
    L_C(r) = { t in [0,1) : min_{v in C} ||v t|| > r },     1_{L_C}(t) = prod_{v in C} 1_{||v t|| > r}.
A beater must COVER L_C with a new speed w; coverage frac_w = meas(L_C ∩ {||wt||<r'})/meas(L_C). S65 found
frac_w = 2r' on AVERAGE (Weyl) with RESONANCES. This script isolates and BOUNDS the SIGNED correction.

THE SIGNED RESONANCE CORRECTION (Fourier identity, the key object):
  1_{L_C} is REAL and EVEN (t -> 1-t symmetry), so hat1(k) := \hat{1_{L_C}}(k) is real & even, hat1(0)=L=meas.
  The danger band D_w = 1_{||wt||<r'} has \hat{D_w}(k) = 0 unless w|k; \hat{D_w}(jw) = sin(2π j r')/(π j)
  (j!=0), = 2r' (j=0). By Parseval,
      frac_w = 2r'  +  (2/L) * SUM_{j>=1} hat1(j w) * sin(2π j r')/(π j).
  So  correction_w := frac_w - 2r'  =  (2/L) SUM_{j>=1} hat1(jw) sin(2π j r')/(π j)   <-- SIGNED sum.

THE RIESZ / THETA FORM (off-diagonal THM-515). Each factor 1_{||vt||>r} has Fourier coeffs c_0 = 1-2r,
c_a = -sin(2π a r)/(π a) (a!=0). The product is a RIESZ-TYPE product, so
      hat1(k) = SUM_{ (a_v): SUM_v a_v v = k }  prod_v c_{a_v}.
k=0 is THM-515's theta form L = meas(L_C). k=jw is the OFF-DIAGONAL extension. Since c_a decays ~1/a,
hat1(k) is dominated by the SHORTEST representation of k as an integer combination of C. Hence:
  * SMALL k (k=1..12 via the core, or k near a multiple of 182) -> SHORT rep -> hat1(k) LARGE -> RESONANCE.
  * k with only LONG reps -> hat1(k) tiny -> equidistribution (correction ~ 0).
FAR ELEMENT: a huge w has jw (small j) with no short rep -> hat1(jw) tiny -> correction -> 0. THIS is the
Riesz-product PROOF of S65's far-element impotence.

This script: (1) FFT to get hat1(k); (2) show peaks at short-rep frequencies; (3) RECONSTRUCT correction_w
from the Fourier identity and MATCH the direct S65 value; (4) the Riesz BOUND |correction_w| <= tail; show
it -> 0 for far w and is O(1) for small resonant w.
"""
import numpy as np
from fractions import Fraction as F

def norms(v, t):
    x = (v*t) % 1.0
    return np.minimum(x, 1-x)

def c_coeff(a, r):  # Fourier coeff of 1_{||.||>r}: c_0=1-2r, c_a=-sin(2πar)/(πa)
    if a == 0: return 1.0 - 2*r
    return -np.sin(2*np.pi*a*r)/(np.pi*a)

if __name__ == "__main__":
    n = 14; C = list(range(1, n-1)) + [n*(n-1)]; Mc = F(n, n*n-n+1)
    r = 0.07; rp = 0.07                      # set-level r; test radius r' (r=r' clean self-consistent)
    N = 600000; t = np.arange(N)/N
    print(f"n={n}, C={C}, M_C={float(Mc):.5f}; L_C level r={r}, test r'={rp}, 2r'={2*rp}")

    # 1_{L_C}
    G = np.full(N, 1.0)
    for v in C: G = np.minimum(G, norms(v, t))
    f = (G > r).astype(np.float64)
    L = f.mean()
    print(f"L = meas(L_C({r})) = {L:.5f}   (#intervals ~ Cantor set)")

    # FFT -> hat1(k), k=0..N/2. \hat{f}(k) = (1/N) sum_n f_n e(-2πi k n/N)
    Fhat = np.fft.rfft(f) / N            # complex; should be ~real (even f)
    hat1 = Fhat.real                     # real part (imag ~ 0 by symmetry)
    print(f"hat1(0)={hat1[0]:.5f} (=L check); max|imag|={np.abs(Fhat.imag).max():.2e} (should be ~0)")

    # (2) PEAKS of |hat1(k)| -- are they short-rep frequencies (core 1..12, and near multiples of 182)?
    K = 400
    mags = np.abs(hat1[:K+1]).copy(); mags[0] = 0
    top = np.argsort(mags)[::-1][:18]
    print("\n(2) TOP |hat1(k)| for 1<=k<=400 (RESONANT frequencies):")
    def rep_hint(k):
        if 1 <= k <= 12: return f"core speed {k} (rep len 1)"
        if abs(k-182) <= 6: return f"near 182 (killer); 182-k={182-k}"
        # short rep via core: k = sum of a few of 1..12
        if k <= 24: return f"short core rep (<=2 terms)"
        return "?"
    for k in sorted(top):
        print(f"    k={k:>4}: hat1={hat1[k]:+.5f}  |.|={mags[k]:.5f}   {rep_hint(k)}")

    # (2b) TWO-ATOM MODEL: L_C concentrates at t*=n/Phi6 and 1-t*, so hat1(k) ≈ L*cos(2π k t*)*env(k).
    tstar = n/(n*n-n+1.0)
    print(f"\n(2b) TWO-ATOM MODEL hat1(k) ≈ L*cos(2π k t*), t*=n/Phi6={tstar:.5f}, 1/t*={1/tstar:.4f}:")
    for k in [13, 26, 39, 52, 65, 8, 61, 91, 182]:
        model = L*np.cos(2*np.pi*k*tstar)
        print(f"    k={k:>4}: hat1={hat1[k]:+.5f}  L*cos(2πk t*)={model:+.5f}   ratio hat1/L={hat1[k]/L:+.3f} vs cos={np.cos(2*np.pi*k*tstar):+.3f}")

    # (3) RECONSTRUCT correction_w from the Fourier identity; MATCH direct frac_w. Sum ALL harmonics <=N/2.
    Kmax = N//2
    print("\n(3) Fourier identity: correction_w = (2/L) Σ_{j>=1} hat1(jw) sin(2π j r')/(π j)  vs DIRECT (all harmonics)")
    Lmask = f > 0
    def direct_corr(w):
        return (( (norms(w,t) < rp) & Lmask ).mean()/L) - 2*rp
    def fourier_corr(w):
        j = np.arange(1, Kmax//w + 1)          # ALL harmonics jw <= N/2
        return (2.0/L)*np.sum(hat1[j*w]*np.sin(2*np.pi*j*rp)/(np.pi*j))
    print(f"    {'w':>5} {'direct':>10} {'fourier':>10} {'|diff|':>9}   note")
    testw = [13, 26, 39, 6, 7, 12, 61, 182, 183, 91, 5000, 5001]
    for w in testw:
        dc = direct_corr(w); fc = fourier_corr(w)
        note = "RESONANT" if dc > 0.05 else ("anti" if dc < -0.02 else "~equidist")
        print(f"    {w:>5} {dc:>+10.4f} {fc:>+10.4f} {abs(dc-fc):>9.4f}   {note}")

    # (4) THE RIESZ BOUND: |correction_w| <= (2/L) Σ_j |hat1(jw)| * min(2r', 1/(π j)).
    #     Show it -> 0 for FAR w (harmonics land off the short-rep lattice) and is O(1) for small resonant w.
    print("\n(4) RIESZ BOUND |corr_w| <= (2/L)Σ_j |hat1(jw)| min(2r',1/(πj)); tightness + far-element decay")
    def riesz_bound(w):
        j = np.arange(1, Kmax//w + 1)
        w_j = np.minimum(2*rp, 1.0/(np.pi*j))
        return (2.0/L)*np.sum(np.abs(hat1[j*w])*w_j)
    print(f"    {'w':>6} {'|corr|':>9} {'bound':>9}  bound>=|corr|?")
    for w in [13, 6, 61, 182, 200, 500, 1000, 5000]:
        b = riesz_bound(w); c = abs(direct_corr(w))
        print(f"    {w:>6} {c:>9.4f} {b:>9.4f}  {'OK' if b>=c-1e-6 else 'VIOLATION'}")

    # (5) RIGOROUS far-element bound. 1_{L_C} is a finite union of I intervals => TV bound |hat1(k)| <= I/(π k)
    #     (each interval endpoint is a jump; FT of an indicator of I intervals decays pointwise like 1/k).
    #     => |correction_w| <= (2 I / (π L w)) * S(r'), S(r')=Σ_j (1/j) min(2r',1/(πj)) < ∞.  O(1/w) -> 0.
    d = np.diff((f>0).astype(int)); I = int((d==1).sum()) + (1 if f[0]>0 else 0)
    # verify the TV bound holds pointwise
    ks = np.arange(1, 401); tv = I/(np.pi*ks)
    viol = int((np.abs(hat1[1:401]) > tv + 1e-9).sum())
    print(f"\n(5) RIGOROUS bound. L_C = {I} intervals => |hat1(k)| <= I/(πk); pointwise violations 1..400: {viol}")
    j = np.arange(1, 2000); S = np.sum((1.0/j)*np.minimum(2*rp, 1.0/(np.pi*j)))
    const = 2*I/(np.pi*L)*S
    print(f"    S(r')=Σ_j(1/j)min(2r',1/(πj)) = {S:.4f};  |correction_w| <= {const:.1f}/w   (O(1/w) far-element decay)")
    print(f"    {'w':>6} {'|corr|':>9} {'|corr|*w':>10} {'C/w bound':>10}  ok?")
    for w in [200, 500, 1000, 2000, 5000, 10000]:
        c = abs(direct_corr(w)); print(f"    {w:>6} {c:>9.4f} {c*w:>10.2f} {const/w:>10.4f}  {'OK' if const/w>=c-1e-6 else 'VIOL'}")
    print(f"    => |corr_w| decays like 1/w (bounded by {const:.0f}/w): the FAR element covers L_C only to O(1/w).")
    print(f"    This makes S65 far-element impotence a QUANTITATIVE bound (rate 1/w), the Riesz/TV certificate.")
    print("\n=> correction_w is a SIGNED sum of hat1(jw); hat1 lives on short reps of C (Riesz/theta, THM-515);")
    print("   far element's harmonics miss the lattice => correction -> 0 = equidistribution (S65) PROVED-via-Fourier.")
