#!/usr/bin/env python3
r"""
lrc_resonance_fourier_kps_S82.py  (kind-pasteur-2026-07-08-S82, HYP-5357)

ATTEMPT: prove the resonance lemma Var(W) <= c*R2 via the Fourier pair-resonance kernel
(derive the triple-overlap mass via THM-641's method, sum the c_m products).

RESULT: a decisive NEGATIVE structural finding -- the c_m-product Fourier series for Var(W)
does NOT converge by truncation, so neither the Fourier route (truncated) nor the low-order
real-space overlap masses give a clean Var<=c*R2 bound.  The ~96% cancellation is ESSENTIAL
(non-perturbative).  This maps WHY brick (B) is the open analytic mile and redirects to the
exact-moment + decorrelation route (opus near/far, klein LEM), NOT the resonance bound.

THE EXACT FOURIER FORM (derived, correct):
  W(x) = int_0^1 prod_i (1 - chi(frac(y - e_i x))) dy,  chi = 1_{(0,1/7)},
  chi(t) = sum_m c_m e(mt),  c_m = int_0^{1/7} e(-2pi i m t) dt,  |c_m|^2 = sin^2(pi m/7)/(pi^2 m^2)
  (VANISHES at m = 0 mod 7 -- the apex-7 invisibility, THM-637/638).
  => W(x) = sum_omega Wh(omega) e(omega x),
     Wh(omega) = sum_{S, m: sum m_i = 0, sum m_i e_i = -omega} (-1)^|S| prod c_{m_i}.
  Var(W) = sum_{omega != 0} |Wh(omega)|^2   (FFT-verified exact).
  NO |S|=1 terms (balanced m needs sum m=0, impossible with one nonzero entry).
  The 'triple-overlap mass' in Fourier form = the |S|=3 balanced-m c_m products (no
  real-space Bernoulli needed) -- BUT the series is NOT truncatable (below).

THE OBSTRUCTION (measured):
  pair-m=1 diagonal = |c_1|^4 * R2 = 0.280 at the block -- OVERSHOOTS true Var=0.047 by 6x.
  |Wh(d)| measured is ~7x SMALLER than the pair prediction r_d|c_1|^2: the higher-|S| terms
  cancel ~83% AT EACH frequency omega=d.  Order-by-order truncation OSCILLATES (0.28/0.72/
  0.31/1.01/... for (ORD,MM)=(2,1)/(2,2)/(3,2)/(4,2)) -- diverges, does not -> 0.047.
"""
import numpy as np, itertools
from collections import defaultdict

def c(m):
    return 1/7 if m == 0 else (1 - np.exp(-2j*np.pi*m/7)) / (2j*np.pi*m)

E = list(range(11)); k = 11
def W_grid(res=8192):
    Ea = np.array(E); xs = (np.arange(res)+0.5)/res
    ph = np.mod(np.outer(xs, Ea), 1.0); ph.sort(1)
    g = np.diff(ph, 1); wrap = (ph[:, 0]+1-ph[:, -1])[:, None]
    return np.maximum(np.concatenate([g, wrap], 1) - 1/7, 0).sum(1)
def var_trunc(ORD, MM):
    Wh = defaultdict(complex)
    for s in range(2, ORD+1):
        for S in itertools.combinations(range(k), s):
            for m in itertools.product([x for x in range(-MM, MM+1) if x != 0], repeat=s):
                if sum(m) != 0: continue
                omega = -sum(m[t]*E[S[t]] for t in range(s))
                amp = (-1)**s
                for mi in m: amp *= c(mi)
                Wh[omega] += amp
    return sum(abs(v)**2 for om, v in Wh.items() if om != 0)

Wtrue = W_grid().var()
print(f"k=11 block: true Var(W) = {Wtrue:.6f}")
print("  |c_m|^2 = sin^2(pi m/7)/(pi^2 m^2), = 0 at m=0 mod 7 (apex-7 invisibility)")
print(f"  pair-m=1 diagonal |c_1|^4*R2 = {(np.sin(np.pi/7)**2/np.pi**2)**2 * 770:.6f} (overshoots 6x)")
print("  Fourier truncation Var_trunc (ORD=max|S|, MM=max|m|) -- OSCILLATES, does NOT converge:")
for MM in [1, 2, 3]:
    for ORD in [2, 3, 4]:
        vt = var_trunc(ORD, MM)
        print(f"    ORD={ORD} MM={MM}: {vt:.6f}  (/true = {vt/Wtrue:.2f})")
print("  => the ~96% cancellation is ESSENTIAL (non-perturbative); no truncated c_m-product")
print("     bound works.  Resonance lemma Var<=c*R2 is the open mile; use exact-moment +")
print("     decorrelation (opus near/far, klein LEM-005/006), NOT a Fourier resonance bound.")
