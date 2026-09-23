#!/usr/bin/env python3
"""Orchestrator's independent audit of the gamma = 7/20 two-disk certificate (Theorem A, C* >= 27/20).
(1) containment: |p| S(p)^g < 1 on the circle |p - c| = rho (then the whole disk, by star-shapedness of
    Omega_0 w.r.t. 0 and convexity of the disk containing 0); dense sampling + explicit Lipschitz margin;
(2) the closed-form conformal radius log R(w(V),0) at (c,rho) and its sanity value at (0,1) = log(8 sqrt3/9).
Independent code (numpy + mpmath); no reuse of the lane's modules."""
import numpy as np, mpmath as mp
g = 7/20; c = 13/200; rho = 39/50
N = 4_000_000
t = np.linspace(0, 2*np.pi, N, endpoint=False)
p = c + rho*np.exp(1j*t)
S = np.abs(p) + np.abs(1-p)
f = np.abs(p) * S**g
fmax = f.max()
# Lipschitz bound of f along the circle w.r.t. arc length: |d/ds f| <= |grad| ; |grad|p|| = 1, |grad S| <= 2,
# so |grad f| <= S^g + |p| g S^(g-1) * 2 ; bound by max over samples plus safety factor 1.5
lip = (S**g + 2*g*np.abs(p)*S**(g-1)).max()*1.5
spacing = 2*np.pi*rho/N
bound = fmax + lip*spacing/2
print(f"(1) containment gamma=7/20: max f on {N} samples = {fmax:.8f}; Lipschitz {lip:.4f}; spacing {spacing:.2e}; certified max <= {bound:.8f} -> {'PASS' if bound < 1 else 'FAIL'}")
mp.mp.dps = 40
def logR(c, rho):
    c = mp.mpf(c); rho = mp.mpf(rho)
    s = mp.sqrt(rho**2 - (mp.mpf(1)/2 - c)**2)
    b0 = mp.atan(s/(mp.mpf(1)/2 + rho - c))
    kap = mp.pi/(2*mp.pi - 4*b0)
    A = 2*kap*(mp.atan(2*s) - b0)
    return mp.log(mp.tan(A)*(mp.mpf(1)/4 + s**2)/(kap*s))
v0 = logR(0, 1); ref = mp.log(8*mp.sqrt(3)/9)
print(f"(2) sanity: logR(0,1) = {mp.nstr(v0,20)} vs log(8 sqrt3/9) = {mp.nstr(ref,20)}  diff {mp.nstr(v0-ref,5)}")
v = logR(mp.mpf(13)/200, mp.mpf(39)/50)
print(f"    logR(13/200, 39/50) = {mp.nstr(v,20)}  -> {'POSITIVE' if v > 0 else 'NOT POSITIVE'} (lane: 0.0059354)")
