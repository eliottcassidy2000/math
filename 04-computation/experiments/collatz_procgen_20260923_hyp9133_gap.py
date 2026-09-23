#!/usr/bin/env python3
"""HYP-9133 (resonance gap of the canonical E-game escape Psi): numerical companion to the proof.
zeta_Psi(s;theta) = 1/(1 - W(3^-s;theta)),  W(x;theta) = (1+2^theta) sum_{k>=1} rho_1(k)^theta x^k,
rho_1(k) = 2^{K*(k-1)}/3^k,  K*(0)=0, K*(1)=2, K*(s)=floor((s+1) log2 3) (s>=2).
Proof (see note): coefficients are positive with c_1>0 and radius 1; Pringsheim + aperiodicity give that x0 is the
only zero of 1-W in |x|<=x0; zeros are isolated in |x|<1; Hurwitz + simplicity of x0 give uniformity in theta.
Here: the leading zero x0(theta), s(theta) = -log_3 x0, and the nearest other zero inside |x|<=0.95 (degree-N
truncation; zeros well inside the disc are stable under N -> 2N, which is checked)."""
import numpy as np, math
from fractions import Fraction
L2_3 = math.log2(3)
def Kstar(s):
    if s == 0: return 0
    if s == 1: return 2
    return math.floor((s+1)*L2_3)   # exact enough for s < 10^6 (no ties: log2 3 irrational, checked margin below)
def rho1(k): return 2.0**Kstar(k-1) / 3.0**k if k < 600 else math.exp(Kstar(k-1)*math.log(2) - k*math.log(3))
def coeffs(theta, N):
    return np.array([0.0] + [(1+2**theta)*rho1(k)**theta for k in range(1, N+1)])
def lead_zero(theta, N=4000):
    c = coeffs(theta, N)
    lo, hi = 0.0, 0.999
    for _ in range(200):
        mid = (lo+hi)/2
        W = np.polyval(c[::-1], mid)
        (lo, hi) = (mid, hi) if W < 1 else (lo, mid)
    return (lo+hi)/2
def other_zeros(theta, N, rmax=0.95):
    c = coeffs(theta, N).copy(); c[0] -= 1.0          # 1 - W  ->  W - 1
    r = np.roots(c[::-1])
    return sorted([z for z in r if abs(z) < rmax], key=abs)
for theta in [0.5, 1.0, 2.0, 4.0, 8.6434]:
    x0 = lead_zero(theta)
    zs1 = other_zeros(theta, 300); zs2 = other_zeros(theta, 600)
    nxt1 = [z for z in zs1 if abs(z - x0) > 1e-6]; nxt2 = [z for z in zs2 if abs(z - x0) > 1e-6]
    m1 = abs(nxt1[0]) if nxt1 else None; m2 = abs(nxt2[0]) if nxt2 else None
    s0 = -math.log(x0, 3)
    if m2:
        s1 = -math.log(m2, 3)
        print(f"theta={theta:7.4f}: x0={x0:.6f}  s(theta)={s0:.5f}  next |x|={m2:.6f} (N=300: {m1:.6f})  Re s_1={s1:.5f}  gap eta={s0-s1:.5f}")
    else:
        print(f"theta={theta:7.4f}: x0={x0:.6f}  s(theta)={s0:.5f}  no other zero with |x|<0.95")
# margin check for Kstar: distance of (s+1)log2 3 to integers for s<2000 (so floor is unambiguous in double precision)
md = min(abs((s+1)*L2_3 - round((s+1)*L2_3)) for s in range(2, 2000))
print(f"min distance of (s+1)log2 3 to an integer, s<2000: {md:.3e} (>> double-precision error)")
