#!/usr/bin/env python3
"""
klein-2026-07-08-S192: the a-priori V_i total-variation bound for the k=11
longest-AP=10 tail decorrelation rate (removes the last certification in the
THM-663 density floor).

Setup: family E_d = d*{0,...,9} u {p}, gcd(d,p)=1. Condition on a=frac(dx);
the 10 AP phases are the Steinhaus orbit {frac(j*a): j=0..9}; the outlier phase
u=frac(px) sweeps [0,1) in d equally-spaced conditional values (klein-S191).
Rate: |E[W^i](E_d) - L_i| <= V_i/d,  V_i = E_a[ TV_u W(a;.)^i ].

This script:
 (1) VERIFIES the exact structure: W(a;u) is continuous piecewise-linear in u
     with slopes in {-1,0,+1}, continuous across AP-phase crossings (both sides
     = W_full), so TV_u[W] = sum over AP-gaps G of  min( 2*(len_G - 1/7)_+ , 2/7 ).
 (2) Establishes the ANALYTIC a-priori bounds:
        TV_u[W]   <= 2*W_full(a) <= 2*(6/7) = 12/7,   and <= 7*(2/7)=2  -> <=12/7
        TV_u[W^i] <= i*(6/7)^{i-1} * TV_u[W]  (chain rule, W in [0,6/7])
     hence V_i <= i*(6/7)^{i-1} * (12/7).
 (3) Measures the ACTUAL V_i = E_a[TV_u W^i] to confirm the analytic bound and
     see the true (much smaller) value.
 (4) Computes the D3 sensitivity |dD3/dm_i| at the limit moments (L1,L2,L3),
     forms C = sum_i |dD3/dm_i| * V_i, and checks C <= 3.47 (required for the
     d>=26 tail: |D3-D3_inf|<=C/d, need C/26 <= D3_inf - bar).
 (5) Cross-checks with the direct interval (box) bound and with D2 (PZ).
"""
import numpy as np
from math import gcd
from fractions import Fraction

INV7 = 1.0/7.0
M = 6.0/7.0
bar = 83549.0/252252.0          # honest k=11 (A') bar
print(f"bar_11 = {bar:.6f}   D3_inf(target) ~ 0.4646")

def W_from_phases(phases):
    """uncovered measure W = sum_gaps (gap - 1/7)_+ for phases on the circle."""
    p = np.sort(np.mod(phases, 1.0))
    gaps = np.diff(p)
    gaps = np.append(gaps, 1.0 - p[-1] + p[0])
    return np.maximum(gaps - INV7, 0.0).sum()

def AP_phases(a, K=10):
    return np.mod(a*np.arange(K), 1.0)

# ---------- (1)+(3) measure TV_u[W^i] and compare to analytic bound ----------
def TV_and_Wfull(a, i, Nu=20000):
    """TV_u[W(a;.)^i] on a fine u-grid, plus W_full(a) (no outlier splitting)."""
    ap = AP_phases(a)
    Wfull = W_from_phases(ap)             # outlier sitting on an AP phase
    us = (np.arange(Nu)+0.5)/Nu
    Wv = np.array([W_from_phases(np.append(ap, u)) for u in us])
    Wi = Wv**i
    tv = np.abs(np.diff(Wi)).sum()
    # close the loop u: Nu-1 -> 0 (periodic)
    tv += abs(Wi[0]-Wi[-1])
    return tv, Wfull, Wv

# analytic gap-based TV (exact formula) to VERIFY the fine-grid TV for i=1
def TV_gapformula(a):
    ap = np.sort(AP_phases(a))
    gaps = np.diff(ap); gaps = np.append(gaps, 1.0-ap[-1]+ap[0])
    return np.minimum(2*np.maximum(gaps-INV7,0.0), 2*INV7).sum()

print("\n--- (1) verify TV_u[W] = gap-formula  (a sample) ---")
rng_a = np.linspace(0.003, 0.997, 40)
maxerr = 0.0
for a in rng_a[:8]:
    tv_grid,_,_ = TV_and_Wfull(a,1,Nu=40000)
    tv_form = TV_gapformula(a)
    maxerr = max(maxerr, abs(tv_grid-tv_form))
print(f"max |TV_grid - TV_gapformula| over sample = {maxerr:.5f}  (grid discretization)")

print("\n--- (3) measured V_i = E_a[TV_u W^i]  vs analytic bound i*(6/7)^(i-1)*(12/7) ---")
Na = 400
avec = (np.arange(Na)+0.5)/Na
Vmeas = {}
for i in (1,2,3):
    tvs = []
    for a in avec:
        tv,_,_ = TV_and_Wfull(a,i,Nu=6000)
        tvs.append(tv)
    Vmeas[i] = np.mean(tvs)
    Vbound = i*(M**(i-1))*(12.0/7.0)
    print(f"  i={i}:  V_i(measured)={Vmeas[i]:.4f}   analytic bound={Vbound:.4f}   "
          f"max_a TV={max(tvs):.4f}")

Vbound = {i: i*(M**(i-1))*(12.0/7.0) for i in (1,2,3)}

# ---------- (4) limit moments L1,L2,L3 and D3 sensitivity ----------
print("\n--- (4) limit moments (block_10 (+) independent outlier) ---")
# L_i = E_{a,u iid}[W(a;u)^i]
Na2, Nu2 = 700, 700
av = (np.arange(Na2)+0.5)/Na2
uv = (np.arange(Nu2)+0.5)/Nu2
S1=S2=S3=0.0
for a in av:
    ap = AP_phases(a)
    Wv = np.array([W_from_phases(np.append(ap,u)) for u in uv])
    S1+=Wv.sum(); S2+=(Wv**2).sum(); S3+=(Wv**3).sum()
tot = Na2*Nu2
L1,L2,L3 = S1/tot, S2/tot, S3/tot
print(f"  L1={L1:.6f}  L2={L2:.6f}  L3={L3:.6f}")

def D3_of(m1,m2,m3):
    N = m1 - m2/M
    D = m2 - m3/M
    return m1/M + N*N/D
def D2_of(m1,m2):
    return m1*m1/m2

D3inf = D3_of(L1,L2,L3)
D2inf = D2_of(L1,L2)
print(f"  D3_inf = {D3inf:.6f}   D2_inf = {D2inf:.6f}   (bar={bar:.4f})")

# partials of D3 at the limit
N = L1 - L2/M; D = L2 - L3/M
dD3_dm1 = 1.0/M + 2*N/D
dD3_dm2 = -2*N/(M*D) - N*N/(D*D)
dD3_dm3 = N*N/(M*D*D)
print(f"  N={N:.5f}  Delta={D:.5f}")
print(f"  |dD3/dm1|={abs(dD3_dm1):.4f}  |dD3/dm2|={abs(dD3_dm2):.4f}  |dD3/dm3|={abs(dD3_dm3):.4f}")

# C from analytic V bound and from measured V
C_bound = abs(dD3_dm1)*Vbound[1] + abs(dD3_dm2)*Vbound[2] + abs(dD3_dm3)*Vbound[3]
C_meas  = abs(dD3_dm1)*Vmeas[1] + abs(dD3_dm2)*Vmeas[2] + abs(dD3_dm3)*Vmeas[3]
req = (D3inf - bar)*26
print(f"\n  C (analytic V bound, triangle-ineq) = {C_bound:.3f}")
print(f"  C (measured  V,       triangle-ineq) = {C_meas:.3f}")
print(f"  required C for d>=26 closure         = (D3inf-bar)*26 = {req:.3f}")
print(f"  ==> triangle-ineq bound closes? analytic:{C_bound<=req}  measured:{C_meas<=req}")

# ---------- (5) direct interval (box) bound over d>=26 ----------
print("\n--- (5) direct box bound: m_i in [L_i - V_i/d, L_i + V_i/d], d>=26 ---")
def D3_box_min(L, Vb, d, ngrid=25):
    """min D3 over the box using analytic V bound (rigorous outer box)."""
    lo = [L[i]-Vb[i+1]/d for i in range(3)]
    hi = [L[i]+Vb[i+1]/d for i in range(3)]
    # clamp to physical: W in [0,6/7] => m1>=0, m2>=m1^2? just clamp >0
    best = 1e9
    g = np.linspace(0,1,ngrid)
    for x in g:
        m1 = lo[0]+(hi[0]-lo[0])*x
        for y in g:
            m2 = lo[1]+(hi[1]-lo[1])*y
            for z in g:
                m3 = lo[2]+(hi[2]-lo[2])*z
                Del = m2 - m3/M
                if Del <= 1e-9: continue
                v = D3_of(m1,m2,m3)
                if v < best: best = v
    return best
for d in (26, 30, 40, 60, 100):
    b_an  = D3_box_min([L1,L2,L3], Vbound, d)
    print(f"  d={d:3d}: box-min D3 (analytic V) = {b_an:.4f}   (bar={bar:.4f})  clears:{b_an>=bar}")

# ---------- (5b) D2 (Paley-Zygmund) a-priori for the tail ----------
print("\n--- (5b) D2 (PZ) tail: only needs V_1,V_2 ---")
dD2_dm1 = 2*L1/L2; dD2_dm2 = -L1*L1/(L2*L2)
C2_bound = abs(dD2_dm1)*Vbound[1] + abs(dD2_dm2)*Vbound[2]
req2 = (D2inf - bar)*26
print(f"  |dD2/dm1|={abs(dD2_dm1):.4f} |dD2/dm2|={abs(dD2_dm2):.4f}")
print(f"  C2 (analytic V) = {C2_bound:.3f}   required = (D2inf-bar)*26 = {req2:.3f}   "
      f"clears:{C2_bound<=req2 and D2inf>=bar}")
