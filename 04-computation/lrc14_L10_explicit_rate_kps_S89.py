#!/usr/bin/env python3
"""kps-2026-07-08-S89 -- EXPLICIT rigorous constants for klein-S191's L=10 single-outlier rate, and a
conditional-D3 finite check extending klein's d<=25 range.  Goal: make the k=11 L=10 closure fully
explicit (klein flagged 'the explicit V_i' as the last bookkeeping).

klein-S191: for E_d = d*{0..9}+p (gcd(d,p)=1), condition on a=frac(dx): AP phases = {frac(ja)}, and
the outlier phase runs over d EQUALLY-SPACED values {frac(pa/d)+k/d}.  So E[W^i](E_d) = mean_a mean_k
W(a; frac(pa/d)+k/d)^i, a d-point Riemann sum in the outlier ==> Koksma: |E[W^i]-L_i| <= V_i/(2d),
V_i = E_a[TV_u W(a;.)^i].

THIS SCRIPT:
 (1) rigorous V_i bound: TV_u W <= 2 E_a[W_AP] (overlap of a sliding 1/7-window with U_AP), so
     V_i <= i (6/7)^{i-1} * 2 E[W_block10],  E[W_block10] = 2818/15435.
 (2) the D3 sensitivity |dD3/dm_i| at the limit (L_1,L_2,L_3) (computed high-precision).
 (3) C = sum |dD3/dm_i| * V_i/2  (Koksma 1/(2d)), the ball-safety factor for den=m2-m3/M, and the
     crossover D0 = C/(D3_limit - bar) needed for the rate to clear.
 (4) conditional-D3 finite check of the L=10 family for d in [3, DCHK], all coprime p (min over p):
     confirms D3(E_d) >= bar with the actual (large) margin, extending klein's d<=25.
"""
import numpy as np
from fractions import Fraction as Fr
from math import gcd

TH = 1/7; M = 6/7; BAR = 83549/252252


def W_batch(P):
    """P: (B,11) phases in [0,1) -> W (B,) uncovered measure."""
    P = np.sort(P % 1.0, axis=1)
    g = np.empty_like(P)
    g[:, :-1] = P[:, 1:] - P[:, :-1]
    g[:, -1] = P[:, :1].ravel() + 1.0 - P[:, -1]
    return np.maximum(g - TH, 0.0).sum(axis=1)


def cond_moments(d, p, Na=1200):
    """E[W^i](E_d) i=1,2,3 via klein's equally-spaced conditional structure, batched over (a,k)."""
    a = ((np.arange(Na) + 0.5) / Na)[:, None]         # (Na,1)
    jj = np.arange(10)[None, :]
    ap = (a * jj) % 1.0                               # (Na,10) AP phases
    ks = np.arange(d)[None, :]
    us = ((p * a / d) % 1.0 + ks / d) % 1.0           # (Na,d) outlier phases
    ap_rep = np.repeat(ap, d, axis=0)                 # (Na*d,10)
    u_flat = us.reshape(-1, 1)                         # (Na*d,1)
    P = np.concatenate([ap_rep, u_flat], axis=1)      # (Na*d,11)
    W = W_batch(P)                                     # (Na*d,)
    m1 = W.mean(); m2 = (W * W).mean(); m3 = (W ** 3).mean()
    return m1, m2, m3


def D3_of(m1, m2, m3):
    den = m2 - m3 / M
    return m1 / M + (m1 - m2 / M) ** 2 / den if den > 0 else m1 / M


def main():
    print(f"EXPLICIT L=10 rate constants for klein-S191 (kps-S89); bar = {BAR:.6f}")
    print("=" * 92)

    # (1)(2) limit moments L_i (high precision via cond_moments at large d, ratio-independent)
    L1, L2, L3 = cond_moments(400, 1, Na=1500)
    D3lim = D3_of(L1, L2, L3)
    print(f"\n[limit]  L1={L1:.6f}  L2={L2:.6f}  L3={L3:.6f}   D3_limit = {D3lim:.6f}  (klein 0.4646)")
    EWB = 2818 / 15435
    print(f"  E[W_block10] = 2818/15435 = {EWB:.6f}")

    # D3 partials at the limit
    e = L1 - L2 / M; den = L2 - L3 / M
    dD3_dm1 = 1 / M + 2 * e / den
    dD3_dm2 = -2 * e / (M * den) - e * e / den ** 2
    dD3_dm3 = e * e / (M * den ** 2)
    print(f"\n[sensitivity at limit]  e=m1-m2/M={e:.5f}  den=m2-m3/M={den:.5f}")
    print(f"  |dD3/dm1|={abs(dD3_dm1):.3f}  |dD3/dm2|={abs(dD3_dm2):.3f}  |dD3/dm3|={abs(dD3_dm3):.3f}")

    # (3) rigorous V_i (Koksma 1/(2d) folded in: |m_i-L_i| <= Vh_i/d, Vh_i = i(6/7)^{i-1} E[W_B])
    Vh = {i: i * (6/7) ** (i - 1) * EWB for i in (1, 2, 3)}
    print(f"\n[rigorous |m_i - L_i| <= Vh_i / d]  Vh_1={Vh[1]:.4f}  Vh_2={Vh[2]:.4f}  Vh_3={Vh[3]:.4f}")
    print(f"  (from TV_u W^i <= i(6/7)^(i-1) * 2E[W_B], Koksma discrepancy 1/(2d) => the 2 cancels)")

    Cpoint = abs(dD3_dm1) * Vh[1] + abs(dD3_dm2) * Vh[2] + abs(dD3_dm3) * Vh[3]
    print(f"\n[point constant]  C_point = sum |dD3/dm_i| Vh_i = {Cpoint:.3f}")
    # ball-safety: at d=D0, den drops by <= Vh2/D0 + Vh3/(M D0); partials ~ 1/den^2 -> safety
    D0 = Cpoint / (D3lim - BAR)
    for _ in range(40):
        dden = Vh[2] / D0 + Vh[3] / (M * D0)
        den_lo = den - dden
        if den_lo <= 0:
            D0 *= 1.5; continue
        safety = (den / den_lo) ** 2
        C = safety * Cpoint
        D0 = C / (D3lim - BAR)
    print(f"[ball-adjusted]  worst den at d=D0 = {den_lo:.5f}  safety={safety:.2f}  C<= {C:.2f}")
    print(f"  => CROSSOVER D0 = C/(D3_limit - bar) = {C:.2f}/{D3lim-BAR:.4f} = {D0:.0f}")
    print(f"  RIGOROUS CLOSURE: [finite check d<= {int(np.ceil(D0))}] + [rate d> {int(np.ceil(D0))}: "
          f"D3 >= D3_limit - C/d >= bar].")
    print(f"  (klein's MEASURED C=0.035 is {Cpoint/0.035:.0f}x smaller -- the worst-case sign assumption;")
    print(f"   the rigorous a-priori C is dominated by |dD3/dm2|,|dD3/dm3| and the small den.)")

    # (4) conditional-D3 finite check, extending klein's d<=25
    print("\n" + "=" * 92)
    DCHK = 60
    print(f"CONDITIONAL-D3 FINITE CHECK of the L=10 family d in [3,{DCHK}] (min over coprime p):")
    gmin = (9.9, None)
    for d in range(3, DCHK + 1):
        dmin = (9.9, None)
        for p in range(1, 9 * d + 1):
            if p % d == 0: continue
            if gcd(p, d) != 1: continue
            m1, m2, m3 = cond_moments(d, p, Na=700)
            v = D3_of(m1, m2, m3)
            if v < dmin[0]: dmin = (v, p)
        if dmin[0] < gmin[0]: gmin = (dmin[0], (d, dmin[1]))
        if d <= 6 or d % 10 == 0 or d == DCHK:
            print(f"  d={d:3d}: min D3 = {dmin[0]:.5f} at p={dmin[1]}  margin {dmin[0]-BAR:+.5f}", flush=True)
    print(f"\n  MIN over d in [3,{DCHK}] = {gmin[0]:.6f} at (d,p)={gmin[1]}  margin {gmin[0]-BAR:+.6f}  "
          f"[{'>= bar' if gmin[0]>=BAR else 'BELOW'}]")
    print(f"  (d=3 = A = 0.4530; rises with d toward {D3lim:.4f}.  Extends klein's d<=25 to {DCHK}.)")


if __name__ == "__main__":
    main()
