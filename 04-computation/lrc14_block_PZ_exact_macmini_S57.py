"""
mac-mini-2026-07-07-S57 -- EXACT block PZ = E[W]^2/E[W^2] for the covering floor (THM-660,
klein-S177). Settles the THIN k=11 leg (klein: margin +0.016, "an exact block_11 PZ would
settle it"). Cross-checks block_13 = monad's 221828403/815409784.

W(x) = sum_i (g_i(x) - 1/7)_+, g_i = circular gaps of {frac(j x): j=0..k-1} (the block).
EXACT via Farey-cell integration:
 - phase order of {frac(jx)} is constant between breakpoints x = m/d, d=1..k-1 (floor(jx)
   also constant there); on each cell frac(jx) = j x - floor(j x_mid) is LINEAR in x.
 - each circular gap g(x) = alpha + beta x is linear; it crosses 1/7 at x* = (1/7-alpha)/beta.
 - split each cell at those crossings; on each sub-cell W = sum_{gaps>1/7}(g-1/7) is linear;
   integrate W (->quadratic) and W^2 (->cubic) EXACTLY in Fractions.
 => E[W], E[W^2], PZ = E[W]^2/E[W^2] exact rationals.
"""
from fractions import Fraction as F
from math import floor

TH = F(1, 7)

def block_moments(k):
    # primary breakpoints: m/d, d=1..k-1 (phase-order / floor changes), plus 0 and 1
    bps = set([F(0), F(1)])
    for d in range(1, k):
        for m in range(0, d+1):
            bps.add(F(m, d))
    bps = sorted(bps)
    EW = F(0); EW2 = F(0)
    for c in range(len(bps)-1):
        x0, x1 = bps[c], bps[c+1]
        xm = (x0 + x1) / 2
        # phases as linear a_j + b_j x with b_j = j, a_j = -floor(j*xm)
        lin = []
        for j in range(k):
            cj = floor(j * xm)          # floor(j x) constant on (x0,x1)
            lin.append((F(-cj), F(j)))  # frac(jx) = j x - cj  = a + b x
        # sorted order by value at xm
        order = sorted(range(k), key=lambda j: lin[j][0] + lin[j][1]*xm)
        sp = [lin[j] for j in order]    # sorted phases (a,b)
        # circular gaps: g_i = p_{i+1}-p_i (i=0..k-2), g_{k-1} = (1 - p_{k-1}) + p_0
        gaps = []
        for i in range(k-1):
            a = sp[i+1][0]-sp[i][0]; b = sp[i+1][1]-sp[i][1]
            gaps.append((a, b))
        aw = F(1) + sp[0][0] - sp[k-1][0]; bw = sp[0][1] - sp[k-1][1]
        gaps.append((aw, bw))
        # sub-breakpoints where any gap crosses 1/7
        subs = set([x0, x1])
        for (a, b) in gaps:
            if b != 0:
                xs = (TH - a) / b
                if x0 < xs < x1: subs.add(xs)
        subs = sorted(subs)
        for s in range(len(subs)-1):
            u0, u1 = subs[s], subs[s+1]
            um = (u0+u1)/2
            # W(x) = sum over gaps with g(um)>1/7 of (a-1/7 + b x); linear A + B x
            A = F(0); B = F(0)
            for (a, b) in gaps:
                if a + b*um > TH:
                    A += (a - TH); B += b
            # integrate (A+Bx) and (A+Bx)^2 over [u0,u1] exactly
            def i1(A, B, lo, hi):  # int (A+Bx) dx
                return A*(hi-lo) + B*(hi*hi-lo*lo)/2
            def i2(A, B, lo, hi):  # int (A+Bx)^2 dx = A^2 x + A B x^2 + B^2 x^3/3
                return (A*A*(hi-lo) + A*B*(hi*hi-lo*lo) + B*B*(hi**3-lo**3)/3)
            EW += i1(A, B, u0, u1)
            EW2 += i2(A, B, u0, u1)
    return EW, EW2

BARS = {11: F(33121, 100000), 12: F(19934, 100000), 13: F(5649, 100000)}
# exact honest bars (m_P + 1 - min meas): use the rational m_P-derived; klein's decimals as guide
print("=== EXACT block PZ = E[W]^2/E[W^2] (covering floor, THM-660) ===\n")
for k in (11, 12, 13):
    EW, EW2 = block_moments(k)
    PZ = EW*EW/EW2
    print(f"k={k}: E[W] = {EW} = {float(EW):.6f}")
    print(f"      E[W^2] = {EW2} = {float(EW2):.6f}")
    print(f"      PZ = {PZ} = {float(PZ):.6f}")
    print(f"      (7/6)E[W] = {float(F(7,6)*EW):.6f};  klein bar ~ {float(BARS[k]):.5f}; "
          f"margin {float(PZ)-float(BARS[k]):+.5f}  {'CLEARS' if PZ>BARS[k] else 'BELOW'}")
    if k == 13:
        mon = F(221828403, 815409784)
        print(f"      cross-check vs monad PZ-on-V {float(mon):.6f}: {'MATCH' if abs(float(PZ)-float(mon))<1e-4 else 'DIFF (block vs AP-min shape)'}")
    print()
