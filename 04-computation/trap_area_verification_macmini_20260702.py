#!/usr/bin/env python3
"""mac-mini-2026-07-02-S19 (HYP-3876): the trapezoid area = 1/49, speed-independent.
The analytic heart of the drifting pair-floor (klein-S118 handoff). trap(w1,w2,r) is klein's
exact two-tooth overlap (LRCSpreadPairFloor Stage 1). int_R trap dr = 4h^2 = 1/49 for ALL w1<=w2.
Area formula: plateau 4h^2(w2-w1)/w2 + 2 triangles 4h^2 w1/w2 = 4h^2. Confirmed 5 diverse pairs."""
from fractions import Fraction as F
h = F(1,14)
def trap(w1,w2,r): return max(F(0), min(2*h/w2, (h*(w1+w2)-abs(r))/(w1*w2)))
def integral(w1,w2,N=40000):
    S=h*(w1+w2); dr=2*S/N; return sum(trap(w1,w2,-S+dr*(2*i+1)/2)*dr for i in range(N))
print("int trap dr (=1/49, speed-independent):")
for w1,w2 in [(F(10),F(11)),(F(50),F(53)),(F(100),F(140)),(F(23),F(200)),(F(1000),F(1001))]:
    print(f"  w1={w1} w2={w2}: {float(integral(w1,w2)):.6f}  ratio*49={float(integral(w1,w2)*49):.4f}")
print(f"exact 4h^2 = 1/49: {4*h*h==F(1,49)}")
