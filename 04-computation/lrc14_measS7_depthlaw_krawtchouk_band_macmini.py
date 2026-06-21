"""
The DECISIVE structural point: EVERY S_J truncation of the Z/7 coverage IE is a
function of the MARGINAL DEPTH h(x) = #distinct residues hit at x ALONE:
   S_J(E) = integral_0^1  g_J(h(x)) dx,   g_J(h) = sum_{j=0}^J (-1)^j C(7-h, j).
So measS7 = integral g_6(h(x)) dx, g_6(h)=[h=7].

This means: measS7 depends on E ONLY through the law of h(x) -- the DEPTH PROFILE
  pi_E(h) = meas{x : exactly h residues hit}.
THE ENTIRE PROBLEM collapses to: which E maximizes pi_E(7) = P(h=7)?
And every truncation S_J = sum_h pi_E(h) g_J(h) is a LINEAR functional of the depth law.

This is the cleanest possible statement of the crux in spectral terms:
  measS7 is the TOP component of the depth-law vector pi_E in R^8 (h=0..7),
  read against the Krawtchouk-like weight vectors g_J.
"""
from math import comb
print("g_J(h) = sum_{j<=J} (-1)^j C(7-h,j)   [coverage IE weight, depends on h only]")
print(f"{'h':>3} | " + " ".join(f"g{J}".rjust(5) for J in range(7)))
for h in range(8):
    row=[sum((-1)**j*comb(7-h,j) for j in range(J+1)) for J in range(7)]
    print(f"{h:>3} | " + " ".join(str(v).rjust(5) for v in row))
print()
print("g_6(h) = [h==7]  (the indicator).  Each g_J is a Krawtchouk-type partial alternating row.")
print("EVEN columns g_2,g_4,g_6 are >=0 weight vectors that UPPER-bound [h=7] termwise;")
print("ODD columns g_1,g_3,g_5 can dip negative -> overshoot -> Bonferroni LOWER bounds.\n")

# Verify measS7 is a function of depth-law only (two shapes with same depth-law -> same measS7)
from fractions import Fraction as F
def depth_law(E):
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*e+1): bps.add(F(a,7*e))
    bps=sorted(bps); law=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2
        h=len(set(int((e*xm % 1)*7) for e in E))
        law[h]+=hi-lo
    return law
for E in [list(range(8)),[0,2,3,4,5,6,7,8]]:
    law=depth_law(E)
    print(f"  E={E}: depth-law pi[h=0..7]={[float(round(x,4)) for x in law]}  measS7=pi[7]={float(law[7]):.5f}")
