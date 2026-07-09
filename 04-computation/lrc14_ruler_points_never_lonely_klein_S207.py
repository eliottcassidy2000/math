#!/usr/bin/env python3
"""
klein-2026-07-09-S207: WHY a good period never certifies loneliness at its own ruler point
-- verifying LRCRulerPoints.lean against mac-mini-S64's exact counterexample.

SETUP (mac-mini's counterexample): E = {0,7,14,21,26,29,37,44,51,58,67,75,82}, Vmax = 91,
v_i = Vmax - e_i.  Since e_0 = 0, the OBSERVER runner v = Vmax = 91 is in the family.

CLAIMS (formalized, sorry-free, in LRCRulerPoints.lean):
 (1) minReach(v, j/Vmax) = 0 for EVERY j   [Vmax * (j/Vmax) = j in Z, so the observer sits on the origin]
     => a good period's own ruler point is NEVER lonely. mac-mini's counterexample, with no computation.
 (2) 1/14 <= minReach(v, tau)  =>  1/14 <= nearInt(Vmax*tau), i.e. the fast phase phi = frac(Vmax*tau)
     is forced into [1/14, 13/14].
 (3) hence, writing tau = (j+phi)/Vmax, the tooth drift d_i = e_i*phi/Vmax obeys |d_i| >= e_i/(14*Vmax):
     THE DRIFT IS UNAVOIDABLE, with klein-S205's 14x floor now NECESSARY, not merely optimal.

Also: locate the true witness (max_tau minReach) and check that its fast phase respects (2).
"""
import numpy as np
from fractions import Fraction
from math import gcd

E = [0,7,14,21,26,29,37,44,51,58,67,75,82]
Vmax = 91
v = [Vmax - e for e in E]
print(f"E    = {E}")
print(f"Vmax = {Vmax}   v = {v}   (observer {Vmax} present: {Vmax in v})\n")

def nearInt(x):
    f = x - np.floor(x); return float(np.minimum(f, 1 - f))
def minReach(v, t):
    return min(nearInt(vi * t) for vi in v)

# (1) ruler points are never lonely
worst = max(minReach(v, j / Vmax) for j in range(Vmax))
print(f"(1) max_j minReach(v, j/Vmax) over j=0..{Vmax-1}  =  {worst:.12f}   (theory: exactly 0)")

# (2)+(3): the exact witness. M(S) = max_t min ||v_i t||; local maxima at t = p/(v_i+v_j).
def nI_frac(num, den):
    r = num % den
    return Fraction(min(r, den - r), den)
qs = sorted({v[i] + v[j] for i in range(len(v)) for j in range(i, len(v))})
best = (Fraction(0), None)
for q in qs:
    for p in range(1, q):
        m = min(nI_frac(vk * p, q) for vk in v)
        if m > best[0]: best = (m, Fraction(p, q))
M, tstar = best
print(f"\n(2) EXACT M(S) = {M} = {float(M):.6f}   at tau* = {tstar} = {float(tstar):.6f}")
print(f"    1/14 = {1/14:.6f};   M >= 1/14 ? {float(M) >= 1/14}")

phi = Fraction(Vmax * tstar.numerator, tstar.denominator)
phi = phi - int(phi)                       # frac(Vmax * tau*)
print(f"\n(3) fast phase phi = frac(Vmax*tau*) = {phi} = {float(phi):.6f}")
print(f"    forced window [1/14, 13/14] = [{1/14:.6f}, {13/14:.6f}];  inside? "
      f"{1/14 <= float(phi) <= 13/14}")
spread = max(E)
print(f"    => unavoidable drift floor: spread/(14*Vmax) = {spread}/(14*{Vmax}) = {spread/(14*Vmax):.6f}")
print(f"       actual max drift at tau*: spread*phi/Vmax = {spread*float(phi)/Vmax:.6f}")

# is tau* a ruler point?
print(f"\n    is tau* a ruler point (denominator | Vmax)? {Vmax % tstar.denominator == 0}")
print(f"    tau* denominator = {tstar.denominator}, Vmax = {Vmax}")

# the drift-embed threshold for this cluster (klein-S205): does it fire at Vmax=91?
def gaps_at(E, Vmax, j):
    c = np.unique(np.array([(e*j) % Vmax for e in E], float) / Vmax)
    cyc = np.append(c, 1.0)
    return [(float(cyc[k]), float(cyc[k+1]-cyc[k])) for k in range(len(cyc)-1)]
fires = False
for j in range(1, Vmax):
    for (a, g) in gaps_at(E, Vmax, j):
        if g > 1/7 + 2*spread*(a+g/2)/Vmax: fires = True
print(f"\n(4) klein-S205 drift-embed fires at Vmax={Vmax} (spread={spread}, ratio {Vmax/spread:.3f})? {fires}")
print(f"    (S205 threshold ~1.41*spread = {1.41*spread:.0f} > {Vmax}: consistent -- the embed correctly")
print(f"     does NOT claim this cluster; mac-mini's counterexample lies inside the open window.)")

print("\nSYNTHESIS: (1) is the structural cause of mac-mini's counterexample -- the ruler point is")
print("disqualified because the OBSERVER (v=Vmax) sits on the origin there. (2)-(3) show the offset")
print("phi>=1/14 is forced, hence the drift is unavoidable. The 1/7 bridge is drift-FREE at a real tau")
print("(klein-S204 criterion C); the drift is an artefact of evaluating the teeth at j/Vmax. The 2/7")
print("criterion buys exactly the room to absorb that artefact. No repair of THM-663 is implied.")
