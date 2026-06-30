"""
THE BRIDGE: rho_j (genuine decorrelation) >= g(O mod 7) (apex skeleton). For 14=2*7 there is ONE odd
prime (7), so the bridge is a SINGLE-prime reduction. Key: the FEJER-BOCHNER MINORANT (min eigenvalue =
the gap) is the right bound; the Z_7*-AVERAGE OVERSHOOTS (klein-S9: Jensen gives mean, not min).
"""
import cmath, math
from itertools import combinations
w=cmath.exp(2j*math.pi/7)
def eigs(O): return [abs(sum(w**(k*x) for x in O))**2 for k in range(1,7)]
def gap(O): return min(eigs(O))      # the MINORANT (Fejer-Bochner)
def avg(O): return sum(eigs(O))/6    # the Z_7*-AVERAGE (Reynolds)
# (1) Z_7*-invariant cores (the only ones where avg=gap trivially)
units=[1,2,3,4,5,6]
inv=[]
for r in range(8):
    for O in combinations(range(7),r):
        O=frozenset(O)
        if all(frozenset((u*x)%7 for x in O)==O for u in units): inv.append(tuple(sorted(O)))
print(f"(1) Z_7*-INVARIANT cores (avg=gap automatically): {inv}")
print(f"    => only {len(inv)} of 128! ({{}}, {{0}}, {{1..6}}, Z_7). The BINDING DOUBLETS are NOT invariant.")
# (2) minorant vs average: the overshoot
print("\n(2) MINORANT (gap) vs Z_7*-AVERAGE for sample cores (avg >= gap always; OVERSHOOTS for non-invariant):")
print(f"    {'core':>14} {'gap (MINORANT)':>15} {'avg (overshoot)':>16} {'invariant?':>11}")
over=0; tot=0
for r in range(1,7):
    for O in combinations(range(7),r):
        O=frozenset(O); tot+=1
        if avg(O)>gap(O)+1e-9: over+=1
for O in [{0,1},{0,1,2},{1,2,4},{0,1,2,3,4},{1,2,3,4,5,6}]:
    O=frozenset(O); iv=all(frozenset((u*x)%7 for x in O)==O for u in units)
    print(f"    {str(tuple(sorted(O))):>14} {gap(O):>15.5f} {avg(O):>16.5f} {str(iv):>11}")
print(f"    => avg > gap (OVERSHOOT) for {over}/{tot} nonempty cores. The AVERAGE is INVALID for the floor;")
print(f"       the floor needs the MIN (worst Fourier mode) = the MINORANT = g(O) >= 4cos^2(3pi/7) (THM-590).")
# (3) the bridge: single odd prime; gap is Z_7*-invariant even when the core isn't
print("\n(3) THE BRIDGE (14 = 2*7, ONE odd prime):")
print("    - the gap g(O) IS Z_7*-invariant (g(uO)=g(O)) even though the doublet core is NOT -- so the")
print("      SKELETON value 4cos^2(3pi/7) is well-defined on the orbit regardless of the averaging.")
print("    - rho_j >= g(O mod 7) is the FEJER-BOCHNER MINORANT bound (decorrelation >= min Fourier coeff),")
print("      NOT the Reynolds average (which overshoots). This is the correct reduction mechanism.")
print("    - 7 is the ONLY odd prime of 14, so there is NO product over primes: rho_j = the apex factor")
print("      alone (mod-2 mirror peeled by THM-580). The bridge is a SINGLE cyclotomic minorant.")
# verify gap Z_7*-invariance on the doublet orbit
d={tuple(sorted(frozenset((u*x)%7 for x in {0,1}))):gap(frozenset((u*x)%7 for x in {0,1})) for u in units}
print(f"\n    doublet {{0,1}} under Z_7*: images {list(d.keys())}, all gap = {set(round(v,5) for v in d.values())}")
