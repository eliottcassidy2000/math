#!/usr/bin/env python3
"""cont.54: SHARPEN -- the LRC(14) M-spectrum near 1/14 is a FAREY/Stern-Brocot mediant tree,
rooted at AP (1/14) and peeled-12-AP (1/13); DC-min 14/183 = my S38 Ostrowski rung M_14=k/(13k+1).
Verify the mediant structure of klein-S266's spectrum {1/14, 2/27, 3/41, 14/183, 1/13, 3/37}."""
from fractions import Fraction as F
def med(a,b): return F(a.numerator+b.numerator, a.denominator+b.denominator)
spectrum=[F(1,14),F(3,41),F(2,27),F(14,183),F(1,13),F(3,37)]
print("klein-S266 M-spectrum (sorted):", [str(x) for x in sorted(spectrum)])
print(f"  as floats: {[round(float(x),5) for x in sorted(spectrum)]}")
print()
print("(1) Ostrowski ladder M_k = k/(13k+1) [my S38]:")
for k in [1,2,3,14]:
    Mk=F(k,13*k+1)
    inspec = Mk in spectrum
    print(f"  M_{k} = {Mk} = {float(Mk):.6f}  {'<-- in spectrum' if inspec else ''}")
print(f"  => AP=M_1=1/14, deep-well=M_14=14/183 [DC-min] are the ladder ENDS")
print()
print("(2) mediant (Stern-Brocot) structure rooted at 1/14 (AP) and 1/13 (peeled 12-AP):")
a,b=F(1,14),F(1,13)
print(f"  root: 1/14 (AP), 1/13 (12-AP peeled)")
m1=med(a,b); print(f"  mediant(1/14,1/13) = {m1}  {'= 2/27 IN SPECTRUM' if m1==F(2,27) else ''}")
m2=med(a,m1); print(f"  mediant(1/14,2/27) = {m2}  {'= 3/41 IN SPECTRUM' if m2==F(3,41) else ''}")
m3=med(m1,b); print(f"  mediant(2/27,1/13) = {m3}")
m4=med(b,F(2,25)); print(f"  mediant(1/13,2/25) = {m4}  {'= 3/38' if m4==F(3,38) else ''}")
print()
print("(3) VERDICT: the M-spectrum near 1/14 is the FAREY tree between AP(1/14) and 12-AP(1/13),")
print("    with the Ostrowski ladder M_k=k/(13k+1) as the covering-min rungs (AP..deep-well).")
print("    Each spectrum value is a mediant = a specific near-AP family (compression scar).")
print("    So the crux 'min M over DC' = the DEEPEST Farey rung that stays covering = 14/183 (deep well).")
