"""lrc14_momentfloor_bars_check_opus_S186.py (opus-S186)
Verify THM-661's min D3(k) clears each (A') bar momentBar(k)=witnessMP+1-capRat(k) for k=8..13,
and the definitional Bonferroni identity momentBar+capRat-1=witnessMP. Routes hfloor through the
PROVED moment floor (THM-661), retiring the open Lemma A. See LRCWitnessMomentFloor.lean."""
from fractions import Fraction as F
mP=F(14249,252252)
cap={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}
minD3={8:0.675608,9:0.567746,10:0.480,11:0.400,12:0.355876,13:0.308844}  # THM-661 (10,11 approx)
print(f"witnessMP={float(mP):.6f}")
for k in range(8,14):
    bar=mP+1-cap[k]
    print(f"k={k}: momentBar={float(bar):.6f} minD3={minD3[k]:.6f} margin{minD3[k]-float(bar):+.6f}  "
          f"arith(bar+cap-1)={float(bar+cap[k]-1):.8f}=={float(mP):.8f}? {(bar+cap[k]-1)==mP}")
