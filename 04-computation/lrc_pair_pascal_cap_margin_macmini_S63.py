#!/usr/bin/env python3
"""The covering-bound cap is a PAIR-NORMALIZED PASCAL MASS, with a higher-moment dip at k=8,9 (S63).

Merge of the owner's 'pair-normalized pascal mass' with codex HYP-3090 (cap_k = C(k+1,2)/C(14,2)) and
mac-mini's gK8 (binding = the pairwise S2 moment).

  cap_k = C(k+1,2)/C(14,2) = C(k+1,2)/91   EXACTLY for k = 10,11,12,13   (codex, verified)
        = [the second factorial moment of the block-occupancy: P(a random pair of the 14-clock lies
           inside the (k+1)-point block)]  -- a PAIR-NORMALIZED PASCAL MASS (hypergeometric C(k+1,2)/C(14,2))
  cap_8, cap_9 sit a small DIP below the pair-Pascal mass  (1081/76440 at k=8, 1/4004 at k=9).

Reading: C(k+1,2)/91 is the degree-2 (pairwise / Krawtchouk j=2) value; the dip at the SPARSE binding
rows k=8,9 is the higher-Pascal (j>=3) tightening = exactly the gK8 -9 S3 + 6 S4 correction (mac-mini S60).
The covering-bound MARGIN is then a clean pair-complement:
    1 - cap_k = (C(14,2) - C(k+1,2))/C(14,2) = (# pairs OUTSIDE the (k+1)-block)/91     [exact for k>=10]
"""
import sys
from fractions import Fraction as F
from math import comb
sys.stdout.reconfigure(line_buffering=True)

caps = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91), 11: F(66,91), 12: F(6,7), 13: F(1)}
C142 = comb(14, 2)  # 91

print("=" * 84)
print(" THE COVERING-BOUND CAP AS A PAIR-NORMALIZED PASCAL MASS  cap_k = C(k+1,2)/C(14,2)")
print("=" * 84)
print(f"{'k':>3} | {'cap_k (true)':>14} | {'C(k+1,2)/91 (pair-Pascal)':>26} | {'dip = tri - cap':>16} | {'margin 1-cap':>14}")
print("-" * 84)
for k in sorted(caps):
    tri = F(comb(k+1, 2), C142)
    dip = tri - caps[k]
    margin = 1 - caps[k]
    pair_complement = F(C142 - comb(k+1, 2), C142)
    tag = "EXACT" if dip == 0 else "dip"
    print(f"{k:>3} | {str(caps[k]):>14} | {str(tri):>20} {tag:>5} | {str(dip):>16} | {str(margin):>14}")
print("-" * 84)
print(" margin 1-cap_k = (91 - C(k+1,2))/91 = # pairs OUTSIDE the (k+1)-block / 91  [exact for k>=10]:")
for k in sorted(caps):
    print(f"   k={k:2d}: 1-cap={float(1-caps[k]):.5f}   pair-complement (91-C(k+1,2))/91 = "
          f"{F(C142-comb(k+1,2),C142)} = {float(F(C142-comb(k+1,2),C142)):.5f}"
          f"   {'==' if (1-caps[k])==F(C142-comb(k+1,2),C142) else '(dip)'}")

print("\n" + "=" * 84)
print(" THE DIP = the higher-Pascal (Krawtchouk j>=3) tightening at the SPARSE binding rows")
print("=" * 84)
print(" The pair-Pascal mass C(k+1,2)/91 is the degree-2 (pairwise) value -- EXACT once the config is")
print(" dense (k>=10). For the sparse rows k=8,9 the true cap dips below it; the dip is the cubic/quartic")
print(" correction = exactly gK8's -9 S3 + 6 S4 (mac-mini S60: k=8 is the S2/S3/S4-balance row).")
print(f"   dip_9 = {caps[9]} vs {F(comb(10,2),91)}  ->  {F(comb(10,2),91)-caps[9]} = 1/4004 = 1/(44*91)")
print(f"   dip_8 = {F(comb(9,2),91)-caps[8]} = {float(F(comb(9,2),91)-caps[8]):.6f}  (the largest dip = the hardest row)")
print("\n PAIR-NORMALIZED PASCAL MASS as a probability (the hypergeometric reading):")
print("   cap_k = P(a uniformly random pair {i,j} of the 14 clock-positions both lie in a (k+1)-block)")
print("         = C(k+1,2)/C(14,2)  = the SECOND factorial moment of the block-occupancy indicator.")
print(" So the covering bound is 'pairwise occupancy <= pairwise capacity', a pair-Pascal inequality;")
print(" the binding-row obligation is the finite higher-Pascal dip (k=8,9). [ties HYP-3085 gK8, HYP-2716 Krawtchouk]")
