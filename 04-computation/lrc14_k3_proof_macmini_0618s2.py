# k=3 UNCONDITIONAL: any 3 points on the circle have circular max-gap >= 1/3 > 2/7.
# Proof: 3 gaps are nonnegative and sum to 1; by pigeonhole the max gap >= 1/3.
# 1/3 = 0.3333 > 2/7 = 0.2857.  So EVERY x is good (mu=1) for |E|=3, regardless of E.
# Margin: max-gap >= 1/3, threshold 2/7, free margin = 1/3 - 2/7 = 1/21 > 0.
from fractions import Fraction as F
print("k=3: gaps sum to 1, 3 gaps => maxgap >= 1/3 =", float(F(1,3)))
print("     threshold 2/7 =", float(F(2,7)), " margin 1/3-2/7 =", F(1,3)-F(2,7), "=", float(F(1,3)-F(2,7)))
print("     => mu(E)=1 for ALL |E|=3.  PROVED unconditional.")
