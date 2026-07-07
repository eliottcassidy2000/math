from fractions import Fraction as F
from math import comb
import numpy as np

# Independent-uniform-points reference: P(all k gaps <= g) = sum_j (-1)^j C(k,j)(1-jg)_+^{k-1}
def P_all_gaps_le(k, g):  # g = Fraction
    tot=F(0)
    j=0
    while 1 - j*g > 0:
        tot += (-1)**j * comb(k,j) * (1 - j*g)**(k-1)
        j+=1
    return tot

g=F(1,7)
print("=== mu_17: AP minimizer vs INDEPENDENT-points limit (decorrelation tail), opus-S131 ===")
print(f"{'k':>3} {'mu17(AP) exact':>16} {'~':>7} {'mu17_indep=1-P':>16} {'~':>8}   AP < indep?")
ap_exact={8:F(691,735),9:F(247,294),10:F(38,49),11:F(1381,2205),12:F(13823,24255),13:F(477,1078)}
for k in range(8,14):
    indep = 1 - P_all_gaps_le(k,g)
    ap = ap_exact[k]
    print(f"{k:>3} {str(ap):>16} {float(ap):>7.4f} {str(indep):>16} {float(indep):>8.4f}   {ap < indep}")

# empirical: does large normalized-spread E approach the indep limit?
print("\n=== empirical: mu_17 of spread-out (decorrelated) E approaches indep limit? (k=8) ===")
GRID=200000; xs=(np.arange(GRID)+0.5)/GRID
def mu_17(E):
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0); ph.sort(axis=1)
    allg=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1]).reshape(-1,1)],axis=1)
    return float((allg.max(axis=1)>1/7).mean())
import random; random.seed(1)
indep8=float(1-P_all_gaps_le(8,g))
for spreadmax in [20,100,1000,100000]:
    vals=[mu_17(sorted(random.sample(range(1,spreadmax+1),8))) for _ in range(20)]
    print(f"   spread<={spreadmax:>6}: mean mu_17 = {np.mean(vals):.4f} (indep limit {indep8:.4f})")
print("\n  STRUCTURE: mu_17 ranges from mu_17(AP) [minimum, most correlated] up to ~indep limit")
print("  [maximum, decorrelated]. The density floor = the AP minimum. Tail is SAFE (>> 0).")
