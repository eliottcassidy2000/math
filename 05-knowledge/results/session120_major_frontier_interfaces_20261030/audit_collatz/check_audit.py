"""Independent scalar-envelope audit; no producer imports or claims of native attainment."""
from fractions import Fraction as F
import json


def require(value, message):
    if not value:
        raise RuntimeError(message)

D = 8 * 5**10
gamma = F(5878847, 19531250)
l11 = F(12207031,19531250)
l12 = F(5,8) - F(1,5*D)
require(l11 == F(5,8)-F(1,D), 'lambda11 units')
require(l12-l11 == F(1,2*5**11), 'increment units')
require(F(1,4)<gamma<F(3,8), 'eighth threshold')

def upper(j):
    cap = (D*(F(j,8)-gamma)).__floor__()
    residue = 5*j % 8
    return cap - (cap-residue)%8

def fee(j):
    return F(j,8)-F(upper(j),5*D)

A = upper(3)
S = F(5*A-3,8)
kappa = fee(3)
require(A==5781487 and S==3613429, 'sharp scalar pair')
require(S.denominator==1 and S <= l11*A-gamma, 'pair admissibility')
require(S > l11*(A+8)-gamma-5, 'next clock violates inherited bound')
require(kappa==F(17587861,48828125), 'sharp scalar fee')
coarse=(F(3,2)+gamma)/5
require(kappa>=coarse and kappa < (2+gamma)/5, 'universal j>=4 separation')
require(3*kappa-1==F(3935458,48828125)>0, 'three-fee closure')
# Universal proof is in AUDIT.md; these are independent exact regression checks.
for j in range(3,10003):
    U=upper(j)
    require((5*U-j)%8==0,'clock congruence')
    require(F(j,8)>=gamma+F(U,D),'upper endpoint feasible')
    require(F(j,8)<gamma+F(U+8,D),'next endpoint infeasible')
    require(fee(j)>=kappa,'envelope lower bound')
    require(upper(j+8)==U+D,'periodic clock envelope')
    require(fee(j+8)==fee(j)+F(4,5),'periodic fee envelope')
# The improved uniform fee alone does not close any of the other six families.
remaining=[(1,0),(1,1),(1,2),(1,3),(2,0),(2,1)]
for b,a in remaining:
    require(b*kappa+F(a,5)<1,'remaining uniform-price obstruction')
print(json.dumps({'status':'PASS','conditional_dependency':'pinned E11+ and native return-kernel criterion',
 'D':D,'gamma':str(gamma),'coarse_fee':str(coarse),'sharp_scalar_fee':str(kappa),
 'scalar_pair_only':{'A':A,'S11':int(S),'native_word_attainment_claimed':False},
 'three_fee_excess':str(3*kappa-1),'regression_j_count':10000,
 'remaining_families_b_a':remaining,'E12':'OPEN','universal_Collatz':'OPEN'},indent=2))
