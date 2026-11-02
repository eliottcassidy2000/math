#!/usr/bin/env python3
"""Exact endpoint audit for PROOF_N526.md; active under python -O."""
from fractions import Fraction as F
from math import factorial
import json
N=526
R=263
c=F(249999,1000000)
d=F(487,2000)
t=d/8
lam=d/64
h=d/16
a=F(3,4)
checks={}
def pos(name,value):
    if not value>0:
        raise RuntimeError((name,value))
    checks[name]={"positive_exact":str(value)}
def yes(name,value):
    if not value:
        raise RuntimeError(name)
    checks[name]=True
s=sum((F(4)**j)/factorial(j) for j in range(8))
yes('partial_exp4_value',s==F(16319,315))
pos('partial_exp4_exceeds50',s-50)
pos('c0_lower_exceeds_kappa_upper',F(49,200)-F(4,R))
pos('V_lower_exceeds4',F(49*(N-1),1600)-4)
yes('exp16_lower_exceeds250000',50**4>250000)
yes('c_value_from_250000',c==F(1,4)*(1-F(1,250000)))
heavy_lower=c*(N-1)/2-F(16*N,N-1)-d*N/4
pos('heavy_count_lower_minus17',heavy_lower-17)
yes('heavy_count_endpoint_value',heavy_lower-17==F(964339,1680000))
pos('heavy_count_derivative_constant',c/2-d/4)
Q=F(R)+F(33,2)+F(8,R)
yes('Q_odd_n_identity',Q==F(2*R+1,2)+F(16*(2*R+1),2*R))
pos('Q_even_n_dominance',Q-F(2*R,2)-F(16*(2*R),2*R-1))
pos('Q_n_derivative',F(1,2)-F(16,(N-1)**2))
yes('tail_interval_nonempty',R-18>=1)
pos('first_tail_slack',c-d-Q*a**(R-1))
slack=c-d-Q/F(R-18)*a**18
pos('second_tail_slack',slack)
yes('second_tail_slack_value',slack==F(1846127022469099,27674621771776000000))
pos('first_tail_ratio_slack',1-a*(1+F(1,R)))
yes('second_tail_ratio_identity',Q/F(R-18)==1+(F(69,2)+F(8,R))/F(R-18))
# Polynomial identity coefficients in (w,S), represented by exponent tuples.
left={(0,2):F(1,4),(1,1):F(1,2),(0,1):-t,(0,0):lam*d,(1,0):-4*lam}
# Expansion of (S/2-t+w/2)^2+w*(2t-w)/4.
right={(0,2):F(1,4),(1,1):F(1,2),(0,1):-t,(0,0):t*t,(1,0):-t/2}
yes('dual_polynomial_identity',left==right)
yes('dual_remainder_nonnegative_on_light_interval',t>0)
yes('h_identity',h==t-4*lam)
pos('dual_contradiction',h-F(4,R))
yes('dual_contradiction_value',h-F(4,R)==F(81,8416000))
pos('heavy_weights_exceed_kappa',t-F(4,R))
# Deliberately weaker n=525 has R=262; these constants must not certify it.
yes('negative_control_R262_fails',h-F(4,262)<0)
print(json.dumps({'status':'PASS','N':N,'R_min':R,'constants':{'c':str(c),'d':str(d),'t':str(t),'lambda':str(lam),'h':str(h)},'checks':checks},sort_keys=True,indent=2))
