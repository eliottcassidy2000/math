#!/usr/bin/env python3
"""Exact endpoints for the new nonlinear-input cutoff; active under -O."""
from fractions import Fraction as F
from math import factorial
import argparse,json
parser=argparse.ArgumentParser();parser.add_argument('--threshold',type=int,default=140)
N=parser.parse_args().threshold;R=N//2
checks={}
def pos(name,x):
    if not x>0:raise RuntimeError((name,str(x)))
    checks[name]=str(x)
def yes(name,x):
    if not x:raise RuntimeError(name)
    checks[name]=True
pos('R_positive',F(R));pos('binary_bootstrap',F(49,200)-F(4,R))
pos('variance_above4',F(49*R,800)-4)
yes('exp4_partial_sum',sum(F(4)**j/factorial(j) for j in range(8))==F(16319,315)>50)
yes('exp16_comparison',50**4>250000)
cs=[F(249999,1000000),F(7,25),F(57,100),F(99,100)]
bootstrap=[]
for i,(c,cp,k) in enumerate(zip(cs,cs[1:],[2,3,4])):
    m=c*R/4-4-F(2,R)
    pos('bootstrap_mass_%d'%i,m)
    T=sum(m**j/F(factorial(j)) for j in range(k+1))
    pos('bootstrap_remaining_probability_%d'%i,1-cp-F(1,4**R))
    slack=T*(1-cp-F(1,4**R))-1
    pos('bootstrap_step_%d'%i,slack)
    bootstrap.append({'from':str(c),'to':str(cp),'mass':str(m),'Taylor_degree':k,'slack':str(slack)})
c=cs[-1];d=F(19,20);t=d/8;lam=d/64;h=d/16
heavy=c*(N-1)/2-F(16*N,N-1)-d*N/4
pos('heavy_matching_exceeds19',heavy-19)
pos('heavy_matching_derivative',c/2-d/4)
Q=F(R)+F(33,2)+F(8,R);a=F(3,4)
pos('Q_monotone_n',F(1,2)-F(16,(N-1)**2))
pos('Q_even_envelope',Q-F(2*R,2)-F(16*(2*R),2*R-1))
yes('Q_odd_envelope',Q==F(2*R+1,2)+F(16*(2*R+1),2*R))
pos('tail_interval',F(R-20)-1)
pos('first_tail_slack',c-d-Q*a**(R-1))
pos('second_tail_slack',c-d-Q/F(R-20)*a**20)
pos('first_tail_decreasing',1-a*(1+F(1,R)))
yes('second_tail_ratio_identity',Q/F(R-20)==1+(F(73,2)+F(8,R))/F(R-20))
left={(0,2):F(1,4),(1,1):F(1,2),(0,1):-t,(0,0):lam*d,(1,0):-4*lam}
right={(0,2):F(1,4),(1,1):F(1,2),(0,1):-t,(0,0):t*t,(1,0):-t/2}
yes('dual_polynomial_identity',left==right)
yes('h_identity',h==t-4*lam)
pos('strict_final_contradiction',h-F(4,R))
pos('heavy_weights_exceed_kappa',t-F(4,R))
if N==140:
    yes('mass_endpoint_values',[x['mass'] for x in bootstrap]==['969951/2800000','61/70','333/56'])
    yes('heavy_endpoint_value',heavy-19==F(12229,27800))
    yes('strict_final_value',h-F(4,R)==F(1,448))
print(json.dumps({'status':'PASS','threshold':N,'R_min':R,'bootstrap':bootstrap,'checks':checks},sort_keys=True,indent=2))
