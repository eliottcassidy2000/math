#!/usr/bin/env python3
"""Exact finite controls for the PPT elliptic positivity classification."""
from fractions import Fraction as F
import json
from math import gcd
from verify_seeded_positive_elliptic_paths_20261019 import add

def R(t): return t**6+4*t**5+3*t**4-8*t**3+35*t**2-12*t+1
def G(t): return t**3+6*t**2+14*t-2
checks=0
def check(ok):
    global checks
    checks+=1
    if not ok: raise RuntimeError('positivity threshold control failed')
l,h=F(13488,100000),F(13489,100000)
check(G(l)<0<G(h));check(R(l)<0 and R(h)<0)
check(R(F(1,5))==F(-904,15625));check(R(F(23,100))>0)
for t in [-1,1]: check(R(F(t))!=0)
for t in [-2,-1,1,2]: check(G(F(t))!=0)
a,b,c=427,1836,1885;t=F(a,b)
check(a*a+b*b==c*c and gcd(a,b)==1)
check(t>F(23,100));check(t*t+4*t-1<0);check(t*t-4*t+1>0)
A=2*(a*a+b*b-4*a*b);Y=16*a*a*b;U=4*Y*(b-a)-A*A;V=Y*Y-A*U
P=(F(A),F(Y));Q=None;points=[]
for n in range(1,10):
 Q=add(Q,P,(-A,U,V));points.append(Q)
 check(Q[0]>0 if n<9 else Q[0]<0)
check(points[-1][0]==F(-126100899082622355694,917204544598225))
print(json.dumps({'status':'FINITE-EXACT controls; infinite proof in companion note', 'checks':checks,
 'rho_bracket':list(map(str,[l,h])), 'R_bracket_values':list(map(str,[R(l),R(h)])),
 'counterexample_ppt':[a,b,c], 'abscissas':list(map(lambda Q:str(Q[0]),points))},indent=2,sort_keys=True))
