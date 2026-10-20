#!/usr/bin/env python3
"""Exact rational certificates for logarithms and the spacing/fixed-h constants."""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json
checks=0

def need(x,s):
    global checks
    checks+=1
    if not x:raise ArithmeticError(s)


def log_interval(a,b,N=48):
    t=Q(a-b,a+b);power=t;lower=Q(0)
    for j in range(N):
        lower+=2*power/(2*j+1);power*=t*t
    return lower,lower+2*power/((2*N+1)*(1-t*t))

l2,u2=log_interval(2,1)
l98,u98=log_interval(9,8)
l1716,u1716=log_interval(17,16)
l3,u3=(3*l2+l98)/2,(3*u2+u98)/2
l34,u34=5*l2+l1716,5*u2+u1716
need(u34<Q(18,5),'log34<18/5')
need(2*u2<Q(7,5),'log4<7/5')
need(l3>Q(109,100),'log3>109/100')
spacing=4*Q(18,5)/64+4*Q(1,8)*Q(7,5)-Q(109,100)
need(spacing==-Q(33,200) and spacing<-Q(1,9),'spacing gap constant')
difference=Q(109,100)-Q(1,9)
need(difference==Q(881,900),'fixed-h gap budget')
rows=[]
for h in [3,4,5,8,10,20,100]:
    threshold=720*Q(18+7*h,881)
    N=(threshold.numerator+threshold.denominator-1)//threshold.denominator
    need(N>=32 and 4*Q(18+7*h,5)/N-Q(109,100)<=-Q(1,9),'fixed-h threshold')
    rows.append({'h':h,'N_h':N,'integer_outputs_require_p_below':N**3})
need(64**3>=100,'arithmetic-gap transfer range')
print(json.dumps({'status':'FINITE-EXACT rational constants',
                  'checks':checks,'spacing_upper_bound':str(spacing),
                  'fixed_h_cutoff_formula':'ceil(720*(18+7*h)/881)^3',
                  'fixed_h_examples':rows,
                  'log_bounds':{'log2':[str(l2),str(u2)],'log3':[str(l3),str(u3)],
                                'log34':[str(l34),str(u34)]},
                  'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))
