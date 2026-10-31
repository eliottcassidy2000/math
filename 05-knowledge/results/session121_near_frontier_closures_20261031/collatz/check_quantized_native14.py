#!/usr/bin/env python3
from fractions import Fraction as F
import json

def require(test, message):
    if not test:
        raise RuntimeError(message)

def lam(e):
    return F(5,8)-F(1,8*5**(e-1))

def one_charge(e, length):
    return sum(length//(2*5**(i-1)) for i in range(1,e+1))

eta=F(183105426,48828125)
gamma=eta-1
D=8*5**11
require(D==390625000,'D')
require(D*(F(22,8)-gamma)==342,'exceptional cutoff')
require(gamma>F(21,8),'lower eighth')
q=5**11
pow3=[pow(3,L,q) for L in range(343)]
count=0
for A in range(1,343):
    target=pow(2,A,q)
    for L in range(1,A):
        count+=1
        require(pow3[L]!=target,f'small-clock native multiplier at {(A,L)}')
require(count==58311,'pair count')
r=F(24,25)*F(23,8)+gamma/25-2
require(r==F(1062011676,1220703125),'new fee')
require(r>F(4,5),'minus-v closure')
for e,price in [(11,F(31,125)),(12,F(6,25))]:
    T=2*5**(e-1)
    require(lam(14)*T-one_charge(e,T)==price,'native all-one fee')
    require(r+price>1,'decoration exclusion')
require(2*(1+r)-1>1,'two-core exclusion')
require(F(24,25)*F(25,8)+gamma/25-2>1,'j25 and larger exclusion')
bands=[]
for j,cap,residue in [(23,48828467,3),(24,97656592,0)]:
    require(D*(F(j,8)-gamma)==cap,'clock band')
    require((5*j)%8==residue,'clock residue')
    require(cap%8==residue,'endpoint congruence')
    bands.append({'j':j,'maximum_clock':cap,'clock_mod8':residue})
print(json.dumps({'status':'PASS','inherited_certificate_rerun':False,'small_clock_pairs_checked':count,'primitive13_nonone_fee':str(r),'minus_v_price_inequality':'PROVED conditional on inherited certificate','decorated_candidates':'EXCLUDED','zero_tangent_remaining_bands':bands,'full_E14':'OPEN'},sort_keys=True,indent=2))
