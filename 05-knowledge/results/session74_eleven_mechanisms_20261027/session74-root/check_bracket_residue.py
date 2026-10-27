#!/usr/bin/env python3
"""Independent exact check of proposed bracket/residue transfer; no producer imports."""
from fractions import Fraction
from math import isqrt
import json

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def prime(p): return p>=2 and all(p%d for d in range(2,isqrt(p)+1))
def direct(q):
    # Lower endpoint excluded, upper included; q*x <= (q-1)*upper.
    return sum(max(0, ((q-1)*(2*n+1)**2)//q-(2*n-1)**2) for n in range(1,2*q+3))
def floor_formula(q):
    T=sum((-((2*n+1)**2))%q for n in range(q))
    return Fraction(16*q*q-24*q+8,3)-Fraction(2*T,q)
def pp_formula(p,e):
    q=p**e
    S=sum(a*(1 if pow(a,(p-1)//2,p)==1 else -1) for a in range(1,p))
    return Fraction(16*q*q-27*q+8,3)+p**(e//2)+Fraction(2*S,p)*((p**((e+1)//2)-1)//(p-1))

cases=[]
for q in range(2,201):
    count=direct(q)
    require(floor_formula(q)==count,('general',q,count,floor_formula(q)))
    # All possible nonempty brackets, checked directly beyond cutoff too.
    positive=[n for n in range(1,2*q+3) if ((q-1)*(2*n+1)**2)//q > (2*n-1)**2]
    require(positive==list(range(1,2*q-1)),('cutoff',q,positive))
for p in range(3,102):
    if not prime(p): continue
    e=1
    while p**e<=2000:
        q=p**e; count=direct(q)
        require(pp_formula(p,e)==count,('prime-power',p,e))
        cases.append([p,e,q,count])
        e+=1
# A separate literal integer set test on a small universe.
for q in range(2,11):
    actual=[]
    for x in range(2,(4*q+5)**2+1):
        n=(isqrt(x-1)+1)//2
        if q*x<=(q-1)*(2*n+1)**2: actual.append(x)
    require(len(actual)==direct(q),('literal',q))
# Negative controls: erase the nonunit stratum or use residue counts without weights.
wrong_zero_count=Fraction(16*121*121-27*121+8,3)+1+Fraction(-22,11)
require(wrong_zero_count==76998 and direct(121)==77008,'omitting nonunit zero roots must fail by ten')
require(sum(a*(1 if pow(a,5,11)==1 else -1) for a in range(1,11))==-11,'QR11 weighted defect')
require(sum(1 if pow(a,5,11)==1 else -1 for a in range(1,11))==0,'unweighted cancellation')
print(json.dumps({'status':'PASS','general_q':[2,200],'literal_q':[2,10], 'prime_power_cases':cases,'eleven_cases':[c for c in cases if c[0]==11], 'controls':['excluded lower endpoint','included upper endpoint','exact bracket cutoff','nonunit stratum required','weighted character not unweighted support']},sort_keys=True,indent=2))
