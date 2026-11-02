#!/usr/bin/env python3
"""Independent residue-block audit of the stated affine-period laws."""
from fractions import Fraction as F
from math import gcd,lcm
import json


def need(x,msg):
    if not x:raise RuntimeError(msg)


def values(q,A,h,m):
    d=next((b for b in range(1,q) if b not in A),q);n=q+h
    E=[u+q*(u in A) for u in range(1,q) if gcd(u,q)==1]
    units=[u for u in range(1,n) if gcd(u,n)==1]
    w={u:F(1,pow(u,-1,n))+F(1,n-pow(u,-1,n)) for u in units}
    N=0; scaled=F(0)
    for e in E:
        K=m*h//e;blocks,rem=divmod(K,n)
        N+=2*(blocks*len(units)+sum(u<=rem for u in units))
        scaled+=F(2,q)*(blocks*sum(w.values(),F(0))+sum((z for u,z in w.items() if u<=rem),F(0)))
    if h==d and gcd(q*m,n)==1:
        N+=len(units)-2
        H=sum((F(1,u) for u in units),F(0))
        scaled+=F(2,q)*(H-w[m%n])
    return N,scaled,E,d


def main():
    cases=boundary_cases=0
    for q in range(3,13):
        for mask in range(1<<(q-1)):
            A={u for u in range(1,q) if mask>>(u-1)&1}
            d=next((b for b in range(1,q) if b not in A),q)
            for h in range(1,d+1):
                n=q+h
                _,_,E,_=values(q,A,h,1)
                L=n*lcm(*E);HE=sum((F(1,e) for e in E),F(0))
                units=[u for u in range(1,n) if gcd(u,n)==1]
                HU=sum((F(1,u) for u in units),F(0))
                s=min(E)
                for m in (s,s+1,s+5):
                    N,B,_,_=values(q,A,h,m);N2,B2,_,_=values(q,A,h,m+L)
                    need(N2-N==2*L*h*len(units)*HE/n,('count period',q,A,h,m))
                    need(B2-B==4*L*h*HE*HU/(q*n),('scaled mass period',q,A,h,m))
                    cases+=1;boundary_cases+=h==d
    print(json.dumps({'status':'FINITE-EXACT / PASS','q':[3,12],'all_attachment_sets':4092,
                      'm_choices':'s,s+1,s+5','all_h':'1..d','period_cases':cases,
                      'boundary_correction_period_cases':boundary_cases},indent=2,sort_keys=True))

if __name__=='__main__':main()
