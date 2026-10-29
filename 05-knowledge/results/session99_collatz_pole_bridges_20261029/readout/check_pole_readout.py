#!/usr/bin/env python3
"""Exact pole-order/amplitude clock bridge, with nonconvergence hostiles."""
from math import factorial
import json

def require(b,m):
    if not b:raise RuntimeError(m)
def v2(n):
    require(n>0,'positive valuation input')
    return (n&-n).bit_length()-1

def primitive(k):return factorial(3*k)//(3**k*factorial(k))
def amplitude(k):return primitive(k)*2**(3*k*k+k)
def syr(k,b=3):
    n=b*k+1
    return n//2**v2(n)
def digit_defect(k):return k.bit_count()-(3*k).bit_count()

def main():
    checks=0
    for k in range(301):
        D=primitive(k);Dp=primitive(k+1)
        require(Dp==D*(3*k+1)*(3*k+2),'primitive ladder');checks+=1
        require(v2(D)==2*k+digit_defect(k),'Legendre digit defect');checks+=1
        require(v2(Dp)-v2(D)==v2(3*k+1)+v2(3*k+2),'valuation increment');checks+=1
        if k%2:
            a=v2(Dp)-v2(D)
            require(syr(k)==(3*k+1)//2**a,'clock from amplitude');checks+=1
            require(a==2+digit_defect(k+1)-digit_defect(k),'digital clock');checks+=1
        else:
            require(v2(Dp)-v2(D)==1+v2(3*(k//2)+1),'even guard');checks+=1
    for k in range(51):
        require(amplitude(k+1)==amplitude(k)*2**(6*k+4)*(3*k+1)*(3*k+2),'leading pole amplitude');checks+=1
        if k%2:
            clock=v2(amplitude(k+1))-v2(amplitude(k))-(6*k+4)
            require(clock==v2(3*k+1),'intrinsic amplitude clock');checks+=1
    # Compare pole-index readout with ordinary Collatz, keeping each ordinary step.
    paths={}
    for n0 in (1,3,5,7,11,27,97,871):
        n=n0;path=[]
        for _ in range(500):
            a=v2(3*n+1);m=(3*n+1)//2**a
            raw=3*n+1;halves=[]
            for j in range(a):
                raw//=2;halves.append(raw)
            require(raw==m and raw%2==1,'native clock replay');checks+=1
            # The exponent index is used instead of constructing enormous pole locations.
            require(m==syr(n),'pole-order readout');checks+=1
            path.append([n,a,m]);n=m
            if n==1:break
        require(n==1,'selected finite orbit control')
        paths[str(n0)]={'odd_edges':len(path),'ordinary_clock':sum(1+a for _,a,_ in path),'first_edges':path[:5]}
    # Identical style of meromorphic pole growth does not imply 3x+1 convergence.
    cyc=[13,33,83]
    require([syr(x,5) for x in cyc]==cyc[1:]+cyc[:1],'5x+1 cycle hostile')
    require(syr(3)==5 and amplitude(5)>amplitude(3),'amplitude not Lyapunov')
    require(syr(5)==1 and amplitude(1)<amplitude(5),'descending control')
    require(syr(1)==1,'fixed pole')
    require((3*2+1)//2**(v2(primitive(3))-v2(primitive(2)))!=syr(2),'even-source misuse hostile')
    print(json.dumps({'status':'PASS','checks':checks,'factorial_checks_through':300,'amplitude_checks_through':50,'finite_orbits':paths,'five_x_plus_one_cycle':cyc,'no_claim':'Finite controls and readout conjugacy do not prove Collatz convergence.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
