#!/usr/bin/env python3
from math import gcd
from census import table

def infer(target, a, b, c, d, fixed_h, fixed_p):
    U,V=a*fixed_h+b*fixed_p,c*fixed_h+d*fixed_p
    if not (a>0 and d>0 and b*c>a*d and fixed_h>0 and fixed_p>=0):
        raise ValueError('require a physical 2-visit hyperbolic form and fixed module profile')
    if target<U or target%2==0:return []
    g=gcd(U,2*V)
    if target%g:return []
    mod=2*V//g
    residue=target//g*pow(U//g,-1,mod)%mod
    # Q=p+2h>=2^n and p>=2^n−2.
    max_n=(2+(target-U)//V).bit_length()-1
    return [(n,h,p,z)for n,h,p,z,_ in table()if n<=max_n and h%mod==residue and U*h+V*p==target]

if __name__=='__main__':
    # H=495 is the first n<=9 leading-"11" example for fixed singleton.
    assert infer(495,11,6,6,3,1,0)
    print('N=495:',infer(495,11,6,6,3,1,0))
