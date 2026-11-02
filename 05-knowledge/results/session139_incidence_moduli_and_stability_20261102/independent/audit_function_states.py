#!/usr/bin/env python3
"""Independent full-function orbit BFS; no parity tags or compressed evaluation vectors."""
from itertools import product
from math import factorial
import json,time

def req(q,m):
    if not q:raise RuntimeError(m)
def gates(f,p):
    P=tuple(x*f[x]%p for x in range(p))
    C=tuple((x*x*f[x]+x*(f[(x+1)%p]-f[(x-1)%p]))%p for x in range(p))
    D=tuple(((x*x*x-2*x)*f[x]+(2*x*x+x)*f[(x+1)%p]+(x-2*x*x)*f[(x-1)%p])%p for x in range(p))
    return P,C,D
def orbit(p):
    f=tuple(range(p));seen={f};queue=[f];i=0
    while i<len(queue):
        f=queue[i];i+=1
        for g in gates(f,p):
            if g not in seen:seen.add(g);queue.append(g)
    zero=(0,)*p;expected=2*p**((p-1)//2)-1
    if p in (5,11):
        req(zero not in seen,('unexpected unary zero',p))
        req(len(seen)==expected-1,('nonzero orbit incomplete',p,len(seen)))
        # Direct full vectors of disjoint signed point-pair supports must be reached.
        a=[0]*p;a[1]=1;a[-1]=p-1
        b=[0]*p;b[2]=1;b[-2]=p-1
        req(tuple(a) in seen and tuple(b) in seen,('missing zero factors',p))
        req(all(x*y%p==0 for x,y in zip(a,b)),('join zero',p))
    else:req(len(seen)==expected and zero in seen,('full unary orbit incomplete',p))
    return {'p':p,'unary_full_functions':len(seen),'zero_in_unary_orbit':zero in seen,'all_context_classes_with_one_join':expected}
def psi(profile,x):
    z=0;t=1
    for c,a in enumerate(profile,1):t*=x-c+1;z+=factorial(c-1)*a*t
    return z
def universal_function_checks():
    n=0
    for M in range(1,33):
        f2s={};s2f={}
        for profile in product((-1,0,1),repeat=5):
            sig=tuple(factorial(c)*factorial(c-1)*a%M for c,a in enumerate(profile,1))
            fun=tuple(psi(profile,x)%M for x in range(M))
            req(fun not in f2s or f2s[fun]==sig,('Psi false merge',M,profile))
            req(sig not in s2f or s2f[sig]==fun,('Psi false split',M,profile))
            f2s[fun]=sig;s2f[sig]=fun;n+=1
            # Recover the actual normalized signature by forward difference.
            for c in range(1,6):
                from math import comb
                d=sum((-1)**(c-j)*comb(c,j)*psi(profile,j) for j in range(c+1))%M
                req(d==sig[c-1],('forward difference',M,c,profile))
    req(pow(3,4,16)!=pow(3,6,16),'ordinary F composite hostile')
    return {'moduli':32,'signed_profiles_per_modulus':243,'function_signature_checks':n,'profile_degree':5,'ordinary_F_composite_hostile':'3^4=1,3^6=9 mod16'}
if __name__=='__main__':print(json.dumps({'status':'PASS','Psi':universal_function_checks(),'orbits':[orbit(p) for p in (3,5,7,11)]},indent=2,sort_keys=True))
