#!/usr/bin/env python3
"""Independent exact audit: direct binomial blocks versus modular identities.
No sequence generators, polynomial helpers or formulas imported from parent proof script.
Reproduction: python3 this_file.py [prime_bound=397]
"""
from math import comb
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ordinary_detector_and_mechanism import mechanism
from sympy import Poly, gcd, symbols


def primes(n):
    a=bytearray(b'\1')*(n+1)
    a[:2]=b'\0\0'
    for i in range(2,int(n**.5)+1):
        if a[i]: a[i*i:n+1:i]=b'\0'*len(range(i*i,n+1,i))
    return [i for i in range(3,n+1) if a[i]]
def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a
def add(a,b,p):
    c=[0]*max(len(a),len(b))
    for i,x in enumerate(a):c[i]=(c[i]+x)%p
    for i,x in enumerate(b):c[i]=(c[i]+x)%p
    return trim(c)
def scale(a,b,p):return trim([(x*b)%p for x in a])
def mul(a,b,p):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%p
    return trim(c)
def der(a,p):return trim([i*a[i]%p for i in range(1,len(a))] or [0])
def compose(a,b,p):
    c=[0]
    for x in a[::-1]:c=add(mul(c,b,p),[x],p)
    return c
def direct(p,r):
    u=[comb(2*n,n)**r for n in range(2*p)]
    assert all((u[p+k]-u[1]*u[k])%p==0 for k in range(p))
    f=[u[k]*pow(4**r,-k,p)%p for k in range(p)]
    g=[((u[p+k]-u[1]*u[k])//p)*pow(u[1],-1,p)*pow(4**r,-k,p)%p for k in range(p)]
    return trim(f),trim(g),u

def low_hypergeo(p,numerator,denominator,a):
    h=[1];j=[0];hk=1;logsum=0
    for k in range(1,a+1):
        aa=(denominator*(k-1)+numerator)*pow(denominator,-1,p)%p
        hk=hk*aa*aa*pow(k,-2,p)%p
        logsum=(logsum+2*(pow(aa,-1,p)-pow(k,-1,p)))%p
        h.append(hk);j.append(hk*logsum%p)
    return trim(h),trim(j)

def one(p):
    f,g,u=direct(p,2)
    W=add(mul(f,der(g,p),p),scale(mul(der(f,p),g,p),-1,p),p)
    lhs=add(mul(f,f,p),[0]+W,p)
    rhs=[comb(p-1,k)*(-1)**k%p for k in range(p)]
    assert lhs==trim(rhs),(p,'Wronskian')
    assert sum(f)%p==pow(-1,(p-1)//2,p)
    fc,gc,uc=direct(p,3)
    aa=1 if p%4==1 else 3
    m=(p-aa)//4
    h,j=low_hypergeo(p,aa,4,m)
    ph=[1] if aa==1 else [1,p-1]
    assert fc==mul(ph,mul(h,h,p),p),(p,'Clausen')
    assert gc==mul(ph,mul(h,j,p),p),(p,'log Clausen')
    quarter_w=add(mul(h,der(j,p),p),scale(mul(der(h,p),j,p),-1,p),p)
    quarter_lhs=add(mul(h,h,p),[0]+quarter_w,p)
    assert quarter_lhs==trim([comb(2*m,k)*(-1)**k%p for k in range(2*m+1)]),(p,'quarter Wronskian')
    assert sum(h)%p==comb(2*m,m)%p!=0
    eps=pow(-1,(p-1)//2,p)
    qrev=compose(g,[1,p-1],p)
    bracket=add(add(g,scale(qrev,eps,p),p),scale(f,-eps*sum(g)%p,p),p)
    assert compose(gc,[0,4,p-4],p)==mul(f,bracket,p),(p,'quadratic log gluing')
    # Existing detector receives unnormalized integer sequences, independently.
    leg=mechanism(u,p,False);cube=mechanism(uc,p,True)
    assert leg['squarefree'] and leg['gcd(F,V)=1']
    assert cube['rad|V']
    zz=symbols('z'); FP=Poly(fc[::-1],zz,modulus=p); GP=Poly(gc[::-1],zz,modulus=p)
    radical=FP.exquo(gcd(FP,FP.diff())).monic()
    assert gcd(FP,GP).monic()==radical,(p,'gcd exactly radical')
    return {'p':p,'degree_f2':len(f)-1,'degree_f3':len(fc)-1,'quarter':aa,'degree_factor':m,'legendre':leg,'cube':cube}

def hostile_controls():
    p=7; f,g,_=direct(p,2)
    W=add(mul(f,der(g,p),p),scale(mul(der(f,p),g,p),-1,p),p)
    P=mul([1,p-1],add(mul(f,f,p),[0]+W,p),p)
    assert der(P,p)==[0] and P!=[1] and P==[1]+[0]*6+[6]
    h,j=low_hypergeo(7,3,4,1)
    dropped=mul(h,j,7)
    assert sum(dropped)%7==6
    return {'zero_derivative_not_constant_p7':P,'dropped_singular_prefactor_g_at_1_p7':sum(dropped)%7}

if __name__=='__main__':
    bound=int(sys.argv[1]) if len(sys.argv)>1 else 397
    out={'scope':'All odd primes 3<=p<=bound; direct integer binomial blocks and existing mechanism; not a proof.','bound':bound,'hostile_controls':hostile_controls(),'rows':[one(p) for p in primes(bound)]}
    print(json.dumps(out,indent=2))
