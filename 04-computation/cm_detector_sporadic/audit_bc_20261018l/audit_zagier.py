#!/usr/bin/env python3
"""Independent exact recurrence/block/Wronskian audit for six Zagier families."""
from math import comb
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_identities import primes,trim,add,scale,mul,der
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ordinary_detector_and_mechanism import mechanism
PARAMS={'A':(7,2,-8),'B':(9,3,27),'C':(10,3,9),'D':(11,3,-1),'E':(12,4,32),'F':(17,6,72)}

def sequence(a,b,c,N):
    u=[1,b]
    for n in range(2,N+1):
        num=(a*(n-1)*n+b)*u[n-1]-c*(n-1)**2*u[n-2]
        q,r=divmod(num,n*n);assert r==0,(a,b,c,n)
        u.append(q)
    return u

def theta(f,p):return trim([i*x%p for i,x in enumerate(f)])
def shifted(f,k):return [0]*k+f

def run(p,a,b,c,u):
    assert c%p and (a*a-4*c)%p and b%p
    f=trim([u[k]%p for k in range(p)])
    dif=[u[p+k]-b*u[k] for k in range(p)]
    assert all(d%p==0 for d in dif),(p,a,b,c,'Lucas')
    v=trim([d//p%p for d in dif]);g=scale(v,pow(b,-1,p),p)
    h=[1,-a%p,c%p]
    def L(y):
        t=theta(y,p);tt=theta(t,p)
        return add(add(tt,scale(shifted(add(add(scale(tt,a,p),scale(t,a,p),p),scale(y,b,p),p),1),-1,p),p),scale(shifted(add(add(tt,scale(t,2,p),p),y,p),2),c,p),p)
    rhs=add(add(scale(theta(f,p),-2,p),shifted(scale(add(scale(theta(f,p),2,p),f,p),a,p),1),p),scale(shifted(add(theta(f,p),f,p),2),-2*c,p),p)
    assert L(f)==[0] and L(g)==rhs,(p,a,b,c,'inhomogeneous equation')
    w=add(mul(f,der(g,p),p),scale(mul(der(f,p),g,p),-1,p),p)
    R=add(mul(f,f,p),shifted(w,1),p)
    P=mul(h,R,p)
    assert der(P,p)==[0]
    expected=[0]*(2*p+1);expected[0]=1;expected[p]=-a%p;expected[2*p]=c%p
    assert P==trim(expected),(p,a,b,c,'Wronskian with Frobenius term')
    mech=mechanism(u,p,False)
    assert mech['squarefree'] and mech['gcd(F,V)=1']
    return {'p':p,'degF':len(f)-1,'degV':len(v)-1,'Lf_zero':True,'Lg_exact':True,'hR_equals_h_frobenius':True,'mechanism':mech}

if __name__=='__main__':
    bound=int(sys.argv[1]) if len(sys.argv)>1 else 397
    ps=[p for p in primes(bound) if p>=7]
    out={'scope':'All six listed recurrences; every prime7<=p<=bound; exact all-coefficient block equations including terminal coefficients. Finite evidence only.','bound':bound,'families':{}}
    for name,(a,b,c) in PARAMS.items():
        u=sequence(a,b,c,2*bound)
        if name=='A':
            assert all(u[n]==sum(comb(n,k)**3 for k in range(n+1)) for n in range(101))
        if name=='C':
            assert all(u[n]==sum(comb(n,k)**2*comb(2*k,k) for k in range(n+1)) for n in range(101))
        if name=='D':
            assert all(u[n]==sum(comb(n,k)**2*comb(n+k,k) for k in range(n+1)) for n in range(101))
        if name=='E':
            assert all(u[n]==sum(comb(n,k)*comb(2*k,k)*comb(2*n-2*k,n-k) for k in range(n+1)) for n in range(101))
        out['families'][name]={'parameters':[a,b,c],'discriminant':a*a-4*c,'rows':[run(p,a,b,c,u) for p in ps]}
    hostile=sequence(11,3,-1,10)
    hostile_m=mechanism(hostile,5,False)
    assert [x%5 for x in hostile[:5]]==[1,3,4,2,1]
    assert not hostile_m['squarefree'] and not hostile_m['gcd(F,V)=1']
    out['hostile_controls']={'Zagier_D_p5':{'F_coefficients':[x%5 for x in hostile[:5]],'h_factor':'(1+2x)^2','F_factor':'(1+2x)^4','mechanism':hostile_m}}
    print(json.dumps(out,indent=2))
