from area import cycleprofile, invariant
from itertools import product
from math import gcd
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
# (a,d,b,c) abbreviates [[1,a,c],[0,d,b],[0,0,1]].
def mul(u,v,m):
    a,d,b,c=u
    e,f,g,h=v
    return ((e+a*f)%m,d*f%m,(b+d*g)%m,(c+h+a*g)%m)
def ident(r):
    return tuple((0,1,0,0) for _ in range(r))+ (0,)
def mm(u,v,m):
    return tuple(mul(a,b,m) for a,b in zip(u[:-1],v[:-1]))+((u[-1]+v[-1])%m,)
def inv(u,m):
    parts=[]
    for a,d,b,c in u[:-1]:
        e=pow(d,-1,m)
        parts.append((-a*e%m,e,-b*e%m,(a*b*e-c)%m))
    return tuple(parts)+(-u[-1]%m,)
def pp(u,k,m):
    v=ident(len(u)-1)
    while k:
        if k&1:v=mm(v,u,m)
        u=mm(u,u,m)
        k//=2
    return v
def block(w,ts,m):
    parts=[]
    for t in ts:
        z=(0,1,0,0)
        for c in w:z=mul(z,(int(c=='A'),t,int(c=='B'),0),m)
        parts.append(z)
    return tuple(parts)+(1,)
def filt(P,lam,j,m):
    c=[1]
    for i,u in enumerate(lam):
        if i==j:continue
        new=[0]*(len(c)+1)
        d=pow(lam[j]-u,-1,m)
        for k,v in enumerate(c):
            new[k]=(new[k]-u*v*d)%m
            new[k+1]=(new[k+1]+v*d)%m
        c=new
    v=ident(len(P[0])-1)
    for u,k in zip(P,c):
        v=mm(v,pp(u,k,m),m)
    return v
def setup(m):
    if gcd(m,2*3*5*7*13*19*37)!=1:raise ValueError(m)
    half=pow(2,-1,m)
    ts=(2,3*half%m,half)
    X,Y=block('AAB',ts,m),block('ABB',ts,m)
    V=mm(X,inv(Y,m),m)
    C=[V]
    for _ in range(3):C.append(mm(mm(X,C[-1],m),inv(X,m),m))
    lams=[pow(t,e,m) for t in ts[:2] for e in (-3,3)]
    E=[filt(C,lams,j,m) for j in range(4)]
    def comm(a,b):
        return mm(mm(mm(a,b,m),inv(a,m),m),inv(b,m),m)
    return ts,X,Y,comm(E[0],E[1]),comm(E[2],E[3])
def literal(w):
    n=len(w)
    v=1<<n
    r=tuple(v-(1<<(i+1))+(v if c=='B' else 0) for i,c in enumerate(w))
    return r+(sum(1<<i for i,c in enumerate(w) if c=='A'),)
def dp(rows):
    n=len(rows)
    f=(1<<n)-1
    E=[[0]*n for _ in range(f+1)]
    h=[1]+[0]*f
    for s in range(1,f+1):
        for i in range(n):
            if not (s>>i)&1:continue
            t=s^(1<<i)
            E[s][i]=1 if not t else sum(E[t][j] for j in range(n) if (rows[i]>>j)&1 and (t>>j)&1)
        h[s]=sum(E[s])
    return h[f],sum(h[s]*h[f^s] for s in range(f+1))
def main():
    for n in range(1,8):
        for bits in product('AB',repeat=n):
            w=''.join(bits)
            if dp(literal(w))!=cycleprofile(w):raise AssertionError(w)
    print('254 independently counted literal one-apex profiles')
    for m in (11,17,23,121,187,289):
        ts,X,Y,Z,T=setup(m)
        if any(z[1]!=1 or z[0] or z[2] for u in (Z,T) for z in u[:-1]) or Z[-1] or T[-1]:raise AssertionError(m)
        if gcd(Z[0][-1],m)!=1 or gcd(T[1][-1],m)!=1 or Z[1][-1] or T[0][-1] or T[2][-1]:raise AssertionError(m)
        def modprofile(u):
            return ((2+u[0][-1])%m,(16+4*(u[1][-1]+u[2][-1]))%m)
        for a,b in product(range(m),repeat=2):
            u=mm(X,mm(pp(Z,a,m),pp(T,b,m),m),m)
            p,q=modprofile(u)
            if p!=(cycleprofile('AAB')[0]+a*Z[0][-1])%m or q!=(cycleprofile('AAB')[1]+4*(a*Z[2][-1]+b*T[1][-1]))%m:raise AssertionError((m,a,b))
        print('Joint profile plane modulo',m,'size',m*m)
    wa,wb='AABAABB','AAABBAB'
    for s in range(15):
        a,b=('A'*s+w+'B'*s for w in (wa,wb))
        ha,qa=cycleprofile(a)
        hb,qb=cycleprofile(b)
        if invariant(a)!=invariant(b) or (ha-hb,qa-qb)!=(6*2**s,60*6**s+12*2**s):raise AssertionError(s)
    print('15 amplified exact Heisenberg-coordinate collisions')
if __name__=='__main__':main()
