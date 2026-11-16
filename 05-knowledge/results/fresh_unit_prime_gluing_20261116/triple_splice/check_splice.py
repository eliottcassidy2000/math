#!/usr/bin/env python3
"""Symbolic run compression, plus independent Hamilton subset DP controls."""
from collections import Counter
from itertools import product

C = (30,60,8,16,36,13)
A = (158,252,8,16,36,13,189,60)
B = (222,252,8,16,36,13,188,60)
ROOT=4
U,V=2,3

def require(test, msg):
    if not test:
        raise ValueError(msg)

def validate(a):
    require(all(not(a[i]>>i&1) and a[i] < 1<<len(a) for i in range(len(a))),"bad row")
    require(all((a[i]>>j&1)+(a[j]>>i&1)==1 for j in range(len(a)) for i in range(j)),"not a tournament")

def skeletons(a,first=None,omit=0):
    """Path table by (ordinary vertex mask, U runs, V runs)."""
    n=len(a)
    F=((1<<n)-1)^((1<<U)|(1<<V)|omit)
    out=Counter({(0,0,0):1}) if first is None else Counter()
    def walk(i,S,p,q):
        out[S,p,q]+=1
        for j in range(n):
            if not (a[i]>>j&1):
                continue
            b=1<<j
            if j==U and p<5:
                walk(j,S,p+1,q)
            elif j==V and q<3:
                walk(j,S,p,q+1)
            elif F&b and not S&b:
                walk(j,S|b,p,q)
    for i in range(n) if first is None else [first]:
        if i==U or i==V or F>>i&1:
            walk(i, (1<<i)&F, int(i==U), int(i==V))
    return F,out

def coefficients(a):
    F,K=skeletons(a)
    _,D=skeletons(a,omit=1<<ROOT)
    _,R=skeletons(a,first=ROOT)
    ham={(i,j):c for (m,i,j),c in K.items() if m==F and min(i,j)>0}
    deleted={(i,j):c for (m,i,j),c in D.items() if m==F^(1<<ROOT) and min(i,j)>0}
    uns,sig=Counter(),Counter()
    for (S,i,j),c in R.items():
        for k,l in product(range(6-i),range(4-j)):
            if min(i+k,j+l)>0 and (z:=K[F^S,k,l]):
                uns[i,j,k,l]+=c*z
                sig[i,j,k,l]+=(-1)**S.bit_count()*c*z
    return ham,deleted,uns,sig

def fill(X,sgn,i,j):
    """Unsigned and parity-signed number of fillings of the runs."""
    require(1<=i+j<=2,"unexpected run degree")
    total=1 if i+j==1 else X-2
    if j==0: return total,sgn*total
    if i==0: return total,total
    return total,-1-sgn

def from_coeff(coeff,X,Y,sx,sy):
    ham,deleted,uns,sig=coeff
    h,g=(sum(c*fill(X,sx,i,0)[0]*fill(Y,sy,j,0)[0] for (i,j),c in d.items()) for d in (ham,deleted))
    tot,diff=0,0
    for (i,j,k,l),c in uns.items():
        u,v=fill(X,sx,i,k),fill(Y,sy,j,l)
        tot+=c*u[0]*v[0]
        diff+=sig[i,j,k,l]*u[1]*v[1]
    return h,(tot+diff)//2,(tot-diff)//2,g

def inflate(a,s,t):
    def tt(m):
        return tuple(((1<<m)-1)^((1<<(i+1))-1) for i in range(m))
    return substitute(a,tt(s),tt(t))

def substitute(a,M,N):
    local={U:M,V:N}
    w=[(i,j) for i in range(len(a)) for j in range(len(local.get(i,(0,))))]
    root=w.index((ROOT,0))
    b=tuple(sum(1<<j for j,(v,l) in enumerate(w)
      if (u==v and u in local and local[u][k]>>l&1) or (a[u]>>v&1))
      for u,k in w)
    return b,root

def covers(a,root):
    n=len(a);F=(1<<n)-1
    G=[0]*(F+1)
    H=[0]*(F+1)
    dp=[[0]*n for _ in range(F+1)]
    H[0]=1
    for S in range(1,F+1):
        for u in range(n):
            if not(S>>u&1):continue
            M=S^(1<<u)
            if not M:dp[S][u]=1
            else:dp[S][u]=sum(dp[M][v] for v in range(n) if M>>v&1 and a[u]>>v&1)
        H[S]=sum(dp[S]);G[S]=dp[S][root]
    return H,G

def exact(a,root):
    F=(1<<len(a))-1
    H,G=covers(a,root)
    E,O=0,0
    for S in range(1,F+1):
        m=G[S]*H[F^S]
        if S.bit_count()%2: O+=m
        else: E+=m
    return H[F],E,O,H[F^(1<<root)]

def lookup(X,Y):
    q=3*(2*X-1)*(Y-1)+8
    e=3*(2*X-1)*Y//2
    return q,e,e+3,3

def profile(base,d):
    return tuple(3*x+(d if i==0 else 0) for i,x in enumerate(base))

def main():
    for a in C,A,B:validate(a)
    require(all(tuple(x&63 for x in h[:6])==C for h in (A,B)),"source not induced")
    require(tuple(x^y for x,y in zip(A,B))==(64,0,0,0,0,0,1,0),"not one arc")
    tables=[coefficients(a) for a in (C,A,B)]
    for i,k in enumerate(tables):
        print("CAB"[i], *(sorted(v.items()) for v in k))
    # Each expression is multi-affine in X,Y for fixed size parities. Four
    # interpolation points suffice; do not silently drop any extra run term.
    for sx,sy,X,Y in product((-1,1),(-1,1),(2,4),(2,4)):
        old=lookup(X,Y)
        for tab,ans in zip(tables,(old,profile(old,0),profile(old,2))):
            require(from_coeff(tab,X,Y,sx,sy)==ans,"polynomial identity")
    for s,t in ((1,1),(1,2),(2,1),(2,2),(3,2),(2,3)):
        old=lookup(1<<s,1<<t)
        for a,ans in zip((C,A,B),(old,profile(old,0),profile(old,2))):
            require(exact(*inflate(a,s,t))==ans,"subset DP disagreement")
        print(s,t,old,profile(old,0),profile(old,2))
    def params(a):
        H,_=covers(a,0)
        F=len(H)-1
        P=sum(H[S]*H[F^S] for S in range(F+1))
        Z=sum((-1)**S.bit_count()*H[S]*H[F^S] for S in range(F+1))
        require(Z==0,"alternating cover identity")
        return H[F],P
    for M,N in (((2,4,1),(2,4,1)),((2,4,1),(6,12,8,1)),((2,4,1),(0,))):
        validate(M);validate(N)
        u,R=params(M);v,T=params(N)
        e=3*(2*R-u)*T//2
        old=(3*(2*R-u)*(T-v)+8*u*v,e,e+3*u*v,3*u*v)
        for a,ans in zip((C,A,B),(old,profile(old,0),profile(old,2*u*v))):
            require(exact(*substitute(a,M,N))==ans,"general module formula")
        print("modules",M,N,old,profile(old,0),profile(old,2*u*v))

if __name__=="__main__":
    main()
