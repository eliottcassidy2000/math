#!/usr/bin/env python3
"""Symbolic run coefficients, and independent subset DP."""
from collections import Counter
from itertools import product

C=(30,60,24,16,32,13)
A=(94,124,24,80,32,77,20)
B=(94,60,88,16,32,13,58)
M,ROOT=4,5

def require(ok,label):
    if not ok: raise ValueError(label)

def dp(a,start=None):
    n=len(a);N=1<<n
    H=[1]+[0]*(N-1)
    D=[[0]*n for _ in range(N)]
    for S in range(1,N):
        for i in range(n):
            if not(S>>i&1):continue
            T=S^(1<<i)
            D[S][i]=sum(D[T][j] for j in range(n) if T>>j&1 and a[j]>>i&1) if T else int(start is None or i==start)
        H[S]=sum(D[S])
    return H

def state(a,r):
    H,P=dp(a),dp(a,r); F=len(H)-1
    E=sum(P[S]*H[F^S] for S in range(1,F+1) if S.bit_count()%2==0)
    O=sum(P[S]*H[F^S] for S in range(1,F+1) if S.bit_count()%2==1)
    return H[F],E,O,H[F^(1<<r)]

def validate(a,strong=True):
    n=len(a)
    require(all((a[i]>>j&1)+(a[j]>>i&1)==1 for i in range(n) for j in range(i)),"invalid arcs")
    require(all(not(r>>i&1) and r<1<<n for i,r in enumerate(a)),"invalid row")
    if strong:
        for step in (lambda x:a[x],lambda x:sum(1<<j for j in range(n) if a[j]>>x&1)):
            S=1
            for _ in a:
                for i in range(n):
                    if S>>i&1:S|=step(i)
            require(S==(1<<n)-1,"non-strong")

def slots(a,drop=0,r=None):
    K=((1<<len(a))-1)^((1<<M)|drop)
    out=Counter({(0,0):1}) if r is None else Counter()
    def rec(i,S,k):
        out[S,k]+=1
        for j in range(len(a)):
            if not a[i]>>j&1:continue
            if j==M and k<3:rec(j,S,k+1)
            elif K>>j&1 and not S>>j&1:rec(j,S|(1<<j),k)
    for i in range(len(a)) if r is None else (r,):
        if K>>i&1 or i==M:rec(i,(1<<i)&K,int(i==M))
    return K,out

def cert(a):
    F,H=slots(a); G,Z=slots(a,drop=1<<ROOT); _,P=slots(a,r=ROOT)
    ham={k:c for (S,k),c in H.items() if S==F and k}
    g={k:c for (S,k),c in Z.items() if S==G and k}
    U,D=Counter(),Counter()
    for (S,i),v in P.items():
        for j in range(4-i):
            if i+j and (w:=H[F^S,j]):
                U[i,j]+=v*w; D[i,j]+=(-1)**S.bit_count()*v*w
    return ham,g,U,D

def tt(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def inflate(a,B):
    v=[(i,j) for i in range(len(a)) for j in range(len(B) if i==M else 1)]
    return tuple(sum(1<<t for t,(j,l) in enumerate(v)
      if a[i]>>j&1 or (i==j==M and B[k]>>l&1)) for i,k in v),v.index((ROOT,0))

def packet():
    for a,h,g,f in ((C,15,1,1),(A,45,3,3),(B,49,3,3)):
        validate(a)
        require(cert(a)==({1:h,2:11*f},{1:g},{(0,1):12*f,(1,0):11*f,(1,1):11*f},
          {(0,1):0,(1,0):f,(1,1):f}),"skeleton mismatch")
    for old in (A,B):require(tuple(r&63 for r in old[:6])==C,"not induced")
    for v in (tt(1),tt(2),tt(3),(2,4,1),(6,12,8,1)):
        H=dp(v); F=len(H)-1; u=H[F]; R=sum(H[S]*H[F^S] for S in range(F+1))
        profile=(11*R-7*u,11*R//2,11*R//2+u,u)
        for a,f,delta in ((C,1,0),(A,3,0),(B,3,4*u)):
            require(state(*inflate(a,v))==tuple(f*x+(delta if i==0 else 0) for i,x in enumerate(profile)),
                "arbitrary module DP")
    print("0/+4 run table and controls pass")

def from_words(left,k,right,cycle=True):
    # Old vertices in SCC order, then the root.
    b=len(left); m=3 if cycle else 0
    w=left+((1,)*k+(0,)*(m-k))+right
    n=len(w)
    a=list(tt(n))
    if cycle:
        a[b]^=1<<(b+2)
        a[b+2]|=1<<b
    for i,h in enumerate(w):
        if not h:a[i]|=1<<n
    a.append(sum(h<<i for i,h in enumerate(w)))
    validate(a)
    return a,n

def strong_g3_contacts(J,D):
    """Return (prefix bits,k,suffix bits) for strong g=3 at the specified J,D."""
    if J<1:return []
    ans=[]
    for b in range(J.bit_length()):
        c,V=divmod(J,1<<b)
        U,k=divmod(c,4)
        left=tuple(map(int,bin(U)[2:])) if U else ()
        right=tuple(map(int,f"{V:0{b}b}")) if b else ()
        if not left and k==0:continue  # root must reach first SCC
        if (right and right[-1]) or (not right and k==3):continue  # last SCC must reach root
        if D!=(-3+3*(V&1) if b else -3+k):continue
        ans.append((left,k,right))
    return ans

def decoder_grid():
    Z=lambda n:(0,)*n
    for a,b in product(range(3),repeat=2):
        s=a+b+2; q=11*2**s-7; E=33*2**(s-1)
        require(state(*from_words((1,0,1,1)+Z(a),0,Z(b)))==(3*q,E,E+3,3),"two-block grid")
    for s in (1,2,3):
        q=11*2**s-7; E=33*2**(s-1)
        for l,k,r,d in (((1,0,1),2,Z(s-1),-16),((1,0),3,Z(s),0),
          ((1,),1,(1,)+Z(s),4),((),2,(1,1)+Z(s),-4),((),1,(0,1,1)+Z(s),0)):
            if s==1 and d==-16:
                got=state(*from_words(l,k,r))
                require(got==(27,E+1,E+2,3),"excluded boundary")
                continue
            require(state(*from_words(l,k,r))==(3*q+d,E,E+3,3),"one-block grid")
        got=Counter(state(*from_words(*x))[0] for x in strong_g3_contacts(11*2**s,-3))
        want=Counter({3*q:(2 if s==1 else s+1),3*q+4:1,3*q-4:1})
        if s>=2:want[3*q-16]=1
        require(got==want,"complete strong decoder")
        require(state(*from_words((1,0,1,1)+Z(s),0,(),False))==(q,E//3,E//3+1,1),
           "g=1 normalization")
    print("all coefficient grids and the exceptional boundary pass")

if __name__=="__main__":
    packet()
    decoder_grid()
