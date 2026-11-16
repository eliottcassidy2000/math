"""Independent mixed-radix contact decoder for the (1203,237,8051) key.
Consumes the literal deletion core bank produced by the canonical C++ generator.
"""
from pathlib import Path

def bits(S):
    while S:
        p=S&-S
        yield p.bit_length()-1
        S^=p

def read(a):
    n=len(a);F=(1<<n)-1
    # Independent path-ending DP by deleting the last vertex.
    I=[F^x^(1<<v) for v,x in enumerate(a)]
    P=[[0]*n for _ in range(F+1)]
    h=[1]+[0]*F
    for s in range(1,F+1):
        for v in bits(s):
            t=s^(1<<v)
            z=sum(P[t][u] for u in bits(t&I[v])) if t else 1
            h[s]+=z
            P[s][v]=z
    assert h[F]==237
    T=[0]*n
    Eo=[0]*n
    for s in range(1,F+1):
        for u in bits(s):
            z=P[s][u]*h[F^s]
            T[u]+=z
            Eo[u]+=z*(s.bit_count()%2)
    # To reuse independent ending DP for starts, reverse all arcs.
    L=[0]*n
    So=[0]*n
    for s in range(1,F+1):
        for v in bits(s):
            t=s^(1<<v)
            z=sum(P[t][u] for u in bits(t&a[v])) if t else 1
            P[s][v]=z
            z*=h[F^s]
            L[v]+=z
            So[v]+=z*(s.bit_count()%2)
    return h,L,So,T,Eo

def local_cut(a,h,Q):
    n=len(a); F=(1<<n)-1
    C=0
    for u in Q:
        P=[[0]*n for _ in range(F+1)]
        P[1<<u][u]=1
        for S in range(1<<u,F+1):
            if S>>u&1:
                for v in bits(S):
                    z=P[S][v]
                    if z==0: continue
                    if S.bit_count()%2==0 and v not in Q:
                        C+=z*h[F^S]
                    for w in bits(a[v]&(F^S)):
                        P[S|(1<<w)][w]+=z
    return C

def absolute_cut(a,h,Q, So, Eo, L,T,R,Y,P):
    p,q=P.bit_length(),len(Y)
    x=[int(b) for b in bin(P)[2:]] if p else []
    c=local_cut(a,h,Q)
    # Both endpoints in one transitive pad, or in opposite pads.
    for z in (x,Y):
        for i,A in enumerate(z):
            if A:
                for j in range(i+1,len(z)):
                    if not z[j]: c+=237*(1<<max(0,j-i-2))
    for i,A in enumerate(x):
        if A:
            k=p-1-i # optional prefix vertices past this endpoint
            for j in range(len(a)):
                if j not in Q: c+= Eo[j] if k==0 else T[j]*(1<<(k-1))
            for j,B in enumerate(Y):
                if not B: c+=(R//2)<<(k+j)
    for i in Q:
        for j,B in enumerate(Y):
            if not B: c+=(So[i] if j==0 else L[i]*(1<<(j-1)))
    return c

def decode(bank, budget):
    by_order={}; suffix={}
    for line in bank.splitlines():
        n,*a=map(int,line.split()); assert n==len(a)
        by_order[n]=by_order.get(n,0)+1
        h,L,So,T,Eo=read(a)
        R=237+sum(L)
        x=R//(R&-R)
        if 8051%x==0: suffix[tuple(a)]=R
        assert R==sum(h[s]*h[len(h)-1^s] for s in range(len(h)))
        for q in range(budget+1):
            limit=(1<<q)*R*(1<<(budget-q))
            X=8051
            while X<=limit:
                if q and X!=8051: break # final retained singleton B forces X odd
                Y=(X*pow(237,-1,1<<q)-1)%(1<<q)
                if q and Y&1: break
                M=(X-237*(Y+1))//(1<<q)
                y=[(Y>>i)&1 for i in range(q-1,-1,-1)]
                for mask in range(1<<n):
                    Q=set(bits(mask))
                    P,rem=divmod(M-sum(L[i] for i in Q),R)
                    if rem or P<0 or P.bit_length()+q>budget: continue
                    if not P and not mask: continue # first component must meet Q
                    if q==0 and len(Q)==n: continue # terminal core must meet B
                    yield (absolute_cut(a,h,Q,So,Eo,L,T,R,y,P), n,P,Y,mask,X,q,tuple(a))
                X*=2
    print('strong H237 bank counts', sorted(by_order.items()))
    print('suffix P2 possibilities at odd phase 8051:',suffix)
    assert not suffix

records=list(decode(Path(__file__).with_name('BANK_237.txt').read_text(),13))
print('phase addresses before even-cut restriction:',len(records))
if records:
    best=min(records)
    print('minimum even cut/row:',best)
    c,n,P,Y,mask,X,q,a=best
    assert c==1232
    p=P.bit_length()
    old=p+n+q
    F=(1<<(old+1))-1
    contact=mask<<p
    for i in range(p):
        contact|=((P>>(p-1-i))&1)<<i
    for i in range(q):
        contact|=((Y>>(q-1-i))&1)<<(p+n+i)
    T=[]
    for i in range(old):
        if i<p or i>=p+n: outgoing=((1<<old)-1)^((1<<(i+1))-1)
        else: outgoing=(a[i-p]<<p)|(((1<<old)-1)^((1<<(p+n))-1))
        T.append(outgoing|((1^((contact>>i)&1))<<old))
    T.append(contact)
    I=[F^x^(1<<v) for v,x in enumerate(T)]
    D=[[0]*(old+1) for _ in range(F+1)]
    h=[1]+[0]*F
    for S in range(1,F+1):
        for v in bits(S):
            R=S^(1<<v)
            D[S][v]=sum(D[R][u] for u in bits(R&I[v])) if R else 1
            h[S]+=D[S][v]
    assert h[F]==237+2*c and h[F^(1<<old)]==237
    for S in range(1,F+1):
        for v in bits(S):
            R=S^(1<<v)
            D[S][v]=sum(D[R][u] for u in bits(R&T[v])) if R else 1
    assert sum(D[S][old]*h[F^S] for S in range(1<<old,F+1))==X
    print('directly replayed least-cut rooted tournament:',T)

print('no reduced signature')
