"""Exact tournament counters and quadratic cut enumeration, standard library only.

The run-assignment transfer is inherited from Session191; attachment
quadraticity is inherited from Session203. See PROOF.md for scope.
"""
from itertools import product

def count_order(a,b,ab):
    # Runs I,a,J,b,K, 0=empty; nonzero status for valid latest end.
    state={(0,0,0):1}
    for ai,bi in zip(a,b):
        out={}
        for (i,j,k),v in state.items():
            # ai=1 iff a->x; bi=1 iff b->x
            t=(1 if not ai else 2,j,k); out[t]=out.get(t,0)+v
            if j or ai:
                t=(i,1 if not bi else 2,k); out[t]=out.get(t,0)+v
            if k or bi:
                t=(i,j,1); out[t]=out.get(t,0)+v
        state=out
    return sum(v for (i,j,k),v in state.items() if i!=2 and j!=2 and (j or ab))
def count(a,b,ab=1):return count_order(a,b,ab)+count_order(b,a,1-ab)
def coeff(L,r,ab):
    a=[int(i==r) for i in range(L)]; z=[0]*L
    c=count(a,z,ab); lin=[]; quad={}
    for i in range(L):
        b=z.copy();b[i]=1;lin.append(count(a,b,ab)-c)
    for i in range(L):
        for j in range(i+1,L):
            b=z.copy();b[i]=b[j]=1;quad[i,j]=count(a,b,ab)-c-lin[i]-lin[j]
    return c,lin,quad

def scan(L,r,ab,cap):
    base,lin,quad=coeff(L,r,ab);h=base;d=lin.copy();old=0
    for k in range(1<<L):
        q=k^(k>>1)
        if k:
            i=(q^old).bit_length()-1;s=1 if q>>i&1 else -1;h+=s*d[i]
            for j in range(L):
                if j!=i:d[j]+=s*quad[min(i,j),max(i,j)]
        if h<=cap:yield h,q
        old=q
def adj(a,b,ab):
    L=len(a);out=[0]*(L+2)
    for i in range(L):
        for j in range(i+1,L):out[i]|=1<<j
        out[L if a[i] else i] |= 1<<(i if a[i] else L)
        out[L+1 if b[i] else i] |= 1<<(i if b[i] else L+1)
    out[L if ab else L+1] |= 1<<(L+1 if ab else L)
    return out

def dp(out):
    n=len(out);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            v=d[mask][i]
            if v:
                bits=out[i]&~mask
                while bits:
                    lb=bits&-bits;j=lb.bit_length()-1;bits-=lb
                    d[mask|lb][j]+=v
    return sum(d[-1])
