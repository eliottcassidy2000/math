"""Exact source-preserving fresh-block migration matrices. Standard library only."""
from functools import lru_cache
from itertools import permutations

def tournament(n, code):
    out=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            u,v=(i,j) if code>>k&1 else (j,i)
            out[u]|=1<<v;k+=1
    return tuple(out)

def multiply(M,v):return tuple(sum(a*b for a,b in zip(row,v)) for row in M)

def compiler(local, z):
    """Old port is 0; z is a distinct new root. All fresh vertices lose to A.
    Input/output state order (h,F,q,g). Local is an actual tournament.
    """
    n=len(local)
    if n<2 or not 0<z<n:raise ValueError('need old port 0 and a distinct new root')
    for i in range(n):
        if local[i]>>n or local[i]>>i&1:raise ValueError('invalid tournament')
        for j in range(i):
            if ((local[i]>>j&1)+(local[j]>>i&1))!=1:raise ValueError('invalid tournament')
    full=(1<<n)-1; fresh=full^1;B=fresh^(1<<z);R=local[z]& (B|1)
    @lru_cache(None)
    def paths(mask):
        if not mask:return ((),)
        return tuple(p for p in permutations(i for i in range(n) if mask>>i&1)
                     if all(local[a]>>b&1 for a,b in zip(p,p[1:])))
    def subsets(mask):
        s=mask
        while True:
            yield s
            if not s:break
            s=(s-1)&mask
    def H(mask):return len(paths(mask))
    def hp_start(mask, cut):return sum(bool(p) and bool(cut>>p[0]&1) for p in paths(mask))
    def hp_end(mask, cut):return sum(bool(p) and bool(cut>>p[-1]&1) for p in paths(mask))
    @lru_cache(None)
    def Frow(U):
        if not U:return (0,1,0,0)
        h=H(U); s=hp_start(U,local[0]);p=0
        for P in subsets(U):
            if P:
                p+=hp_end(P,U^local[0]) *H(U^P)
        return (H(U|1)-h,h,p,h-s)
    hrow=Frow(B);frow=Frow(fresh)
    if not B:
        grow=(1,0,0,-1) if R&1 else (0,0,0,0)
    else:
        e=hp_end(B,R)
        se=sum(bool(local[0]>>p[0]&1) and bool(R>>p[-1]&1) for p in paths(B))
        pe=0
        for P in subsets(B):
            Q=B^P
            if P and Q:pe+=hp_end(P,B^local[0])*hp_end(Q,R)
        grow=(hp_end(B|1,R)-e,e,pe,e-se)
    qrow=[0,0,0,0]
    for U in subsets(B):
        V=B^U
        s=hp_start(V,R)
        for i,a in enumerate(Frow(U)):qrow[i]+=s*a
        a=0
        for P in subsets(V|1):
            if P&1:
                prefix=sum(p[-1]==0 and bool(R>>p[0]&1) for p in paths(P))
                a+=prefix*H((V|1)^P)
        qrow[2]+=H(U)*a
        qrow[0]+=H(U)*hp_start(V|1,R)
    return tuple(map(tuple,(hrow,frow,qrow,grow)))

def extend(A,cut,local,z):
    """Return new unmarked parent/cut; labels retain A then local 0..m-1."""
    n=len(A);m=len(local);out=list(A)+[x<<n for x in local]
    for i in range(n):
        if cut>>i&1:out[n]|=1<<i
        else:out[i]|=1<<n
        for j in range(1,m):out[i]|=1<<(n+j)
    apex=n+z;keep=[i for i in range(n+m) if i!=apex]
    parent=tuple(sum(1<<j for j,v in enumerate(keep) if out[u]>>v&1) for u in keep)
    newcut=sum(1<<j for j,v in enumerate(keep) if out[apex]>>v&1)
    return parent,newcut

if __name__=='__main__':
    for n in range(2,6):
        bank={compiler(tournament(n,c),n-1) for c in range(1<<(n*(n-1)//2))}
        print(n,len(bank))
