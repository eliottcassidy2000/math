from itertools import product
from collections import defaultdict
def path_arrays(A):
    n=len(A); F=(1<<n)-1
    out=[[0]*n for _ in range(F+1)]; h=[0]*(F+1); h[0]=1
    for i in range(n):out[1<<i][i]=1
    for S in range(1,F+1):
        for i,x in enumerate(out[S]):
            if x:
                B=A[i]&~S
                while B:
                    b=B&-B; B^=b
                    out[S|b][b.bit_length()-1]+=x
        h[S]=sum(out[S])
    return h,out
def tensors(A):
    n=len(A);F=(1<<n)-1
    h,E=path_arrays(A);tot=h[F];f=sum(x*y for x,y in zip(h,h[::-1]))
    _,S=path_arrays([r^F^(1<<i) for i,r in enumerate(A)])
    p=[sum(a[i]*b for a,b in zip(S,h[::-1])) for i in range(n)]
    q=[sum(a[i]*b for a,b in zip(E,h[::-1])) for i in range(n)]
    for U in range(F+1):
        P=sum(p[i] for i in range(n) if U>>i&1)
        Q=sum(q[i] for i in range(n) if not U>>i&1)
        Ss=sum(E[F][i] for i in range(n) if U>>i&1)
        T=sum(S[F][i] for i in range(n) if not U>>i&1)
        K=path_arrays(tuple(r if U>>i&1 else r|(1<<n) for i,r in enumerate(A))+(U,))[0][-1]
        yield (tot,f,P,Q,Ss,T,K,U)
def canonical(r):
    n=len(r);ds=tuple(i.bit_count() for i in r);classes=defaultdict(list)
    for i,d in enumerate(ds):classes[d].append(i)
    from itertools import permutations
    best=None
    for axes in product(*(permutations(classes[d]) for d in sorted(classes))):
        ix=sum(axes,())
        perm=tuple(sum(((r[i]>>j)&1)<<t for t,j in enumerate(ix)) for i in ix)
        if best is None or perm<best:best=perm
    return best
def catalogue(m):
    layer={(0,)}
    for n in range(1,m+1):
        yield n,layer
        if n==m:break
        layer={canonical(tuple(r if U>>i&1 else r|(1<<n) for i,r in enumerate(A))+(U,)) for A in layer for U in range(1<<n)}
def longest(T):
    starts=(x for x in T if x-2 not in T)
    best=[]
    for j in starts:
        L=[]
        while j in T:L.append(j);j+=2
        if len(L)>len(best):best=L
    return best
if __name__=='__main__':
    for n,layer in catalogue(7):
        for rows in sorted(layer):
            print(n,*rows)
