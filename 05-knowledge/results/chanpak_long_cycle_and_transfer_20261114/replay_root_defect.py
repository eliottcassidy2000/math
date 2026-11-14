from itertools import combinations
from binary_cycle import literal
from functools import cache

def rooted(A,v):
    n=len(A);Q=(1<<n)-1
    L=literal(A)
    @cache
    def start(m,x):
        p=m^(1<<x)
        if not p:return 1
        return sum(start(p,j) for j in range(n) if p>>j&1 and A[x]>>j&1)
    k=sum(start(m,v)*L[Q^m] for m in range(1,1<<n) if m>>v&1)
    return L[-1],k

if __name__=='__main__':
    for n in range(1,6):
        edges=tuple(combinations(range(n),2))
        for k in range(1<<len(edges)):
            A=[0]*n
            for j,(a,b) in enumerate(edges):
                if k>>j&1:A[a]|=1<<b
                else:A[b]|=1<<a
            for v in range(n):
                h,r=rooted(A,v);d=A[v].bit_count()
                if d:assert r-h>=(1<<(d-1))-1
    for d in range(1,8):
        v,z=d,d+1
        T=[sum(1<<x for x in range(j+1,d))|(1<<z) for j in range(d)]
        T += [(1<<d)-1,1<<v]
        h,r=rooted(T,v)
        assert h==(1<<d)+1 and r-h==(1<<d)-2
    print('root defect for all labelled inputs through order 5 and exponential-scale examples')
