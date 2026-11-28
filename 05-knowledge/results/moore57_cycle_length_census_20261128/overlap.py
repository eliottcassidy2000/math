"""Independent check of the two-root overlap in the Hoffman–Singleton graph."""
from count import graph
from itertools import combinations, product
from math import prod

def rooted(A,r,m):
    A0=sorted(A[r])
    F={b:tuple(sorted(A[b]-{r})) for b in A0}
    P={(b,c):{x:next(iter(A[x]&set(F[c]))) for x in F[b]}
       for b in A0 for c in A0 if b!=c}
    def canon(V):
        i=V.index(min(V)); V=V[i:]+V[:i]
        return min(V,V[:1]+V[:0:-1])
    out=set()
    for S in combinations(A0,m):
        # m=3, one triangle order up to reversal.
        assert m==3
        for D in product(*(F[b] for b in S)):
            U=[P[S[i-1],S[i]][D[i-1]] for i in range(m)]
            if any(a==b for a,b in zip(U,D)):continue
            V=tuple(x for i in range(m) for x in (U[i],S[i],D[i]))
            out.add(canon(V))
    return out

if __name__=='__main__':
    A=graph(7); r=0; x=min(A[r]); y=min(set(range(len(A)))-A[r]-{r})
    S=rooted(A,r,3)
    assert not (S & rooted(A,x,3))
    pair=S & rooted(A,y,3)
    assert len(pair)==prod(range(4,7))//3
    t,=A[r]&A[y]
    a=sorted(A[r]-{t})
    b=[next(iter(A[u]&(A[y]-{t}))) for u in a]
    W={next(iter(A[b[i]]&A[a[j]])) for i in range(6) for j in range(6) if i!=j}
    assert len(W)==30 and not W.intersection(a+b+[r,y,t])
    assert set(range(50))-W-set(a+b+[r,y,t])==A[t]-{r,y}
    print('Hoffman–Singleton 9-cycle overlap: adjacent 0, nonadjacent',len(pair))
    C=55**21+1
    from fractions import Fraction as F
    singles=3250*(prod(range(37,58))//42)*(C-F(56,37))
    overlaps=(3250*57*56//2)*(prod(range(36,57))//21)
    bnd=singles-overlaps
    print('Degree 57 overlap-corrected branch lower bound',-(-bnd.numerator//bnd.denominator))
