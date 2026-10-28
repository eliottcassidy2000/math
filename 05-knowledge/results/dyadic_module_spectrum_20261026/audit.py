#!/usr/bin/env python3
"""Independent finite baseline enumeration and exact F4 chart controls."""
from functools import lru_cache
from itertools import product
from fractions import Fraction
import json
import check as main

def require(p,m):
    if not p:raise RuntimeError(m)

def baseline_audit():
    n=5;pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    positions={p:k for k,p in enumerate(pairs)}
    def orientation(A,B):
        return sum(1<<positions[i,j] for i in range(n) for j in range(i+1,n) if A>>i&1 and B>>j&1)
    @lru_cache(None)
    def build(S):
        if S.bit_count()<=1:return frozenset([0])
        result=set()
        # All J expressions; this enumeration stores actual graphs, not costs.
        for A in range(1,S):
            if A&S!=A:continue
            B=S^A
            X=orientation(A,B)
            for x in build(A):
                for y in build(B):result.add(x|y|X)
        # All ordered three-part expressions; duplications disappear in the set.
        for A in range(1,S):
            if A&S!=A:continue
            R=S^A
            for B in range(1,R):
                if B&R!=B:continue
                C=R^B;X=orientation(A,B)|orientation(B,C)|orientation(C,A)
                for x in build(A):
                    for y in build(B):
                        for z in build(C):result.add(x|y|z|X)
        return frozenset(result)
    bank=build((1<<n)-1);hist={}
    for code in range(1<<len(pairs)):
        T=[0]*n
        for k,(i,j) in enumerate(pairs):
            if code>>k&1:T[i]|=1<<j
            else:T[j]|=1<<i
        exact=min((code^b).bit_count() for b in bank)
        dp,_=main.parse(T)
        require(dp==exact,'independent baseline bank versus subset-cost recursion')
        hist[exact]=hist.get(exact,0)+1
    return dict(vertices=n,labelled_targets=1<<len(pairs),distinct_baselines=len(bank),distance_histogram=hist)

def chart_audit():
    # Address field and coefficient ring are distinct copies with an explicit identification.
    alpha=(3,1)
    perm=[]
    for i in range(4):
        a,b=main.mul(alpha,(i&1,i>>1));perm.append((a%2)+2*(b%2))
    require(perm==[0,3,1,2],'F4 address permutation')
    pairings=[frozenset([frozenset((0,i)),frozenset(j for j in range(4) if j not in (0,i))]) for i in (1,2,3)]
    for P in pairings:
        require(frozenset(frozenset(perm[x] for x in pair) for pair in P) in pairings,'three pairings permuted')
    V=[(0,0),(1,0),(0,1),(1,2),(-2,-1)]
    checks=0
    for z in product(V,repeat=4):
        w=[None]*4
        for i in range(4):w[perm[i]]=main.mul(alpha,z[i])
        require(sum(main.Q(v) for v in w)==7*sum(main.Q(v) for v in z),'source norm normalization')
        left=tuple(map(sum,zip(*z)));right=tuple(map(sum,zip(*w)))
        require(main.Q(right)==7*main.Q(left),'fixed terminal norm covariance')
        require(5*right[0]-right[1]==7*(2*left[0]-left[1]),'marked source row covariance')
        for P in pairings:
            v=sum(main.Q(tuple(sum(z[i][j] for i in pair) for j in (0,1))) for pair in P)
            transformed=sum(main.Q(tuple(sum(w[perm[i]][j] for i in pair) for j in (0,1))) for pair in P)
            require(transformed==7*v,'whole chart family intertwines')
        checks+=1
    # A hostile with identical allowed source norms and no cancellation:
    # changing the parse cannot erase its endpoint difference.
    aligned=[]
    for K in range(1,13):
        M=2**K;E=[M*M//(2**k) for k in range(K+1)]
        B=main.budgets(E)
        require(all(B[key]==M*M-M for key in ('binary','fourary','eightary','adaptive_1to3')),'real aligned debt survives every controller')
        aligned.append([M,B['adaptive_1to3']])
    centred=0
    beta=(1,2)
    def embed(z):a,b=z;return (a,-b,b-a)
    def skew(v):x,y,z=v;return (y-z,z-x,x-y)
    require(main.mul(beta,beta)==(-3,0),'C3 spectral field')
    for a,b in product(range(-7,8),repeat=2):
        z=(a,b);v=embed(z)
        require(sum(t*t for t in v)==2*main.Q(z),'A2 centred lattice metric')
        require(skew(v)==embed(main.mul(beta,z)),'C3 skew intertwiner')
        rows=(2*a-b,-a-b,-a+2*b)
        require(rows==embed(main.mul((2,1),z)),'actual marked ROOT rows')
        require(sum(t*t for t in rows)==6*main.Q(z),'actual three-row energy')
        centred+=1
    towers=[]
    for k in range(1,9):
        n=3**k
        # Positive singular magnitude 3^j sqrt(3), twice per eigenpair.
        multiplicities=[(3**j,2*3**(k-j-1)) for j in range(k)]
        require(sum(c for m,c in multiplicities)==n-1,'C3 tower rank')
        require(sum(3*m*m*c for m,c in multiplicities)==n*(n-1),'C3 tower Frobenius energy')
        require(sum(m*c for m,c in multiplicities)==2*k*3**(k-1),'C3 tower trace norm coefficient')
        towers.append([k,n,sum(m*c for m,c in multiplicities)])
    return dict(address_permutation=perm,intertwining_cases=checks,aligned_hostile=aligned,
                centred_C3_controls=centred,C3_tower_trace_over_sqrt3=towers)

if __name__=='__main__':print(json.dumps(dict(status='PASS',baseline=baseline_audit(),charts=chart_audit()),indent=2,sort_keys=True))
