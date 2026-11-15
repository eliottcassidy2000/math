#!/usr/bin/env python3
"""Literal subset path and ordered-run controls for the observable-lattice proof."""
from itertools import combinations, permutations, product
from functools import lru_cache
from math import comb, factorial, gcd
from random import Random
from functools import reduce

def masks(n, bits):
    A=[0]*n
    for t,(i,j) in enumerate(combinations(range(n),2)):
        if bits>>t&1: A[i]|=1<<j
        else:A[j]|=1<<i
    return tuple(A)

def literal_cover(A):
    n=len(A); cs=[0]*n
    for w in permutations(range(n)):
        missing=sum((A[i]>>j&1)==0 for i,j in zip(w,w[1:]))
        for r in range(missing+1,n+1):
            cs[r-1]+=comb(n-1-missing,r-1-missing)
    return cs

def probe(r):
    m=r-1;n=2*m+1; A=[0]*n
    for i in range(1,m+1):
        A[0]|=1<<i; A[i+m]|=1; A[i]|=1<<(i+m)
    for i,j in combinations(range(1,m+1),2):
        A[i]|=1<<j; A[i+m]|=1<<(j+m)
    for i,j in product(range(1,m+1),repeat=2):
        if i!=j:A[j+m]|=1<<i
    return tuple(A)

def a(A,z,r):
    n=len(A); S=((1<<n)-1)^(1<<z)
    @lru_cache(None)
    def dfs(last,remain,more):
        if not remain and not more:return 1
        valid=remain if last==-1 else remain&A[last]
        total=0
        while valid:
            bit=valid&-valid;valid^=bit
            total+=dfs(bit.bit_length()-1,remain^bit,more)
        if more and (last<0 or (A[last]>>z&1)):
            total+=dfs(z,remain,more-1)
        return total
    return dfs(-1,S,r)

def join(S,blocks):
    vertices=[(i,j) for i,B in enumerate(blocks) for j in range(len(B))]
    return tuple(sum(1<<v for v,(k,l) in enumerate(vertices)
                      if ((blocks[i][j]>>l)&1 if i==k else (S[i]>>k)&1))
                 for i,j in vertices)

def filled(S,z,U):
    return join(S,[U if i==z else (0,) for i in range(len(S))])

def hp(A):
    n=len(A); N=1<<n
    D=[[0]*n for _ in range(N)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,N):
        for u in range(n):
            k=(N-1-s)&A[u]
            while k and D[s][u]:
                t=k&-k;k^=t
                D[s|t][t.bit_length()-1]+=D[s][u]
    return sum(D[-1])

def aut(A):
    n=len(A)
    return sum(all((A[i]>>j&1)==(A[v[i]]>>v[j]&1) for i,j in combinations(range(n),2))
               for v in permutations(range(n)))

def main():
    tests=0
    for r in range(1,7):
        A=probe(r)
        assert a(A,0,r)==factorial(r-1) and a(A,0,r+1)==0
        assert all(a(A,0,s)%factorial(s-1)==0 for s in range(1,r+1))
    NT=(6,12,8,1)
    # Natural rows: 0->1,2; 3->0; all other arcs forward.
    assert (a(NT,1,1),a(NT,1,2),a(NT,1,3))==(5,2,0)
    U=(16,1,3,7,14)
    C=literal_cover(U)
    assert C==[9,66,186,240,120] and aut(U)==1
    assert reduce(gcd,[factorial(i)*x for i,x in enumerate(C)])==3
    assert hp(filled(probe(2),0,NT))==37
    rng=Random(25511)
    bank=[masks(n,x) for n in range(2,5) for x in range(1<<(n*(n-1)//2))]
    for _ in range(200):
        A=rng.choice(bank); z=rng.randrange(len(A)); B=rng.choice(bank)
        c=literal_cover(B)
        assert hp(filled(A,z,B))==sum(a(A,z,r)*x for r,x in enumerate(c,1))
        assert all(a(A,z,r)%factorial(r-1)==0 for r in range(1,len(B)+1))
        if c[0]%3==0:assert hp(filled(A,z,B))%3==0
    for _ in range(60):
        S=masks(3,rng.randrange(8)); blocks=[rng.choice(bank) for i in range(3)]
        vals=[hp(join(S,p)) for p in itertools_permute(blocks)]
        target=(hp(S)*prod_of(hp(B) for B in blocks)+2*sum(a(S,i,2) for i in range(3)))%4
        assert all(v%4==target for v in vals)
        tests+=len(vals)
    print({'probe_levels':6,'random_literal_inflations':200,'module_permutation_tests':tests,'rigid_certificate':C})
def itertools_permute(x):return permutations(x)
def prod_of(x):
    p=1
    for t in x:p*=t
    return p
if __name__=='__main__':main()
