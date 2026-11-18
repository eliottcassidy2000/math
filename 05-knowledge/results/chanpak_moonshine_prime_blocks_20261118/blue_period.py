"""Literal six-vertex metagraph control and bivariate pure-blue checks."""
from itertools import combinations, permutations, product
from functools import lru_cache
from collections import defaultdict,Counter
from math import comb

def make(n,B,z):
    T=[0]*n
    for i in range(n-1):T[i]|=1<<(i+1)
    for k,(i,j) in enumerate(B):
        if z>>k&1:T[j]|=1<<i
        else:T[i]|=1<<j
    return tuple(T)
def rev(A):
    return tuple((1<<len(A))-1^u^(1<<i) for i,u in enumerate(A))
def restrict(A,I):
    return tuple(sum((A[i]>>j&1)<<k for k,j in enumerate(I)) for i in I)
@lru_cache(None)
def iso(A,merged):
    keys=[]
    for B in (A,rev(A)) if merged else (A,):
        D=defaultdict(list)
        for i,x in enumerate(B):D[x.bit_count()].append(i)
        for W in product(*(permutations(D[d]) for d in sorted(D))):
            keys.append(restrict(B,tuple(i for V in W for i in V)))
    return min(keys)
def H(A):
    n=len(A)
    return sum(all(A[i]>>j&1 for i,j in zip(I,I[1:])) for I in permutations(range(n)))
def plain(n):
    B=[(i,j) for i,j in combinations(range(n),2) if j>i+1]
    m=1<<len(B)
    labels=[]
    Cs=defaultdict(set);oriented=set();lines=Counter()
    for z in range(m):
        A=make(n,B,z)
        C=iso(A,True);D=iso(A,False);oriented.add(D)
        isblue=A==restrict(rev(A),tuple(range(n-1,-1,-1)))
        Cs[C].add(isblue)
        labels.append((C,isblue))
    for z,(C,blue) in enumerate(labels):
        v=z^(m-1)
        if z<v:
            D,b=labels[v];assert b==blue
            lines[tuple(sorted((C,D))),blue]+=1
    pure=sorted(H(X) for X,S in Cs.items() if S=={True})
    T=labels[0][0]
    incidents=[(H(e[1] if e[0]==T else e[0]),c,w)
               for (e,c),w in lines.items() if T in e]
    print('order',n,'oriented,merged',len(oriented),len(Cs),'pure-blue heights',pure,'transitive incidences',incidents)
    return pure

def expected(n):
    D=[]
    for t,h in ((0,1),(1,1),(3,3),(5,15)):
        if n<t or (n-t)&1:continue
        m=(n-t)//2
        for j in range(m//3+1):
            i=m-3*j
            D.extend([h*(9**j)]*comb(i+j,j))
    return sorted(D)

if __name__=='__main__':
    for n in range(3,7):assert plain(n)==expected(n)
    for i in range(64):
        for j in range(64):
            assert comb(i+j,i)%2==int(not i&j)
    print('expected normalized order-six class count:', expected(6))
