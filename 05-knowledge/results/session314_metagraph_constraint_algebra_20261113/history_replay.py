from itertools import permutations, combinations
from math import isqrt, gcd

def orderbits(n,m):
    return {(i,j): ((m>>r)&1) for r,(i,j) in enumerate(combinations(range(n),2))}
def arc(A,i,j):
    return A[i,j] if i<j else 1-A[j,i]
for n in range(2,6):
    for bits in range(1<<((n-1)*(n-2)//2)):
        core={(i,i+1):1 for i in range(n-1)}
        for r,(i,j) in enumerate((i,j) for i,j in combinations(range(n),2) if j>i+1):
            core[i,j]=(bits>>r)&1
        for k in [2,3]:
            N=n+2*k
            W={(i,j): (1 if j==i+1 else (1-core[i-k,j-k] if k<=i<j<k+n else 0)) for i,j in combinations(range(N),2)}
            scores=[sum(arc(W,i,j) for j in range(N) if j!=i) for i in range(N)]
            assert max(scores[:k])<k and k<=min(scores[k:k+n])<=max(scores[k:k+n])<k+n and min(scores[k+n:])>=k+n
            assert [(i,j) for i in range(k) for j in range(k,k+n) if arc(W,i,j)]==[(k-1,k)]
            assert [(i,j) for i in range(k,k+n) for j in range(k+n,N) if arc(W,i,j)]==[(k+n-1,k+n)]

def compose(a,b): return tuple(a[b[i]] for i in range(5))
def inv(g): return tuple(g.index(i) for i in range(5))
id5=tuple(range(5))
def group_single(g):
    D={id5}; h=g
    while h not in D:
        D.add(h);h=compose(h,g)
    return frozenset(D)
G=[g for g in permutations(range(5)) if sum(g[i]>g[j] for i,j in combinations(range(5),2))%2==0]
C3={group_single(g) for g in G if len(group_single(g))==3}
C5=tuple({group_single(g) for g in G if len(group_single(g))==5})
def caction(g,i): return C5.index(frozenset(compose(compose(g,h),inv(g)) for h in C5[i]))
assert len(C3)==10 and len(C5)==6
pairs={}
for C in C3:
    g=next(h for h in C if h!=id5)
    a=frozenset(caction(h,0) for h in C)
    pairs[a]=frozenset(i for i in range(5) if g[i]==i)
assert len(pairs)==len(set(pairs.values()))==10
def edges(a):return {frozenset([i,j]) for i in a for j in range(6) if j not in a}
for a,b in combinations(pairs,2):
    assert len(edges(a)&edges(b))==5
def marked(k):
    F=[a if k in a else set(range(6))-a for a in pairs]
    return frozenset(frozenset([i,j]) for i,j in combinations(range(10),2) if len(F[i]&F[j])==1)
assert len({marked(k) for k in range(6)})==6
for k in range(6):assert all(sum(i in e for e in marked(k))==3 for i in range(10))
D={i:{j for j in range(1,20) if i!=j and isqrt(i+j)**2==i+j} for i in range(1,20)}
assert {i for i in D if len(D[i])<=1}=={16,18}
assert D[16]=={9} and D[18]=={7} and D[9]=={7,16}
for k in range(1,6):
    mod=19**k
    assert pow(2048,18*19**(k-1),mod)==1
    for p in [2,3]+([19] if k>=2 else []):assert pow(2048,18*19**(k-1)//p,mod)!=1
print("Endpoint geometry, A5 frames, Q19 obstruction and full-period padding controls passed")

