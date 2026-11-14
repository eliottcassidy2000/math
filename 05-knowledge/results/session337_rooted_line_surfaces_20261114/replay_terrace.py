from terrace import canonical, pathcount, four_state, pair_inverse
from itertools import product

def make(g):
    n=sum(g)+len(g)+1
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    a=0
    for r in g:
        b=a+1+r
        rows[a]^=1<<b;rows[b]^=1<<a;a=b
    return rows

def literal(M):
    n=len(M);v=[[0]*n for _ in range(1<<n)]
    for i in range(n):v[1<<i][i]=1
    for s in range(1,1<<n):
        for i,k in enumerate(v[s]):
            if not k:continue
            z=M[i]&~s
            while z:
                b=z&-z;z-=b
                v[s|b][b.bit_length()-1]+=k
    return sum(v[-1])

from itertools import permutations, groupby
def canon(M):
    n=len(M);order=sorted(range(n),key=lambda i:M[i].bit_count())
    return min(tuple(sum((M[i]>>j&1)<<k for k,j in enumerate(q)) for i in q)
        for p in product(*(tuple(permutations(g)) for _,g in groupby(order,key=lambda i:M[i].bit_count())))
        for q in [sum(p,())])
for n in range(4,10):
    C=set()
    for g in canonical(n):
        T=make(g);deg=sorted(x.bit_count() for x in T)
        assert deg==[1]+list(range(1,n-1))+[n-2]
        assert pathcount(g)==literal(T)==four_state(g)
        C.add(canon(T))
    assert len(C)==2**(n-4)
    print('order',n,'exact score classes',len(C))
for k in range(1,7):
    for R in product((1,2,3),repeat=k):
        val=1
        for a in R:val*=1+2**a
        assert pathcount(R)==val
for a,b in product(range(1,16),repeat=2):
    N=(1+2**a)*(1+2**b)
    assert pair_inverse(N)==tuple(sorted((a,b)))
    assert pathcount((a,0,b))==2*N-1
    assert pair_inverse(((2*N-1)+1)//2)==tuple(sorted((a,b)))
for n in range(6,26):
    X=2**(n-5)
    assert [pathcount(g) for g in ((n-2,), (2,n-5),(1,0,n-5))]==[1+8*X,5+5*X,5+6*X]
print('Product and shifted-product inverses and nonmonotone chain checked')
assert 11 not in {3}|{pathcount(g) for n in range(4,8) for g in canonical(n)}
# Literal two-hub witness 110, hubs 0,1 with 0->1.
R=[2,0,24,16,0]
for i,d in enumerate((1,1,0),2):
    for h in (0,1):
        if d>>(1-h)&1:R[h]|=1<<i
        else:R[i]|=1<<h
assert literal(R)==11
print('11 is attained but not on any first strong terrace')
