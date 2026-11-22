from itertools import combinations,permutations
from collections import Counter
edges=list(combinations(range(4),2))
orders=list(permutations(range(4)))
def rows(m):
    a=[0]*4
    for k,(i,j) in enumerate(edges):
        if m>>k&1:a[j]|=1<<i
        else:a[i]|=1<<j
    return a
def code(a,v):
    return sum((1-((a[v[i]]>>v[j])&1))<<k for k,(i,j) in enumerate(edges))
def canonical(m):
    return min(code(rows(m),v) for v in orders)
def h(m):
    a=rows(m)
    return sum(all(a[i]>>j&1 for i,j in zip(v,v[1:])) for v in orders)
assert [canonical(m) for m in (0,2,16,18)]==[0,2,5,4]
assert [h(m) for m in (0,2,16,18)]==[1,3,3,5]
R={0:{0:3,2:1,4:1,5:1},2:{0:3,4:3},4:{4:3,5:1,0:1,2:1},5:{0:3,4:3}}
assert all(Counter(canonical(m^(1<<i)) for i in range(6))==R[canonical(m)] for m in range(64))
pg=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
def step(n):return next((a,b) for a,b in pg if n%b==0)
N=2**60*(2**61-1)-1
assert N==(2**60-1)*(2**61+1)
assert all(N%p==0 and N%(p*p) for p in (7,11,13))
assert [step(m) for m in (11,13,N,N+2)]==[(13,11),(11,13),(17,91),(55,1)]
print('four-vertex flip types and PRIMEGAME priority on E61 neighbors pass')
