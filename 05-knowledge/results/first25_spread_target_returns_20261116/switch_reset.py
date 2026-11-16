from itertools import combinations,permutations
V=tuple(range(5))
E=tuple(combinations(V,2))
K=(1<<len(E))-1
def mask(edges): return sum(1<<E.index(tuple(sorted(e))) for e in edges)
def cut(D):return mask((a,b) for a,b in E if (a in D)!=(b in D))
def euler(x):
    d=set()
    for i,e in enumerate(E):
        if x>>i&1:d.symmetric_difference_update(e)
    return x^cut(d)
moves={euler(K^mask(zip(p,p[1:]))) for p in permutations(V)}
space={euler(x) for x in range(K+1) if x.bit_count()%2==0}
def q(x):return (x.bit_count()//2)%2
def b(x,y):return (x&y).bit_count()%2
P={i:{x for x in moves if not any(x>>j&1 for j,e in enumerate(E) if i in e)} for i in V}
assert len(moves)==15 and all(len(o)==3 for o in P.values())
for o in P.values():
    assert all(x^y in o and not b(x,y) for x,y in combinations(o,2))
def rank(x):return 0 if not x else 3 if x==K else 2 if q(x) else 1
for x in space:
    d=x
    while d:
        t=next(t for t in sorted(moves) if rank(d^t)==rank(d)-1)
        d^=t
    assert d==0
layers={0};ball={0}
for i in range(1,4):
    layers={x^y for x in layers for y in moves}-ball
    ball|=layers
    assert layers=={x for x in space if rank(x)==i}
print('reset layers:',*[sum(rank(x)==i for x in space) for i in range(4)])
