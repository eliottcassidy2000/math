from functools import cache
from collections import defaultdict
from itertools import permutations

@cache
def shape(n,d):
    """Integer arithmetic, uncolored and unordered, with 1's edge cut."""
    if not d:return ()
    children=[2*n]
    if n%3==2 and n!=2:
        children.append((2*n-1)//3)
    return tuple(sorted(shape(k,d-1) for k in children))

def unit_controls():
    for r in range(1,5):
        M=3**r
        a=[n for n in range(3,4*M+3) if n%3]
        quotient=defaultdict(set)
        for n in a:quotient[shape(n,4*r-2)].add(n%M)
        assert all(len(residues)==1 for residues in quotient.values())
        # The bound is not sharp in general. Measure earliest successful
        # depth on this finite sample, without promoting it to a theorem.
        for d in range(4*r-1):
            q=defaultdict(set)
            for n in a:q[shape(n,d)].add(n%M)
            if all(len(S)==1 for S in q.values()):
                print('r',r,'modulus',M,'finite minimum distinguishing depth',d)
                break
    assert shape(1,4)!=shape(2,4)
    assert all(shape(a,4)!=shape(b,4) for a in (1,2) for b in range(3,200))

def children(n):
    return [2*n]+([(2*n-1)//3] if n%3==2 and n!=2 else [])

def growth_controls():
    A={1}
    for h in range(14):
        frontier={y for x in A for y in children(x)}-A
        assert len(frontier)==1+sum(len(children(x))==2 for x in A)
        A|=frontier
    A={1,2,4,8,16,5,32}
    parents={n:p for p in A for n in children(p) if n in A}
    schedules=0
    for s in permutations(A):
        rank={x:i for i,x in enumerate(s)}
        schedules+=all(rank[p]<rank[n] for n,p in parents.items())
    from math import factorial,prod
    def size(v):return 1+sum(size(w) for w in children(v) if w in A)
    assert schedules*prod(map(size,A))==factorial(len(A))
    print('example tree schedules',schedules)

if __name__=='__main__':
    unit_controls()
    growth_controls()
