#!/usr/bin/env python3
"""Exact local orbit-quotient checks (no hypothetical graph instantiated)."""
from itertools import combinations, product

def partitions(n, upper):
    if n == 0:
        yield ()
    for j in range(1, min(n, upper)+1):
        for p in partitions(n-j, j):
            yield (j,)+p

def pattern(p, f, diagonal):
    row_degree=57-f-diagonal
    square_excess=56+p-f-diagonal*diagonal-diagonal
    # Everything other than doubles/triples is bounded by available excess.
    out=set()
    for multi in partitions(square_excess-row_degree, square_excess-row_degree):
        if any(x%2 for x in multi):
            continue
        vals=[]
        for x in multi:
            z=2
            while z*(z-1)<x:
                z+=1
            if z*(z-1)!=x:
                break
            vals.append(z)
        else:
            rest=row_degree-sum(vals)
            if rest>=0:
                out.add(tuple(sorted(vals,reverse=True))+(1,)*rest)
    return out

assert pattern(3,10,0)=={(2,)+(1,)*45}
assert pattern(5,50,0)=={(2,2)+(1,)*3}
assert pattern(5,50,2)=={(1,)*5}
assert pattern(11,5,0)=={(2,)*5+(1,)*42,(3,2,2)+(1,)*45}
assert pattern(11,5,2)=={(2,2,2)+(1,)*44,(3,)+(1,)*47}
for A in product(list(combinations(range(5),2)),repeat=4):
    assert {sum(a)%5 for a in product(*A)}==set(range(5))
assert [v for v in range(141) if v%2==0 and (2*v+8*(690-10))%15==0]==[10,40,70,100,130]
print('Proved local row-pattern checks for orders 3, 5, 11; a doubled 4-cycle always lifts modulo 5')
