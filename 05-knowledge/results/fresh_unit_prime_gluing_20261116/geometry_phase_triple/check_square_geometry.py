from itertools import combinations
from math import isqrt
from fractions import Fraction

def graph(n):
    D={v:set() for v in range(1,n+1)}
    for u,v in combinations(D,2):
        if isqrt(u+v)**2==u+v:
            D[u].add(v);D[v].add(u)
    return D
def cycle_lengths_at_n(n):
    D=graph(n)
    def walk(P,A):
        if n in D[P[-1]] and len(P)>=3:
            yield (len(P),5 in A,8 in A,9 in A)
        for q in D[P[-1]]-A:
            yield from walk(P+[q],A|{q})
    return set(walk([n],{n}))

first_lengths={}
first_joint=None
for n in range(1,23):
    for k,_,has8,has9 in cycle_lengths_at_n(n):
        first_lengths.setdefault(k,n)
        if has8 and has9 and first_joint is None:first_joint=n
assert [first_lengths[k] for k in range(5,10)]==[15,19,19,20,22]
assert first_joint==22
assert not any(k in first_lengths for k in (3,4))

D=graph(19)
while any(len(a)<2 for a in D.values()):
    v=next(v for v in D if len(D[v])<2)
    for w in D.pop(v):D[w].remove(v)
paths=[(1,3,6),(1,15,10,6),(1,8,17,19,6)]
assert {frozenset((u,v)) for u in D for v in D[u]}=={frozenset(e) for p in paths for e in zip(p,p[1:])}
chi=Fraction(1,5)+Fraction(1,6)+Fraction(1,7)-Fraction(1,2)
assert chi==Fraction(1,105)
print('Q19 theta face lengths 5,6,7; junction curvature/(2 pi) =',chi)
print('lengths 5 through 9 first appear at', [first_lengths[k] for k in range(5,10)])
print('first common cycle through 8 and 9 at Q',first_joint)
