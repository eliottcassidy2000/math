"""Standalone controls for the labeled 3-to-5 split and the geometry count."""
from collections import Counter
from math import comb, factorial, isqrt

def children(m):
    assert m%2==1 and m%3
    es=[j for j in range(1,7) if pow(2,j,3)*m%3==1]
    assert sorted(pow(2,j,9)*m%9 for j in es)==[1,4,7]
    c=tuple((2**j*m-1)//3 for j in es)
    assert all(u>0 and u&1 for u in c)
    assert sum(u%3==0 for u in c)==1
    return c

for m in range(1,1000,2):
    if m%3:
        children(m)

S={5}
for k in range(1,13):
    S={u for m in S for u in children(m) if u%3}
    assert len(S)==2**k

# Distinct histories may lead to the same marked/labeled frontier.
result=Counter({frozenset({5}):1})
for h in range(8):
    assert len(result)==comb(2*h,h)//(h+1)
    assert sum(result.values())==factorial(h)
    assert all(len(F)==2*h+1 for F in result)
    assert all(sum(u%3==0 for u in F)==h for F in result)
    new=Counter()
    for F, ways in result.items():
        for m in F:
            if m%3:
                T=frozenset((F-{m}).union(children(m)))
                assert len(T)==len(F)+2
                new[T]+=ways
    result=new

for d in range(1,2001):
    marked=[(a+2,d//a+2) for a in range(1,d+1) if d%a==0]
    fix=sum(p==q for p,q in marked)
    assert fix==int(isqrt(d)**2==d)
    assert len({tuple(sorted(z)) for z in marked})==(len(marked)+fix)//2

forms=[(2*a,2*b,2*c) for a in range(2,20)
       for b in range(a+1,20) for c in range(b+1,20) if a*b+a*c+b*c>=a*b*c]
assert forms==[(4,6,8),(4,6,10),(4,6,12)]
print('PASS: 12 binary levels, factorial/Catalan histories to 7 expansions, 2000 Schläfli shells, even nonnegative figures')
