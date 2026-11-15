#!/usr/bin/env python3
"""Exact top-shelf lines using only strong profiles through order 7.

Input: n H P2 canonical-code strong-code TSV, complete at least through 7.
The normalized threshold 4 is strict.
"""
import sys
from fractions import Fraction

atoms=set()
for line in sys.stdin:
    n,h,p,code,sc=line.split()
    n,h,p = map(int, (n,h,p))
    q=Fraction(p+2*h,2**n)
    if 3<=n<=7 and sc!="-" and q<4:
        atoms.add((n,h,q))

# Below 4, triangle cost is at most 5; hence every strong component
# has at most 7 vertices. Enumerate all ordered strings of eligible atoms.
def extend(n,h,q):
    yield n,h,q
    for m,g,r in atoms:
        if q*r<4:
            yield from extend(n+m,h*g,q*r)

shelf=sorted(set(extend(0,1,Fraction(1))),key=lambda s: (s[2],s[1],s[0]))
print("retained vertices, H, normalized Q")
for s,h,q in shelf:
    print(s,h,q)

def targets(n,U,V,cap):
    W=U-2*V
    lower=V+min(Fraction(0),Fraction(W,n+1))
    assert lower>0 and cap < 4*2**n*lower, "cap outside certified shelf"
    return sorted({int(N) for s,h,q in shelf if s<=n and (N:=V*2**n*q+W*h)<=cap})

for n,U,V,want in [(12,11,6,[24575,36861,49147,55287]),(13,11,6,[49151]),
                   (13,5,3,[24575,36861,49147,55287]),(14,5,3,[49151])]:
    got=targets(n,U,V,59999)
    assert got==want,(n,U,V,got)
    print(n,U,V,got)
