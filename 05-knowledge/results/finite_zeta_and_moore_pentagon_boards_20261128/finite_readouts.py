"""Finite controls for the reciprocal arc readouts; not evidence of infinitude."""
from math import comb, isqrt

def prime(p):
    return all(p % d for d in range(2, isqrt(p)+1))

def bernoulli(p):
    B=[1]
    for m in range(1,p-2):
        B.append(-sum(comb(m+1,j)*B[j] for j in range(m))*pow(m+1,-1,p)%p)
    return B

def readouts(p,order):
    out=[]
    for k in [3,5]:
        pref=z=0
        for x in order:
            z+=pref*pow(x,1-k,p)
            pref+=pow(x,-1,p)
        out.append(z%p)
    return out

seen=0
for p in range(7,251):
    if not prime(p):continue
    B=bernoulli(p)
    z=readouts(p,list(range(1,p)))
    assert z==[B[p-k] for k in [3,5]]
    assert readouts(p,list(range(p-1,0,-1)))==[(-x)%p for x in z]
    b=readouts(p,[2*i%p for i in range(1,p)])
    assert b==[pow(2,-k,p)*x%p for k,x in zip([3,5],z)]
    a0,a1=sum(z)%p,sum(b)%p
    y=(a1-pow(2,-3,p)*a0)*pow((pow(2,-5,p)-pow(2,-3,p))%p,-1,p)%p
    assert [(a0-y)%p,y]==z
    if p in [7,11,13]:print(p,z,a0,a1)
    seen+=1
print("Prime controls:",seen)

from itertools import combinations, permutations
def tournament(backwards):
    V=range(6)
    def arc(x,y):
        return (x<y) != (tuple(sorted((x,y))) in backwards)
    def cycles(k):
        return [S for S in combinations(V,k) for L in permutations(S[1:])
                if all(arc(a,b) for a,b in zip((S[0],*L),(*L,S[0])))]
    tri=cycles(3)
    return (len(tri),len(cycles(4)),len(cycles(5)),
            sum(set(a).isdisjoint(b) for a,b in combinations(tri,2)),
            sum(all(arc(a,b) for a,b in zip(W,W[1:])) for W in permutations(V)))
base={(0,1),(0,5),(1,4)}
assert tournament(base)==(6,7,6,2,33)
assert tournament(base|{(1,2)})==(6,7,6,1,29)
print("Six-vertex single flip:",tournament(base),tournament(base|{(1,2)}))
