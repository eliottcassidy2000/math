from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from math import comb
E=list(combinations(range(4),2))
F=list(combinations(range(4),3))
flags=list(permutations(range(4)))
assert set(Counter(w[0] for w in flags).values())=={6}
assert set(Counter(tuple(sorted(w[:2])) for w in flags).values())=={4}
print('four vertex charts of six, six edge charts of four')

balanced={(v,frozenset(P)) for v in range(4) for P in combinations([p for p in flags if p[0]==v],3)}
def relabel(z, p):
    v,P=z
    return p[v],frozenset(tuple(p[x] for x in w) for w in P)
seen=set(); orbits=[]
for z in balanced:
    if z not in seen:
        orbit={relabel(z,p) for p in flags}
        orbits.append(len(orbit));seen |= orbit
assert sorted(orbits)==[8,24,24,24]
print('80 balanced rooted charts, orbit sizes:', sorted(orbits))

def mask(edges):
    return sum(1<<E.index(e) for e in edges)
cut={mask(e for e in E if (S>>e[0]&1) != (S>>e[1]&1)) for S in range(16)}
def faces(x):
    return sum((((x&mask(combinations(f,2))).bit_count()%2)<<j) for j,f in enumerate(F))
fiber=defaultdict(Counter); orbit=Counter()
for x in range(64):
    def arc(a,b):
        return (a<b) != bool(x>>E.index(tuple(sorted((a,b))))&1)
    h=sum(all(arc(a,b) for a,b in zip(p,p[1:])) for p in flags)
    out=tuple(sorted(sum(arc(a,b) for b in range(4) if b!=a) for a in range(4)))
    t=sum(arc(a,b)==arc(b,c)==arc(c,a) for a,b,c in F)
    assert h==1+2*t
    fiber[faces(x)][h]+=1
    orbit[out,h]+=1
    assert all(faces(x^y)==faces(x) for y in cut)
assert len(cut)==8 and len(fiber)==8
assert all(sum(n.values())==8 for n in fiber.values())
print('sorted outdegrees, H, labeled multiplicity:',sorted((out,h,n) for (out,h),n in orbit.items()))
print('H distributions by four face-parity bits:',[(format(f,'04b'),sorted(h.items())) for f,h in sorted(fiber.items())])

# F8 using theta^3+theta^2+1; theta orbit forms a normal basis.
def mul(a,b):
    c=0
    while b:
        if b&1:c ^= a
        b >>= 1; a <<= 1
        if a&8:a ^= 13
    return c
basis=[2, mul(2,2), mul(mul(2,2),mul(2,2))]
assert basis[0]^basis[1]^basis[2]==1
def normal(x):
    r=0
    for i,t in enumerate(basis):
        if x>>i&1:r ^= t
    return r
assert len({normal(x) for x in range(8)})==8
for x in range(8):
    rotate=((x<<1)&7)|(x>>2)
    assert mul(normal(x),normal(x))==normal(rotate)
print('normal basis:',basis, 'fixed point in flag coordinates:',7)
