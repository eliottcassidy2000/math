from itertools import product
from math import prod
P=(3,5,7,13,17,241); M=(2,4,3,12,8,24)
assert prod(P)*3==2**24-1
assert all(next(i for i in range(1,25) if pow(2,i,p)==1)==m for p,m in zip(P,M))
mask=lambda r,m:sum(1<<i for i in range(r,24,m))
D=[[mask(r,m) for r in range(m)] for m in M]
cover=[]
for A in product(*map(range,M)):
    S=[v[r] for v,r in zip(D,A)]
    t=0
    for u in S:t|=u
    if t==(1<<24)-1:
        assert all(s & ~sum(1<<r for r in range(24) if any(j!=i and T>>r&1 for j,T in enumerate(S))) for i,s in enumerate(S))
        cover.append(A)
reps={min(tuple((r+d)%m for r,m in zip(A,M)) for d in range(24)) for A in cover}
assert reps=={(0,1,0,7,3,23),(0,1,0,11,3,7)} and len(cover)==48
for A,k in zip(sorted(reps),(3968308,3177553)):
    assert all(k*pow(2,r,p)%p==1 for p,r in zip(P,A))
assert all(any((509203*pow(2,e,p)-1)%p==0 for p in P) for e in range(24))
assert all(any((78557*pow(2,e,p)+1)%p==0 for p in (3,5,7,13,19,37,73)) for e in range(36))
assert [e for e in range(24) if all((13*pow(2,e,p)+1)%p for p in P)]==[2,4,8,10,16,20]
for p,m,a in [(7,3,0),(9,6,1),(49,21,6),(343,147,69)]:
    assert pow(2,m,p)==1 and [e for e in range(m) if (13*pow(2,e,p)+1)%p==0]==[a]
print('Two normalized CRT wheels, 48 covers, 13-ray phases verified.')
