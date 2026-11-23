from collections import Counter
from itertools import combinations
from math import gcd,lcm

P=(3,5,11,17,41)
need=Counter()
for p in P:
    k=next(k for k in range(1,p) if pow(2,k,p) in (1,p-1))
    assert pow(2,k,p)==p-1
    need[k]+=(p-1)//(2*k)
assert need==Counter({1:1,2:1,4:2,5:1,10:2})
assert sum(k*n for k,n in need.items())==36

def partitions(n,k=1):
    if n==0:
        yield ()
    for first in range(k,n+1):
        for rest in partitions(n-first,first):
            yield (first,)+rest

def profiles(parts):
    stay,twist=Counter(),Counter()
    for a in parts:
        stay[a]+=(a-1)//2
        if a%2==0:
            twist[a//2]+=1
    for a,b in combinations(parts,2):
        stay[lcm(a,b)]+=gcd(a,b)
    return +stay,+twist

answers={}
for n in range(1,16):
    for parts in partitions(n):
        stay,twist=profiles(parts)
        if need<=stay+twist:
            answers.setdefault("unordered",(n,parts))
        if need<=stay:
            answers.setdefault("untwisted",(n,parts))
assert answers=={"unordered":(13,(2,2,4,5)), "untwisted":(15,(1,1,2,2,4,5))}
parts=tuple(sorted(2*k for k,n in need.items() for _ in range(n)))
assert sum(parts)==72 and profiles(parts)[1]==need
print(answers,"signed:",parts)

from itertools import permutations
from math import log2,prod,comb
joint=prod(p**((p-1)//2) for p in P)
assert 2**163<joint<2**164
assert sum(2**comb(n-1,2) for n in range(1,20))<2**154<joint
assert 2**comb(19,2)>joint
assert prod(P)==115005

def max_height(n):
    E=list(combinations(range(n),2))
    edge={(a,b):(1<<k,0) for k,(a,b) in enumerate(E)}
    edge|={(b,a):(0,1<<k) for k,(a,b) in enumerate(E)}
    patterns=[]
    for path in permutations(range(n)):
        yes=no=0
        for a,b in zip(path,path[1:]):
            u,v=edge[a,b]
            yes|=u
            no|=v
        patterns.append((yes,no))
    return max(sum(m&yes==yes and not m&no for yes,no in patterns)
               for m in range(1<<len(E)))

assert [max_height(n) for n in (4,5)]==[5,15]
assert comb(9,4)*5*15==9450<2*11*17*41
print(f"joint bits: {log2(joint):.6f}; exact lifting bound: 9450 < 15334")
