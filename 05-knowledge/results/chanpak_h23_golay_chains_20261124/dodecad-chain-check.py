from collections import Counter
from itertools import permutations
from math import comb

T = (26,52,57,50,32,1)

def arc(i,j):
    return T[i] >> j & 1

assert all(not (m & (m-1)) or any(not (m >> i & 1) and (T[i] & m) not in (0,m)
           for i in range(6)) for m in range(1,63))
assert sorted(x.bit_count() for x in T) != sorted(5-x.bit_count() for x in T)
assert all(any(arc(i,j)!=arc(p[i],p[j]) for i in range(6) for j in range(6))
           for p in permutations(range(6)) if p!=tuple(range(6)))
for i,j,z in ((1,1,23),(1,2,9),(2,1,11),(2,2,4)):
    words=set(permutations([0,1,3,5]+[2]*i+[4]*j))
    assert sum(all(arc(a,b) for a,b in zip(w,w[1:])) for w in words)==z

def read(g):
    n=len(g)
    def h(S):
        return sum(all(g[a]>>b&1 for a,b in zip(w,w[1:]))
                   for w in permutations(i for i in range(n) if S>>i&1))
    return h((1<<n)-1),sum(h(s)*h(((1<<n)-1)^s) for s in range(1<<n))
assert read((6,4,0))==(1,8) and read((2,4,1))==(3,12)

gen=sum(1<<i for i in (0,2,4,5,6,10,11))
C={0}
for i in range(12):
    v=gen<<i
    v|=(v.bit_count()%2)<<23
    C|={x^v for x in C}
d=0x6fdd
assert d in C and d.bit_count()==12
assert Counter(x.bit_count() for x in C)=={0:1,8:759,12:2576,16:759,24:1}
assert 759*comb(8,5)==comb(24,5)
profiles=Counter(((c&d).bit_count(),(c&~d).bit_count()) for c in C)
assert profiles=={(0,0):1,(2,6):132,(4,4):495,(6,2):132,
                  (0,12):1,(4,8):495,(6,6):1584,(8,4):495,(12,0):1,
                  (6,10):132,(8,8):495,(10,6):132,(12,12):1}
for O in C:
    if O.bit_count()==8:
        assert Counter(c&O for c in C)=={S:32 for S in C|{c&O for c in C}
                                             if S&O==S and S.bit_count()%2==0}

def h(a,b): return 4*2**(a+b)+3*2**a+2**b-1

def inverse(H):
    u=H+1
    s=u.bit_length()-3
    t=(u&-u).bit_length()-1
    return next(((a,b) for a,b in ((t,s-t),(s-t,t),(t-2,t-2))
                 if a>=1 and b>=1 and h(a,b)==H),None)

for k in (36,37,100):
    heights=set()
    for a,b in profiles:
        w=a+b
        H=3**w*h(k-a,k-b)
        assert H%16==15 and H//3**w%3!=0
        assert inverse(H//3**w)==(k-a,k-b)
        if w:
            assert H%3==0
        else:
            assert H==(2**(k+1)+1)**2-2
        heights.add(H)
    assert len(heights)==13
print("checked the 13 split profiles, labelled octad differences, and exact inverse")
