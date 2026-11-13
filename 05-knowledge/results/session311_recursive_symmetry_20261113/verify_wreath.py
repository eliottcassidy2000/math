from itertools import permutations, product
from collections import Counter

def require(x):
    if not x: raise ValueError("wreath replay failed")

B=(2,80,3,7,173,143,189,15)
D=(130,156,153,209,225,207,135,0)
def aut(rows):
    n=len(rows)
    return sum(all((rows[i]>>j&1)==(rows[p[i]]>>p[j]&1)
                   for i in range(n) for j in range(i+1,n))
               for p in permutations(range(n)))
require(aut(B)==aut(D)==1)

def canon(c):
    if isinstance(c,int):return c
    t=tuple(map(canon,c))
    return min(t,t[1:]+t[:1],t[2:]+t[:2])
def exponent(c):
    if isinstance(c,int):return 0
    return sum(map(exponent,c))+int(c[0]==c[1]==c[2])
def flip_first(c):
    if isinstance(c,int):return 1-c
    return canon((flip_first(c[0]),c[1],c[2]))
A={0,1}
expected=[2,4,24,4624]
for d,n in enumerate(expected):
    h=Counter(map(exponent,A))
    require(len(A)==n and h[0]==2 and set(h)==set(range((3**d+1)//2)))
    mono=0
    for _ in range(d):mono=(mono,)*3
    require(exponent(mono)-exponent(flip_first(mono))==d)
    if d+1<len(expected):
        AA={canon(c) for c in product(A,repeat=3)}
        hh=Counter()
        for (i,a),(j,b),(k,c) in product(h.items(),repeat=3):hh[i+j+k]+=a*b*c
        for i,a in h.items():hh[3*i]-=a
        require(all(x%3==0 for x in hh.values()))
        hh=Counter({i:x//3 for i,x in hh.items() if x})
        for i,a in h.items():hh[3*i+1]+=a
        require(hh==Counter(map(exponent,AA)))
        A=AA
print('PASS: rigid seeds, four exact class inventories, group interval and sharp single-leaf change')
