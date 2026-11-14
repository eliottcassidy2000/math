from itertools import product
from math import isqrt

def H(T):
    n=len(T)
    e=[[0]*n for _ in range(1<<n)]
    for m in range(1,1<<n):
        z=m
        while z:
            b=z&-z; z-=b; i=b.bit_length()-1; q=m^b
            if not q: e[m][i]=1
            else:
                w=q&T[i]
                while w:
                    k=w&-w;w-=k;e[m][i]+=e[q][k.bit_length()-1]
    return [1]+list(map(sum,e[1:]))

def check(ok):
    if not ok: raise AssertionError("literal symmetry check")

def rotor(A):
    return [sum(1<<((i+a)%7) for a in A) for i in range(7)]

check({H(rotor([j if b else 7-j for j,b in enumerate(z,1)]))[-1]
       for z in product((0,1),repeat=3)}=={175,189})
R=rotor((1,2,3)); h=H(R); full=127
check(h[full]==175 and sum(h[s]*h[full^s] for s in range(128))==1988)
for k,p in ((1,593),(3,2297)):
    check(p>1 and all(p%d for d in range(2,isqrt(p)+1)))
    # Core first, ordinal tape next; portal 0->1->B->0.
    m=7+k
    T=[2,((1<<(m+2))-4)]
    for i in range(m):
        if i<7: x=(R[i]<<2)|(((1<<k)-1)<<9)
        else: x=((1<<m)-(2<<i))<<2
        T.append(1|x)
    check(H(T)[-1]==7*p)
print('Eight 7-rotors, two-cover count, literal 7p witnesses, and trial division passed.')
