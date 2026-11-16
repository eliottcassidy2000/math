from functools import cache
from itertools import combinations

SEEDS=(
 ((42,12,1,4,15,30),21,75,15,51,61),
 ((34,16,3,7,13,30),6,77,33,47,65),
 ((34,12,1,5,15,30),20,79,19,51,61),
 ((48,13,1,5,14,30),10,81,23,49,63)
)
def calc(u,root=None):
    n=len(u);F=(1<<n)-1
    @cache
    def P(s,i):
        t=s^(1<<i)
        return sum(P(t,j) for j in range(n) if (t&u[i])>>j&1) if t else 1
    @cache
    def H(s):
        return sum(P(s,j) for j in range(n) if s>>j&1) if s else 1
    if root is None: out=H(F)
    else:
        q=[0,0]
        for s in range(1,F+1):
            if s>>root&1:
                q[s.bit_count()%2]+=P(s,root)*H(F^s)
        out=(H(F),H(F^(1<<root)),*q)
    P.cache_clear();H.cache_clear()
    return out
def join(A,m,k,ins):
    a=len(A); tail=((1<<k)-1)<<(a+1)
    out=[v|tail|((1<<a) if not (m>>i&1) else 0) for i,v in enumerate(A)]
    out.append(m|(tail^(ins<<(a+1))))
    out.extend((tail&~((1<<(a+2+j))-1))|((1<<a) if ins>>j&1 else 0) for j in range(k))
    return tuple(out)
def tied(A,m):
    a=len(A)
    return tuple(v|((1<<a) if not(m>>i&1) else 0) for i,v in enumerate(A))+(m,)
for A,m,f,h,R,S in SEEDS:
    assert calc(tied(A,m),6)==(f,h,R,S)
    for k in range(2,7):
        assert calc(join(A,m,k,1<<(k-1)))==f+112*2**(k-1)
for k in range(1,5):
    for ins in range(1<<k):
        offsets=[calc(join(A,m,k,ins))-f for A,m,f,h,R,S in SEEDS]
        assert (len(set(offsets))==1)==(ins==0 or ins.bit_count()==1 and ins!=1)
print('All literal four-state windows and passive-shield classification controls passed.')
