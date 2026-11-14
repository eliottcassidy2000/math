from itertools import combinations, product
from random import Random
from frame_decoder import *
from residual_machine import certify

def hk(out,active=None):
    n=len(out)
    target=(1<<n)-1 if active is None else active
    d=[[0]*n for i in range(1<<n)]
    for i in range(n):
        if target>>i&1:
            d[1<<i][i]=1
    for S in range(1<<n):
        for i,count in enumerate(d[S]):
            if count:
                R=out[i]&target&~S
                while R:
                    j=R&-R
                    R^=j
                    d[S|j][j.bit_length()-1]+=count
    return sum(d[target])

def rows_from_bits(n,c):
    o=[0]*n
    for k,(a,b) in enumerate(combinations(range(n),2)):
        if c>>k&1:
            o[a]|=1<<b
        else:
            o[b]|=1<<a
    return o

for n in range(1,6):
    assert not any(hk(rows_from_bits(n,c))==23
        for c in range(1<<(n*(n-1)//2)))
o=adjacency('100110')
assert hk(o)==63
assert [sum(v.bit_count()==i for v in o) for i in range(7)]==[0,1,2,2,0,1,2]
assert [hk(o,255^(1<<i)) for i in range(8)]==[11,37,37,33,33,25,25,1]
for w,h in [('110',11),('1121',25),('1120',33),('1200',37),('1110',23)]:
    assert hk(adjacency(w))==h
assert coefficient_tensor([(1,1),(0,None),(1,2),(0,1)])=={(1,):37,(2,):13,(3,):0}

cases = [
    (385,36,[(1,1),(0,None),(3,1),(0,1),(2,None),(2,1),(1,1)],(1,1),'1030221'),
    (573,108,[(2,1),(0,None),(3,1),(0,1),(2,None),(1,1),(2,1)],(1,2),'20302212')]
for N,M,S,k,w in cases:
    A=(0,2,5,6)
    assert decode(N,S,A,M)==[k]
    assert choose_class(N,[(S,A,M)])==w
    assert hk(adjacency(w))==N
    assert decode(N,S,A,M+2)==[]

# Restricted-radius random frames; every known token passes the target
# inverse after the neutral widths are set to one. If the decoder returns
# more choices outside the sampled radius, check their exact counts too.
r=Random(63)
for _ in range(240):
    S=[(r.randrange(4),r.choice([1,None])) for j in range(r.randint(0,5))]
    d=sum(t is None for _,t in S)
    c=coefficient_tensor(S)
    live=[c[(1,)*j+(2,)+(1,)*(d-j-1)]>0 for j in range(d)]
    k=tuple(r.randint(1,3) if a else 1 for a in live)
    N=H(word(S,k))
    assert evaluate(c,k)==N
    z=decode(N,S)
    assert k in z and all(H(word(S,y))==N for y in z)
certify()
print('23 first order 6, 63 deletions and batch rail checked, 240 block inversions, and residual rank 17.')
