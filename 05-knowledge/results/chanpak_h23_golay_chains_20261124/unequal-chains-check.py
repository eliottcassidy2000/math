from itertools import permutations
from collections import Counter
from math import comb
R=(26,52,57,50,32,1)
counts=[]
for i in [1,2]:
    counts.append([sum(all((R[u]>>v)&1 for u,v in zip(w,w[1:]))
                       for w in set(permutations([0,1,3,5]+[2]*i+[4]*j)))
                   for j in [1,2]])
assert counts==[[23,9],[11,4]]
def H(a,b):return (1<<(a+b+2))+3*(1<<a)+(1<<b)-1
def inv(h):
    s=h+1;k=s.bit_length()-1;t=(s&-s).bit_length()-1
    return next(((a,b) for a,b in [(t,k-2-t),(k-2-t,t),(t-2,t-2)]
                 if min(a,b)>=1 and H(a,b)==h),None)
for k in range(4,100):
    shell={H(a,k-2-a) for a in range(1,k-2)}
    assert len(shell)==k-3
    assert max(shell)+1==(1<<k)+3*(1<<(k-3))+2
    for h in shell:
        assert (h+1).bit_length()-1==k and inv(h) and H(*inv(h))==h
    assert all(inv(h) is None for h in [(1<<k)-1,(1<<k)+1])
print(counts,'dyadic shells passed')
g=sum(1<<i for i in (11,9,7,6,5,1,0))
C=[0]
for i in range(12):C += [w^(g<<i) for w in C]
assert len(set(C))==4096
weight=Counter(w.bit_count() for w in C)
assert weight=={0:1,7:253,8:506,11:1288,12:1288,15:506,16:253,23:1}
assert len({H(1+w,24-w) for w in weight})==8
assert len(C)*sum(comb(23,i) for i in range(4))==1<<23
print(sorted(weight.items()))
expected = {0,1,2,3,5,6}
assert [[H(a,b)%7 for b in (3,1,2)] for a in (3,1,2)]==[[0,5,1],[0,2,6],[0,3,2]]
for j in range(1,5):
    q=7**j
    got=set()
    for a in range(1,4):
        x=1
        for b in range(1,3*7**(j-1)+1):
            x=2*x%q
            got.add(((4*(1<<a)+1)*x+3*(1<<a)-1)%q)
    assert got=={r for r in range(q) if r%7 in expected}
print('7-adic image passed')
