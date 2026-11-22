from math import prod,factorial
from collections import Counter
from itertools import combinations
from random import Random

MODULES=[
((256,289,67,23,295,141,59,95,236),1321,24800),
((6,72,42,49,71,83,141,63),331,5368),
((768,257,3,7,815,399,831,863,524,46),451,21136),
((256,257,3,23,39,399,319,95,156),305,9268),
((2,4,1),3,12),((2,4,1),3,12),((2,8,3,5),5,32),
((256,257,515,263,783,31,831,127,164,427),1057,38656),
((8,5,17,6,11),13,104),((4,33,18,7,11,29),31,328)]

def Hsubs(a,root=None):
    n=len(a);N=1<<n
    d=[[0]*n for _ in range(N)]
    for i in range(n):
        if root is None or root==i:d[1<<i][i]=1
    for s in range(1,N):
        for i,v in enumerate(d[s]):
            if not v:continue
            m=a[i]&~s
            while m:
                b=m&-m;m^=b
                d[s|b][b.bit_length()-1]+=v
    h=[sum(x) for x in d]
    if root is None:h[0]=1
    return h

def strong(a):
    n=len(a)
    for i in range(n):
        s=1<<i
        while True:
            t=s
            for j in range(n):
                if s>>j&1:t|=a[j]
            if t==s:break
            s=t
        if s!=(1<<n)-1:return False
    return True

def wrap(a):
    n=len(a);return tuple(u|(1<<n) for u in a)+(1<<(n+1),(1<<n)-1)

def op(a):
    allv=(1<<len(a))-1
    return tuple(x^allv^(1<<i) for i,x in enumerate(a))

def flip(a,i,j):
    rr=list(a);rr[i]^=1<<j;rr[j]^=1<<i
    return tuple(rr)

for a,h,q in MODULES:
    H=Hsubs(a)
    assert (H[-1],sum(x*H[-1^i] for i,x in enumerate(H)))==(h,q)
    assert strong(a)
h=prod(h for a,h,q in MODULES);q=2**24*prod(q for a,h,q in MODULES)
assert h==2**60-1 and q==2**61*h
assert h+q==(1<<60)*((1<<61)-1)-1
assert sum(len(a) for a,h,q in MODULES)+24+2==93
m=Counter(a for a,h,q in MODULES)
assert factorial(34)//(factorial(24)*prod(factorial(t) for t in m.values()))==237918897216000
assert 2**61>6*(11*factorial(10))**2
print('93-vertex class family, exact counts and size inequality pass')

rng=Random(44)
for n in range(1,8):
    for _ in range(10):
        a=[0]*n
        for i,j in combinations(range(n),2):
            if rng.randrange(2):i,j=j,i
            a[i]|=1<<j
        a=tuple(a);W=wrap(a);h=Hsubs(a);full=(1<<n)-1
        Q=sum(x*h[full^s] for s,x in enumerate(h))
        base=Hsubs(W)[-1]
        assert base==h[-1]+Q
        for u in range(n):
            st=Hsubs(a,u);en=Hsubs(op(a),u)
            S,E=st[-1],en[-1]
            s=sum(x*h[full^i] for i,x in enumerate(st))
            d=sum(x*h[full^i] for i,x in enumerate(en))
            assert Hsubs(flip(W,u,n))[-1]-base==s-d-S-E
            assert Hsubs(flip(W,n+1,u))[-1]-base==d-s-S-E
print('70 literal controls of the two-slot wrapper and every exterior contact pass')

B91=[
 (3963,5889,4,(512,769,523,547,271,919,63,95,717,208)),
 (331,671,3,(6,72,42,49,71,83,141,63)),
 (1661,6355,3,(256,785,523,547,269,919,575,863,76,273)),
 (1525,6605,2,(256,37,25,147,163,333,31,359,94)),
 (1599,6951,2,(258,36,145,279,99,333,399,59,150)),
 (217,1,12,(2,24,19,133,201,31,47,103))]
for h,o,v,a in B91:
    H=Hsubs(a)
    assert H[-1]==h and sum(x*H[-1^i] for i,x in enumerate(H))==(o<<v)
    assert strong(a)
assert sum(len(a) for h,o,v,a in B91)+35+2==91
assert prod(h for h,o,v,a in B91)==prod(o for h,o,v,a in B91)==2**60-1
assert sum(v for h,o,v,a in B91)+35==61
assert factorial(41)//factorial(35)==3237399360
print('91-vertex refined cyclic family passes')
