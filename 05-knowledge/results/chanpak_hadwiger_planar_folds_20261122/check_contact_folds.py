from fractions import Fraction as Q
from itertools import permutations

def H(a):
    n=len(a);M=1<<n
    d=[[0 for _ in range(n)] for _ in range(M)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,M):
        for i,h in enumerate(d[s]):
            if not h:continue
            m=a[i]&~s
            while m:
                t=m&-m;m^=t
                d[s|t][t.bit_length()-1]+=h
    return sum(d[-1])

def flip(a,i,j):
    return tuple(x^(1<<j if k==i else 1<<i if k==j else 0)
                 for k,x in enumerate(a))
def extend(a,s):
    n=len(a);return tuple(x|(((s>>i&1)^1)<<n) for i,x in enumerate(a))+(s,)
def det(a,b):return a[0]*b[1]-a[1]*b[0]
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def mul(z,a):return tuple(z*x for x in a)
def sub(a,b):return add(a,mul(-1,b))

BANK=[
 ((16,1,3,7,14),(0,3),0,9,((0,0),(4,6),(12,9),(14,13)),(-36,-32,-30,-26)),
 ((8,17,3,6,13),(0,1),16,13,((12,9),(12,14),(15,9),(12,11)),(-15,0,-6,9)),
 ((2,24,3,5,13),(0,3),16,11,((13,12),(13,16),(17,12),(13,12)),(-16,0,0,16)),
 ((2,36,17,7,11,29),(4,5),14,29,((43,43),(44,52),(50,42),(43,43)),(-64,0,0,64)),
]
for a,(i,j),s,h,corners,jd in BANK:
    b=flip(a,i,j)
    assert (H(a),H(b))==(h,h+2)
    got=[]
    for m in (s,s|1<<i,s|1<<j,s|1<<i|1<<j):
        got.append(((H(extend(a,m))-h)//2,(H(extend(b,m))-h-2)//2))
    assert tuple(got)==corners
    P,A,B,C=corners[0],sub(corners[1],corners[0]),sub(corners[2],corners[0]),\
            add(sub(sub(corners[3],corners[1]),corners[2]),corners[0])
    jac=lambda x,y:det(add(A,mul(y,C)),add(B,mul(x,C)))
    F=lambda x,y:add(P,add(mul(x,A),add(mul(y,B),mul(x*y,C))))
    assert tuple(jac(x,y) for x,y in ((0,0),(1,0),(0,1),(1,1)))==jd
    if h==13:
        assert F(Q(16,25),Q(2,5))==F(Q(19,25),Q(3,5))==(Q(1554,125),Q(1429,125))
    if h==11:
        assert F(0,0)==F(1,1) and jac(0,0)*jac(1,1)<0
    if h==29:
        assert F(0,0)==F(1,1)==(43,43) and jac(0,0)*jac(1,1)<0
        X,Y=(extend(a,m) for m in (14,62)); U,V=(extend(b,m) for m in (14,62))
        def same(a,b,p):
            return all((a[i]>>j&1)==(b[p[i]]>>p[j]&1)
                       for i in range(7) for j in range(7))
        perms=list(permutations(range(7)))
        assert any(same(X,Y,p) for p in perms) and any(same(U,V,p) for p in perms)
        assert not any(same(X,Y,p) and same(U,V,p) for p in perms)
print('direct Hamilton DP verifies all paired faces, both types of collision, and the signed Jacobians')
