from itertools import product
from math import isqrt
from importlib import import_module
C=import_module('backedge-chart')

def paths(A):
    n=len(A);N=1<<n
    D=[[0]*n for _ in range(N)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,N):
        for v,h in enumerate(D[s]):
            if not h:continue
            e=A[v]&~s
            while e:
                b=e&-e;e-=b;D[s|b][b.bit_length()-1]+=h
    return [1]+list(map(sum,D[1:]))
def profile(A):
    p=paths(A);F=(1<<len(A))-1
    return p[F],sum(h*p[F^s] for s,h in enumerate(p))
def connected(A):
    F=(1<<len(A))-1
    for i in range(len(A)):
        S=1<<i
        while True:
            U=S
            for j,a in enumerate(A):
                if S>>j&1:U|=a
            if U==S:break
            S=U
        if S!=F:return False
    return True
def reverse(A,i,j):
    q=list(A);q[i]^=1<<j;q[j]^=1<<i;return tuple(q)
def trans(w):
    return tuple(((1<<w)-1)^((1<<(i+1))-1) for i in range(w))
def inflate(A,z,B):
    P=[(i,j) for i in range(len(A)) for j in range(len(B) if i==z else 1)]
    return tuple(sum(1<<k for k,(t,a) in enumerate(P)
                 if ((A[i]>>t)&1) or (i==t==z and (B[j]>>a)&1)) for i,j in P)

a=(4,9,34,21,71,27,47)
b=reverse(a,3,5)
assert connected(a) and connected(b)
for R,d in ((a,3),(b,5)):
    assert [v for v,r in enumerate(R) if r>>6&1]==[4]
    for w in range(1,5):
        A=inflate(R,6,trans(w))
        assert connected(A) and profile(A)[0]==(1<<(w+5))+d
    U=(2,4,1)
    assert profile(inflate(R,6,U))[0]==32*12+d*3

D=(
 (5,5,(2,8,3,5)),
 (19,8,(2,16,11,19,37,15)),
 (217,12,(2,24,19,133,201,31,47,103)),
 (223,12,(2,136,11,145,71,31,47,117)),
 (599,13,(28,73,146,148,98,15,45,115)),
)
for q,h,R in D:
    assert connected(R) and profile(R)==(q,1<<h)
    A=inflate((2,4,1),0,R)
    assert connected(A) and profile(A)[0]==q+(1<<h)
assert 217==7*31
for q in (5,19,7,31,223,599):
    assert all(q%d for d in range(2,isqrt(q)+1))

nested=((0,3),(1,2));cross=((0,2),(1,3))
for pattern in (nested,cross):
    P=C.compile(4,pattern)
    for x,y,z in product(range(3),repeat=3):
        H=C.value(P,(x,y,z))
        assert H==C.literal((x,y,z),pattern)
        if pattern==nested and y>=1:
            assert H==1+(1<<y)+(1<<(x+z+1))*(3**y+(1<<(y+1))-1)
        if pattern==cross and not y and min(x,z)>0:
            assert H==2*((1<<x)+1)*((1<<z)+1)-1
base=((0,2),(1,4));quad=base+((3,4),)
for w in range(5):
    gap=(w+1,0,1,1)
    h0=C.literal(gap,base);h1=C.literal(gap,quad)
    assert (h0,h1)==((36<<w)+17,(36<<w)+15)
for pattern,n in ((quad,5),(((0,2),(1,3),(2,4)),5),(((0,1),(1,2),(0,2)),3)):
    T=C.compile(n,pattern)
    for w in product(range(2),repeat=n-1):assert C.value(T,w)==C.literal(w,pattern)
extra=(
 (259, 13, (256, 257, 259, 71, 47, 335, 407, 319, 24)),
 (803, 14, (256, 1, 27, 147, 99, 79, 143, 311, 126)),
 (695, 14, (256, 257, 43, 307, 327, 147, 303, 95, 164)),
 (713, 14, (256, 261, 289, 7, 15, 411, 63, 95, 216)),
 (761, 14, (256, 5, 49, 39, 331, 339, 143, 319, 78)),
 (773, 14, (256, 269, 297, 193, 271, 283, 55, 375, 72)),
 (707, 14, (256, 273, 267, 451, 141, 287, 311, 103, 144)),
 (815, 14, (256, 21, 25, 67, 41, 79, 407, 319, 62)),
 (1957, 15, (258, 36, 401, 295, 75, 405, 47, 91, 210)),
 (1969, 15, (258, 132, 353, 23, 39, 331, 155, 61, 218)),
 (1945, 15, (260, 69, 112, 327, 43, 139, 177, 287, 118)),
)
private=(37,73,139,23,761,773,101,163,103,179,389)
for (h,q,A),p in zip(extra,private):
    assert connected(A) and profile(A)==(h,1<<q)
    assert h%p==0 and all(k%p for k,_,_ in D+extra if k!=h)
    assert all(p%j for j in range(2,isqrt(p)+1))
for h,q,A in D+extra:
    assert h**100 < 1<<(73*q)
for w in range(1,5):
    for A,t in ((a,18),(b,19)):
        expected=8*(11*(3**(w+1)-1)+t*(1<<w))
        assert profile(inflate(A,6,trans(w)))[1]==expected
print('unit rails, sixteen dyadic seeds, two-slot height modes, crossing interpolation, and descent compiler pass')
