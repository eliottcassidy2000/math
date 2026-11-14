from collections import defaultdict

def ending(A):
    n=len(A); F=(1<<n)-1
    out=[[0]*n for _ in range(F+1)]
    for i in range(n):out[1<<i][i]=1
    for S in range(1,F+1):
        for i,x in enumerate(out[S]):
            if x:
                B=A[i]&~S
                while B:
                    b=B&-B;B^=b
                    out[S|b][b.bit_length()-1]+=x
    return out
def count(A):
    return sum(ending(A)[-1]) if A else 1
def rev(x,n):
    return int(f'{x:0{n}b}'[::-1],2)
def g(x):
    h=1
    while x:
        x,b=divmod(x,2)
        if not b:h+=x+(x&1)
    return h

def profile(A,U):
    n=len(A);F=(1<<n)-1
    ee=ending(A)
    ss=ending([r^F^(1<<i) for i,r in enumerate(A)])
    h=[1]+[sum(row) for row in ee[1:]]
    f=sum(h[s]*h[F^s] for s in range(F+1))
    P=sum(ss[s][i]*h[F^s] for s in range(F+1) for i in range(n) if U>>i&1)
    Q=sum(ee[s][i]*h[F^s] for s in range(F+1) for i in range(n) if not U>>i&1)
    S=sum(ee[F][i] for i in range(n) if U>>i&1)
    T=sum(ss[F][i] for i in range(n) if not U>>i&1)
    K=count(tuple(r|((((U>>i)&1)^1)<<n) for i,r in enumerate(A))+(U,))
    E=sum(ee[s][i]*ss[F^s][j] for s in range(F+1) for i in range(n) for j in range(n) if U>>i&1 and U>>j&1)
    assert K-P+S==h[F]-E
    return h[F],f,P,Q,S,T,K

def physical(A,U,x,y):
    left=f'{x:b}'
    right=''.join('1' if b=='0' else '0' for b in f'{y:b}'[::-1])
    l,n,r=len(left),len(A),len(right); F=(1<<(l+n+r))-1
    C=[]
    for i,b in enumerate(left):
        C.append((F^((2<<i)-1))|(((b=='0')*1)<<(l+n+r)))
    for i,R in enumerate(A):
        C.append((R<<l)|(((1<<r)-1)<<(l+n))|((((U>>i)&1)^1)<<(l+n+r)))
    for j,b in enumerate(right,l+n):
        C.append((F^((2<<j)-1))|(((b=='0')*1)<<(l+n+r)))
    return C+[sum((b=='1')<<i for i,b in enumerate(left))|(U<<l)|sum((b=='1')<<(l+n+i) for i,b in enumerate(right))]

from itertools import product

def digits(x):
    return f'{x:b}' if x else ''

def phys(bits):
    m=len(bits)
    rows=[((1<<m)-1)^((2<<i)-1) | (int(not int(w))<<m) for i,w in enumerate(bits)]
    rows.append(sum(int(w)<<i for i,w in enumerate(bits)))
    return rows

for t in range(1000):
    for x,y,want in [(10*t,34*t+1,-2),(10*t+9,34*t+32,2)]:
        counts={}
        for a,b,c in product((0,1),repeat=3):
            z=24+33*a-14*b-4*c
            bits=digits(x)+f'{z:06b}'+''.join(str(1-int(w)) for w in digits(y)[::-1])
            base=64*x*y+57*x+24*y+g(x)+g(y)+21+(x&1)
            v=base+a*(33*(y-x)+(y&1)-(x&1))-(14*b+4*c)*(x+y+1)
            assert g(int(bits,2))==v
            if t<=1:
                assert count(phys(bits))==v
            counts[a,b,c]=v
        assert counts[1,1,1]-counts[0,0,0]==want
        if t<=1:
            print(x,y,counts[0,0,0],counts[1,1,1])
print('Verified interleaved contact switches and both exact arithmetic returns.')
from random import Random
rng=Random(61114)
def tour(n):
    A=[0]*n
    for i in range(n):
        for j in range(i):
            if rng.randrange(2):A[i]|=1<<j
            else:A[j]|=1<<i
    return tuple(A)
def ord_apex(parts):
    F=sum(len(A) for A,U in parts); C=[]; M=0; pos=0
    for A,U in parts:
        for i,r in enumerate(A):
            C.append((r<<pos)|(((1<<F)-1)^((1<<(pos+len(A)))-1))|((1^((U>>i)&1))<<F))
        M|=U<<pos;pos+=len(A)
    return C+[M]
for _ in range(40):
    A,B=tour(rng.randrange(1,4)),tour(rng.randrange(1,4))
    ua,ub=rng.randrange(1<<len(A)),rng.randrange(1<<len(B))
    ha,fa,pa,qa,sa,ta,ka=profile(A,ua)
    hb,fb,pb,qb,sb,tb,kb=profile(B,ub)
    m=rng.randrange(3,6)
    C=tuple(((1<<m)-1)^((2<<i)-1) for i in range(m))
    def literalZ(z):return count(ord_apex([(A,ua),(C,rev(z,m)),(B,ub)]))
    base=literalZ(0)
    c={}
    for z in range(1<<m):
        c[z]=literalZ(z)
        assert c[z]==base+ha*qb*z-pa*hb*rev(z,m)+ha*hb*(g(z)-1)+ha*tb*(z&1)-sa*hb*(z>>(m-1))
    for z in c:
        for i in range(m-1):
            if (z>>i)&3==2:
                assert c[z]>c[z-(1<<i)]
print('Verified exact single-spacer polynomial and strict move order with forty arbitrary marked flanks.')
