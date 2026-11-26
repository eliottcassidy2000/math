from itertools import combinations, product
from functools import lru_cache
from collections import Counter
from math import comb
import random

def S(x): return (3*x+1)//2 if x&1 else x//2
def word(x,k):
    w=[]
    for _ in range(k):
        w.append(x&1);x=S(x)
    return tuple(w),x

for w in product(range(2), repeat=4):
    if sum(w)!=2:continue
    b=0
    for j,d in enumerate(w):b=3**d*b+d*2**j
    r,= [r for r in range(16) if word(r,4)[0]==w]
    for x in range(r,128,16): assert word(x,4)[1]*16==9*x+b
    print('four-bit', ''.join(map(str,w)),b,r,sorted({word(r+16*k,4)[1]%12 for k in range(4)}))
for x in range(3,1000,4):assert 4*(word(x,2)[1]+1)==9*(x+1)
for x in range(1,144,2):
    assert x*x%12==(9 if x%3==0 else 1)
    assert x**4%24==x*x%24

def adj(n,x,cycle=False):
    a=[[0]*n for _ in range(n)]
    for i,j in combinations(range(n),2):
        if cycle and j==i+1:b=0
        elif cycle and (i,j)==(0,n-1):b=1
        else:b=x&1;x>>=1
        a[i][j]=1-b;a[j][i]=b
    return a

def mul(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(a))) for j in range(len(a))] for i in range(len(a))]

def supports(a,k):
    n=len(a)
    def walk(r,s):
        if len(s)==k:return a[s[-1]][r]
        return sum(walk(r,s+[v]) for v in range(r+1,n) if v not in s and a[s[-1]][v])
    return sum(walk(r,[r]) for r in range(n))

def subpaths(a):
    n=len(a);D=[[0]*n for _ in range(1<<n)];h=[1]*(1<<n)
    for s in range(1,1<<n):
        bits=[i for i in range(n) if s>>i&1]
        for j in bits:
            D[s][j]=1 if len(bits)==1 else sum(D[s^(1<<j)][i]*a[i][j] for i in bits)
        h[s]=sum(D[s])
    return h[-1],sum(q*h[-1^i] for i,q in enumerate(h))

def test(a):
    n=len(a);b=mul(a,a);a4=mul(b,b);a5=mul(a4,a)
    assert sum(a5[i][i] for i in range(n))==5*supports(a,5)
    C=[s for s in combinations(range(n),3) if all(sum(a[i][j] for j in s)==1 for i in s)]
    e=[b[j][i] for i in range(n) for j in range(n) if a[i][j]]
    t=[sum(a[i][j]*b[j][i] for j in range(n)) for i in range(n)]
    c=sum(e)//3
    assert c==len(C) and c*3==sum(t)
    assert comb(c,2)-sum(comb(x,2) for x in t)+sum(comb(x,2) for x in e)==sum(not(set(u)&set(v)) for u,v in combinations(C,2))
    h,p=subpaths(a)
    s=2+(n%2==0);k=n if n%2 else n-1
    assert p%(2**(s+1))==2**s*(supports(a,k)%2) and p%3==(-1)**n*h%3
    return h,p,supports(a,k)
rng=random.Random(1209)
total=0
for n in range(3,9):
    xs=range(1<<comb(n,2)) if n<=5 else [rng.randrange(1<<comb(n,2)) for _ in range(100)]
    for x in xs:test(adj(n,x));total+=1
print('independent matrix, cycle, subset path controls',total)
for m in [12,22,25]:print('height 33',m,test(adj(6,m,True)))

@lru_cache(None)
def words(h,lo=2):
    if h==1:return ((),)
    return tuple((d,)+w for d in range(lo,h+1) if h%d==0 and d not in (7,21) for w in words(h//d,d))
for a,b,e in product(range(4),range(3),range(2,6)):
    h=3**a*5**b*7**e;M=a+b;k=e//2;g=M+k
    actual={w for w in words(h) if len(w)>=g}
    base=[3]*a+[5]*b
    def w(s):return tuple(sorted(s))
    if e%2==0: expected={w(base+[49]*k)}
    else:
        expected={w(base+[49]*(k-1)+[343])}
        if b:expected.add(w([3]*a+[5]*(b-1)+[49]*k+[35]))
    assert actual==expected
print('top-grade signatures verified in 48 arithmetic models; uses no assumed prime classes beyond 3 and 5 for those grades')
