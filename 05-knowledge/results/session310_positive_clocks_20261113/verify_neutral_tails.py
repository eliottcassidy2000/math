from itertools import product
from math import isqrt
A=(2,28,25,81,97,79,7)
B=(2,80,3,7,173,143,189,15)
def require(ok):
    if not ok: raise ValueError("certificate failed")
def profile(a):
    n=len(a); m=(1<<n)-1
    h=[[0]*n for _ in range(m+1)]
    for i in range(n):
        require(not (a[i]>>i&1))
        for j in range(i):require((a[i]>>j&1)+(a[j]>>i&1)==1)
        h[1<<i][i]=1
        s=1<<i
        for _ in a:
            for j,z in enumerate(a):
                if s>>j&1:s|=z
        require(s==m)
    for s in range(1,m+1):
        for i,v in enumerate(h[s]):
            if not v:continue
            t=a[i]&~s
            while t:
                b=t&-t;t^=b;h[s|b][b.bit_length()-1]+=v
    hp=list(map(sum,h))
    f=[[0]*(n+1) for _ in range(m+1)];f[0][0]=1
    for s in range(1,m+1):
        t=s
        while t:
            for j in range(1,n+1):f[s][j]+=hp[t]*f[s^t][j-1]
            t=(t-1)&s
    return f[m]
fa,fb=profile(A),profile(B)
require(fa[1:]==[95,1086,5166,13200,19200,15120,5040])
require(fb[1:]==[95,2362,18756,73464,162000,205920,141120,40320])
require(all(fb[j]==j*((fa[j] if j<len(fa) else 0)+fa[j-1]) for j in range(1,9)))
def tt(n):
    f=[1]+[0]*n
    for _ in range(n):f=[0]+[j*(f[j]+f[j-1]) for j in range(1,n+1)]
    return f
def cyc(a,b):
    z=0;f,g=tt(a),tt(b)
    for i,j,k in product(range(1,9),range(1,a+1),range(1,b+1)):
        for start in range(3):
            word=[(start+t)%3 for t in range(i+j+k)]
            if [word.count(x) for x in range(3)]==[i,j,k]:
                z+=fb[i]*f[j]*g[k]
    return z
def prime(n):
    return n>1 and all(n%d for d in range(2,isqrt(n)+1))
require(cyc(1,1)==2647 and prime(2647))
p=228334402816633
require(cyc(5,10)==7*p)
q=64720635719
require(p-1==2**3*3**2*7**2*q and prime(q))
require(pow(5,p-1,p)==1 and all(pow(5,(p-1)//d,p)!=1 for d in (2,3,7,q)))
print('PASS: full cover polynomial identity and prime / seven-prime releases')

from collections import Counter
from itertools import combinations, permutations
def packing_types(a):
    n=len(a);cy=[]
    for j in range(3,n+1,2):
        for v in combinations(range(n),j):
            for tail in permutations(v[1:]):
                w=(v[0],)+tail+(v[0],)
                if all(a[x]>>y&1 for x,y in zip(w,w[1:])):
                    cy.append((sum(1<<x for x in v),j))
    c=Counter()
    def go(i,used,L):
        c[tuple(sorted(L))]+=1
        for z in range(i,len(cy)):
            m,j=cy[z]
            if not m&used:go(z+1,m|used,L+(j,))
    go(0,0,())
    return c
require(packing_types(A)=={():1,(3,):10,(5,):18,(7,):7,(3,3):6})
require(packing_types(B)=={():1,(3,):10,(5,):30,(7,):7})
print('PASS: two pentagons replace one disjoint triangle pair in the full profile')
