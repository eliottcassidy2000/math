"""Independent path DP for the mixed cycle inequality and literal arc port."""
from itertools import product,permutations,combinations
from math import comb,prod
from random import Random

def need(b):
    if not b:raise RuntimeError("cycle/port mismatch")

def subset_h(n,E):
    tab=[[0]*n for _ in range(1<<n)]
    h=[0]*(1<<n);h[0]=1
    for i in range(n):tab[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,v in enumerate(tab[mask]):
            if v:
                for j in range(n):
                    if not mask>>j&1 and (i,j) in E:
                        tab[mask|(1<<j)][j]+=v
        h[mask]=sum(tab[mask])
    return h

fixed={(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)}
cross=list(product((1,2),(3,4)))
def graph(mask):
    return fixed|{p if mask>>i&1 else p[::-1] for i,p in enumerate(cross)}

def inflate(n,E,v,m,D):
    M=[([i] if i<v else list(range(v,v+m)) if i==v else [i+m-1]) for i in range(n)]
    return {(x,y) for i,j in E for x in M[i] for y in M[j]}|{(v+i,v+j) for i,j in D}

def signature(n,E,v):
    out=[0]*(n+1)
    for p in permutations(x for x in range(n) if x!=v):
        b=0
        a=int((v,p[0]) in E)+int((p[-1],v) in E)
        for x,y in zip(p,p[1:]):
            fits=(x,v) in E and (v,y) in E
            if (x,y) not in E:
                if not fits:break
                b+=1
            else:a+=fits
        else:
            for k in range(a+1):out[b+k]+=comb(a,k)
    return tuple(out[:3]),sum(out[3:])

need({subset_h(5,graph(i))[-1] for i in range(16)}=={9,11,13,15})
need(signature(5,graph(3),2)==((5,11,6),0))
need(signature(5,graph(1),2)==((5,13,6),0))
rng=Random(20261115)
def random_t(n):
    return {p if rng.randrange(2) else p[::-1] for p in combinations(range(n),2)}

for m in range(1,6):
    for _ in range(10):
        D=random_t(m)
        h=subset_h(m,D);p=h[-1];q=sum(a*b for a,b in zip(h,reversed(h)))
        need(subset_h(m+4,inflate(5,graph(3),2,m,D))[-1]==6*q-p)
        need(subset_h(m+4,inflate(5,graph(1),2,m,D))[-1]==6*q+p)
for k in (2,3):
    n=1+3*k
    for _ in range(25):
        E=random_t(n);a=[];b=[];e=[]
        for i in range(k):
            u=1+3*i;v=u+1
            for p in ((0,u),(u,v),(v,0)):E.discard(p[::-1]);E.add(p)
        h=subset_h(n,E)
        for i in range(k):
            u=1+3*i;mask=7<<u
            a.append(h[mask]);b.append(h[4<<u]);e.append(h[mask|1]-a[-1])
        rhs=2*prod(x+y for x,y in zip(a,b))-prod(a)
        rhs+=sum((e[i]-2*b[i])*prod(a[j] for j in range(k) if i!=j) for i in range(k))
        rhs+=2*sum(b[i]*b[j]*prod(a[l] for l in range(k) if l not in (i,j))
                   for i,j in combinations(range(k),2))
        need(h[-1]>=rhs)

w="1001100"; n=len(w)
E=set(combinations(range(n),2))|{(n,i) if w[i]=='1' else (i,n) for i in range(n)}
c={}
for i,j in combinations(range(n),2):
    if w[i]=='1' and w[j]=='0':
        for s in range(0,j-i,2):c[s+3]=c.get(s+3,0)+comb(j-i-1,s)
need(c=={3:8,5:17,7:6} and subset_h(n+1,E)[-1]==63)
print("16 local bouquets, 100 literal arc-port counts, 50 larger bouquets, and the 63 cycle ledger")
