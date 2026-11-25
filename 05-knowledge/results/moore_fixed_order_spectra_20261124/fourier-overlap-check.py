"""Exact capacity, channel gluing on Hoffman–Singleton, and Chan–Pak controls."""
from fractions import Fraction
from itertools import combinations
from math import comb, factorial

def capacities():
    out = []
    for c in (25, 100, 175):
        t = c//5
        r = (2065+t)//15
        M = 8395-c
        Q = Fraction((295-c)*49+t*268,225)
        loss = Fraction(175+72*t,15)
        lower = Fraction(11*r*r,100)-4*Q-Fraction(t*268,225) \
            -Fraction(M+4*(1350-2*c),225)+loss
        L = -(-lower.numerator//lower.denominator)
        R = 57*r-L
        slots = 11*r*(589-r)//2
        a = 1
        total = a
        for j in range(1,R+1):
            a = a*(slots-j+1)//j
            total += a
        # Only the 295-bit basis mask and 2-bit branch precede the rank.
        raw_bits = (total-1).bit_length()+297
        gauge_bits = (total//11**294-1).bit_length()+297
        out.append((c,r,L,R,slots,raw_bits,gauge_bits))
    assert max(x[-1] for x in out)<=comb(295,2)
    assert max(x[-2] for x in out)<=comb(298,2)
    assert factorial(295)<1<<2003
    return out

def rank_subset(xs, n):
    return sum(comb(n,j) for j in range(len(xs))) + sum(comb(x,j) for j,x in enumerate(xs,1))

def unrank_subset(code, n):
    k = 0
    while code>=comb(n,k):
        code-=comb(n,k)
        k+=1
    upper = n-1
    ids = []
    for j in range(k,0,-1):
        lo,hi=j-1,upper
        while lo<hi:
            t=(lo+hi+1)//2
            if comb(t,j)<=code:
                lo=t
            else:
                hi=t-1
        ids.append(lo)
        code-=comb(lo,j)
        upper=lo-1
    assert not code
    return tuple(reversed(ids))

def inverse(a, p):
    n=len(a)
    a=[list(row)+[int(i==j) for j in range(n)] for i,row in enumerate(a)]
    for j in range(n):
        k=next((i for i in range(j,n) if a[i][j]%p),None)
        if k is None:
            return None
        a[j],a[k]=a[k],a[j]
        f=pow(a[j][j],-1,p)
        a[j]=[v*f%p for v in a[j]]
        for i in range(n):
            if i!=j:
                f=a[i][j]
                a[i]=[(g-f*h)%p for g,h in zip(a[i],a[j])]
    return [v[n:] for v in a]

def mm(a,b,p):
    return [[sum(x*y for x,y in zip(row,col))%p for col in zip(*b)]
            for row in a]

def hs():
    # Real Moore graph control. All nontrivial roots of order 5 split over
    # the prime field F_11; using every primitive evaluation tests exact
    # recovery of the degree-<5 binary polynomial in this fixture.
    V=[(a,b,c) for a in range(2) for b in range(5) for c in range(5)]
    def move(v):
        a,b,c=v
        return (a,(b+1)%5,c) if a else (a,b,(c-b)%5)
    def edge(u,v):
        a,b,c=u
        d,e,f=v
        if a==d:
            return b==e and (f-c)%5 in (1+a,4-a)
        if a:
            return edge(v,u)
        return (f-c-b*e)%5==0
    O=[]
    seen=set()
    for v in V:
        if v in seen:
            continue
        if move(v)==v:
            continue
        row=[v]
        for _ in range(4):
            row.append(move(row[-1]))
        O.append(row)
        seen.update(row)
    assert len(O)==9
    F=[[tuple(d for d,y in enumerate(b) if edge(a[0],y))
        for b in O] for a in O]
    roots=[z for z in range(2,11) if pow(z,5,11)==1]
    choices={tuple(sum(pow(z,d,11) for d in S)%11 for z in roots):S
             for k in range(6) for S in combinations(range(5),k) if k<5}
    assert len(choices)==31
    blocks=[]
    for z in roots:
        P=[[(2*(i==j)-sum(pow(z,d,11) for d in S))*9%11
            for j,S in enumerate(row)] for i,row in enumerate(F)]
        assert mm(P,P,11)==P
        blocks.append(P)
    for r in range(9,0,-1):
        found=False
        for S in combinations(range(9),r):
            parts=[inverse([[P[i][j] for j in S] for i in S],11)
                   for P in blocks]
            if all(z is not None for z in parts):
                rebuilt=[mm(mm([[P[i][j] for j in S] for i in range(9)],Q,11),
                             [P[i] for i in S],11) for P,Q in zip(blocks,parts)]
                if r==4:
                    assert rebuilt==blocks
                found=True
                break
        if found:
            break
    assert all(P==Q for P,Q in zip(blocks,rebuilt))
    found = [[tuple((2*(i==j)-5*P[i][j])%11 for P in rebuilt)
              for j in range(9)] for i in range(9)]
    assert [[choices[z] for z in row] for row in found]==F
    return r,S

def tour(mask,n=8):
    a=[[False]*n for _ in range(n)]
    for b,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=bool(mask>>b&1)
        a[j][i]=not a[i][j]
    return a

def swap(a,i,j):
    b=[r[:] for r in a]
    b[i][j],b[j][i]=b[j][i],b[i][j]
    return b

def heights(a):
    n=len(a)
    d=[[0]*n for _ in range(1<<n)]
    out=[[j for j,v in enumerate(row) if v] for row in a]
    for i in range(n):
        d[1<<i][i]=1
    for m in range(1,1<<n):
        for i,v in enumerate(d[m]):
            if v:
                for j in out[i]:
                    if not m>>j&1:
                        d[m|1<<j][j]+=v
    return [1]+list(map(sum,d[1:]))

def shadows(a,b,guards,clone=()):
    hs,ht=heights(a),heights(b)
    F=(1<<len(a))-1
    ans=[]
    for k in range(len(guards)+1):
        for G in combinations(guards,k):
            remove=sum(1<<u for u in G)
            if clone and clone[0] in G:
                remove|=sum(1<<u for u in clone[1:])
            ans.append((G,hs[F^remove],ht[F^remove]))
    return ans

def inflate(a,x,m):
    n=len(a)
    orig=list(range(n))+[x]*(m-1)
    return [[(i<j if p==q==x else a[p][q]) for j,q in enumerate(orig)]
            for i,p in enumerate(orig)]

def seams():
    a=tour(113248304)
    b=swap(a,1,3)
    x=shadows(a,b,(0,4,7))
    assert x[0]==((),113,115) and all(s==t for _,s,t in x[1:])
    a=tour(75588141)
    b=swap(a,2,5)
    assert min(sum(a[1]),7-sum(a[1]))==1
    r=[]
    for m in (1,2,3,4):
        w=shadows(inflate(a,1,m),inflate(b,1,m),(3,4))
        assert w[0]==((),54*2**(m-1)-1,54*2**(m-1)-3)
        assert all(s==t for _,s,t in w[1:])
        r.append(w)
    return x,r[:2]

if __name__=='__main__':
    print(capacities())
    for n in range(1,11):
        for k in range(n+1):
            for S in combinations(range(n),k):
                assert unrank_subset(rank_subset(S,n),n)==S
    print('Hoffman–Singleton rank and basis:',hs())
    print('neutral deletions:',seams())
