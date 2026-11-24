from collections import Counter
from itertools import combinations, permutations, product
from math import factorial

def word_h(w):
    if not w: return 1
    ans=w[0]+1-w[-1]
    b=0
    for i,x in enumerate(w):
        if not x: ans+=b
        if i and not w[i-1] and x: ans+=1
        b=2*b+x
    return ans

def value(w):
    b=0
    for x in w:b=2*b+x
    return b

def mate_h(w):
    assert w and w[-1]==0
    empty=(0,0,0,1) # clipped run count, top run singleton, pending bad descent, endpoint permitted
    D={(empty,empty,2):1}
    for i,a in enumerate(w):
        out=1-a if i+1<len(w) else 0 # flipped apex outcontact a -> i
        N=Counter()
        for (l,r,last), num in D.items():
            for j,old in enumerate((l,r)):
                c,s,b,ok=old
                # For left lane keep the low run's last vertex, for right the high run's first.
                good=(1-out) if j==0 else out
                z=[]
                if not c:
                    z.append((1,1,0,good))
                else:
                    if j==last:
                        z.append((c,0,0,good if j==0 and c==1 else ok))
                    if not b:
                        z.append((2,1,int(j==last and s),good if j==1 else ok))
                for new in z:
                    N[(new,r,j) if j==0 else (l,new,j)]+=num
        D=N
    return sum(n for (l,r,j),n in D.items() if not l[2] and not r[2] and l[3] and r[3])

from itertools import combinations, product

def hp(A):
    n=len(A); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for k in range(1,1<<n):
        for i in range(n):
            if dp[k][i]:
                for j in range(n):
                    if not (k>>j)&1 and A[i][j]:
                        dp[k|1<<j][j]+=dp[k][i]
    return sum(dp[-1])

def marked(w,flip):
    n=len(w)+1
    A=[[int(i<j) for j in range(n)] for i in range(n)]
    for i,x in enumerate(w): A[-1][i]=x; A[i][-1]=1-x
    if flip:
        for i,j in combinations(range(n),2):
            if j>i+1: A[i][j],A[j][i]=A[j][i],A[i][j]
    return A

for m in range(2,9):
    for w in ((1,)*m+(0,), (1,)+(0,)*m):
        assert hp(marked(w,False))==2**m+1==word_h(w)
        assert hp(marked(w,True))==mate_h(w)
m=55
b=(1,)*m+(0,); c=(1,)+(0,)*m
hb,hc=mate_h(b),mate_h(c)
for w, expected_selfconverse in [(b,False),(c,True)]:
    deg=Counter(map(sum,marked(w,True)))
    assert max(deg.values())==2
    assert (deg==Counter({m+1-k:v for k,v in deg.items()}))==expected_selfconverse
print("star source and both mates checked by subset DP for 2..8 leaves")
print("56-vertex star plus helper: source",2**m+1,"black mate",hb,"blue mate",hc)
print("merged mask counts:",2**m+1,2*hb,hc)

a_values=[a for a in range(0,351,5) if (335+a)%15==0 and (1729-(335+a)//15)%5==0]
for a in a_values:
    z=2*a//5
    print("order-5 graph fixed HS: a =",a,"circulant pentagons",z,
          "alternating ring vertices",140-z,"E_7 nontrivial multiplicity",(1729-(335+a)//15)//5)
assert a_values==[25,100,175,250,325]

# Positive test on the actual Hoffman–Singleton graph, with fixed pentagon,
# using the same group-ring certificate at k=7, p=5.
G=[set() for _ in range(50)]
def add(G,u,v): G[u].add(v); G[v].add(u)
for h,x in product(range(5),repeat=2):
    add(G,5*h+x,5*h+(x+1)%5)
    add(G,25+5*h+x,25+5*h+(x+2)%5)
    for i in range(5): add(G,5*h+x,25+5*i+(x+h*i)%5)
def t(v):
    h,x=divmod(v%25,5)
    return 5*((h+1)%5)+x if v<25 else 25+5*h+(x+h)%5
assert all({t(w) for w in G[v]}==G[t(v)] for v in range(50))
F=[v for v in range(50) if t(v)==v]
orb=[]; seen=set(F)
for v in range(50):
    if v in seen: continue
    row=[v]
    for _ in range(4): row.append(t(row[-1]))
    assert len(set(row))==5 and t(row[-1])==v
    seen.update(row);orb.append(row)
rep=[x[0] for x in orb]
M=[[{d for d in range(5) if o[d] in G[i]} for o in orb] for i in rep]
S=[[int(i in G[f]) for i in rep] for f in F]
for i in range(len(orb)):
    for j in range(len(orb)):
        h=[int(u in M[i][j]) for u in range(5)]
        for k in range(len(orb)):
            for a,b in product(M[i][k],M[k][j]):h[(a+b)%5]+=1
        c=sum(row[i]*row[j] for row in S)
        assert h==[1-c+6*int(i==j and u==0) for u in range(5)]
for f,x in enumerate(F):
    for j in range(len(orb)):
        assert S[f][j]+sum(S[i][j] for i,y in enumerate(F) if y in G[x])+sum(S[f][i]*len(M[i][j]) for i in range(len(orb)))==1
    for g,y in enumerate(F):
        assert int(y in G[x])+len(G[x]&G[y]&set(F))+5*sum(a*b for a,b in zip(S[f],S[g]))==1+6*int(f==g)
print("group-ring certificate passed on real 50-site HS lift:",len(F),"fixed,",len(orb),"free orbits")

# Lift alternating doubled rings with arbitrary center flux. The obstruction
# for length 4 and the sharp local survival of lengths 6,8,10 are checked.
for n in (4,6,8,10):
    for flux in range(5):
        A=[set() for _ in range(5*n)]
        for i in range(n):
            p=1+i%2; m=flux if i==n-1 else 0
            for x in range(5):
                for s in (-p,p):
                    add(A,5*i+x,5*((i+1)%n)+(x+m+s)%5)
        assert not any(A[i]&A[j] for i in range(len(A)) for j in A[i])
        assert any(len(A[i]&A[j])>=2 for i,j in combinations(range(len(A)),2))==(n==4)
print("all five center fluxes: rings 4 rejected, rings 6,8,10 locally pass girth")
