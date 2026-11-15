"""Exact certificate of both physical residue-fiber routes."""
from itertools import combinations, permutations
from math import comb
def arcs(n,g):
    return {(i,j) if g>>k&1 else (j,i)
            for k,(i,j) in enumerate(combinations(range(n),2))}
def strong(n,A):
    for j in range(n):
        R={j}
        while True:
            S=R|{y for x,y in A if x in R}
            if S==R: break
            R=S
        if len(R)!=n:return False
    return True
def h(n,A):
    r=[{} for _ in range(1<<n)]
    for i in range(n):r[1<<i][i]=1
    for s in range(1,1<<n):
        for i,v in r[s].items():
            for j in range(n):
                if not s>>j&1 and (i,j) in A:
                    u=r[s|(1<<j)]
                    u[j]=u.get(j,0)+v
    return sum(r[-1].values())
def response(n,g):
    A=arcs(n,g);v=n-1
    out=[0]*(n+1)
    for p in permutations(range(v)):
        a=((v,p[0]) in A)+((p[-1],v) in A); b=0
        for i,j in zip(p,p[1:]):
            d=((i,v) in A and (v,j) in A)
            if (i,j) in A:a+=d
            elif d:b+=1
            else:break
        else:
            for j in range(a+1):out[b+j]+=comb(a,j)
    return out
def inflate(n,g,d,c):
    v=n-1; B=arcs(n,g); A={(a+v,b+v) for a,b in arcs(d,c)}
    for i,j in B:
        if v not in (i,j):A.add((i,j))
        elif i==v:A.update((x,j) for x in range(v,v+d))
        else:A.update((i,x) for x in range(v,v+d))
    return A
def slots(d,A,j):
    from itertools import product
    s=0
    for w in product(range(j),repeat=d):
        v=1
        for i in range(j):
            S=[k for k in range(d) if w[k]==i]
            v*=sum(all((a,b) in A for a,b in zip(p,p[1:])) for p in permutations(S))
        s+=v
    return s
def need(v):
    if not v:raise AssertionError("route witness failed")
for n,T,P in ((7,[366909,301373,268605],[(15,27,9),(15,29,10),(17,31,11)]),
              (6,[29130,29002,29000],[(9,27,21,4),(11,29,21,4),(13,31,21,4)])):
    for g,p in zip(T,P):
        need(response(n,g)==list(p)+[0]*(n+1-len(p)))
        need(strong(n,arcs(n,g)))
    for g,k in zip(T,T[1:]):
        q=g^k;need(q.bit_count()==1)
        need(n-1 not in list(combinations(range(n),2))[q.bit_length()-1])
    for d,c in ((1,0),(2,1),(3,5),(3,7),(4,63)):
        A=arcs(d,c); H=h(d,A); q2=slots(d,A,2);q3=slots(d,A,3)
        expected=[9*H+j*q2 for j in [9,10,11]] if n==7 else [4*q3+9*q2+j*H for j in [-3,-1,1]]
        got=[h(n+d-1,inflate(n,g,d,c)) for g in T]
        need(got==expected)
print("Both port triples, all four flips, and thirty literal inflated counts verified.")

from itertools import combinations
from functools import lru_cache
def resp(n,g,D,E):
 es=list(combinations(range(n),2)); a={x if (g>>k)&1 else x[::-1] for k,x in enumerate(es)}
 full=((1<<n)-1)^(1<<D)^(1<<E)
 @lru_cache(None)
 def f(mask,last,p,q):
  out={(p,q):1} if mask==full and p and q else {}
  for j in range(n):
   if (last,j) not in a and last!=-1:continue
   if j==D:
    if p==n-1:continue
    h=f(mask,j,p+1,q)
   elif j==E:
    if q==n-1:continue
    h=f(mask,j,p,q+1)
   else:
    if mask>>j&1:continue
    h=f(mask|1<<j,j,p,q)
   for k,v in h.items():out[k]=out.get(k,0)+v
  return out
 return f(0,-1,0,0)

for e,coeffs in enumerate([(0,1,0,0),(3,1,0,0),(-5,-2,-2,-2)]):
    a,b=(resp(7,g,6,e) for g in [366909,301373])
    S=[(1,2),(2,2),(1,3),(2,3)]
    need(tuple(b.get(x,0)-a.get(x,0) for x in S)==coeffs)
print("Two-port separation failures verified at all three available ports.")

for e,vs in enumerate([(2,0),(2,1),(4,2)]):
    X,Y,Z=(resp(7,t,6,e) for t in (366909,301373,268605))
    key=X.keys()|Y.keys()|Z.keys()
    b={k:Z.get(k,0)-2*Y.get(k,0)+X.get(k,0) for k in key}
    b={k:v for k,v in b.items() if v}
    need(b=={k:v for k,v in zip(((1,2),(2,2)),vs) if v})
print("Exact strict two-port step defects verified.")
