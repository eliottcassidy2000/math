from itertools import combinations,permutations,product
from collections import Counter
from math import prod
def compile(n,arcs):
 out=Counter()
 for d in range(len(arcs)+1):
  for D in permutations(arcs,d):
   forced={};lo=[-1]*(d+1);hi=[n]*(d+1)
   for j,(a,b) in enumerate(D):
    if (a in forced and forced[a]!=j+1) or (b in forced and forced[b]!=j):break
    forced[a]=j+1;forced[b]=j
    lo[j+1]=a;hi[j]=b
   else:
    if any(l>h for l,h in zip(lo,hi)):continue
    U=set(arcs)-set(D)
    T=[(forced[i],) if i in forced else tuple(j for j in range(d+1) if lo[j]<=i<=hi[j]) for i in range(n)]
    for A in product(*T):
     if any(not lo[j]<=i<=hi[j] for i,j in enumerate(A)):continue
     forbidden=[(a,b,A[a]) for a,b in U if A[a]==A[b] and not any(A[k]==A[a] for k in range(a+1,b))]
     choices=[set(j for j in range(d+1) if lo[j]<=p<hi[j]) for p in range(n-1)]
     for mask in range(1<<len(forbidden)):
      V=list(map(set,choices))
      for k,(a,b,j) in enumerate(forbidden):
       if mask>>k&1:
        for p in range(a,b):V[p].discard(j)
      out[tuple(map(len,V))]+=(-1)**mask.bit_count()
 return {k:c for k,c in out.items() if c}
def value(P,g):return sum(c*prod(b**x for b,x in zip(a,g)) for a,c in P.items())
def literal(g,arcs):
 n=sum(g)+len(g)+1;A=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
 mark=[sum(g[:i])+i for i in range(len(g)+1)]
 for a,b in arcs:
  u,v=mark[a],mark[b];A[u]^=1<<v;A[v]^=1<<u
 H=[[0]*n for _ in range(1<<n)]
 for i in range(n):H[1<<i][i]=1
 for S in range(1,1<<n):
  for i,x in enumerate(H[S]):
   if not x:continue
   e=A[i]&~S
   while e:
    w=e&-e;e-=w
    H[S|w][w.bit_length()-1]+=x
 return sum(H[-1])
if __name__=='__main__':
 for a in (((0,1),(2,3)),((0,3),(1,2)),((0,2),(1,3))):
  F=compile(4,a);print(a,sorted(F.items()))
  for v in product(range(3),repeat=3):assert value(F,v)==literal(v,a),(a,v)
 for a in (((0,1),(1,2)),((0,2),(1,2)),((0,1),(0,2))):
  F=compile(3,a);print(a,sorted(F.items()))
  for v in product(range(4),repeat=2):assert value(F,v)==literal(v,a),(a,v)
