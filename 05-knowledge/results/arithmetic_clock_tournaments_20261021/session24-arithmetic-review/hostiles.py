from itertools import product,combinations
from math import gcd
import json
checks=0
def need(x,s):
 global checks
 checks+=1
 if not x: raise RuntimeError(s)
def module(adj,M):
 return all(len({adj[z][v] for v in M})==1 for z in range(len(adj)) if z not in M)
def skew(n):
 for b in product((0,1),repeat=(n-1)//2):
  yield {i if v else n-i for i,v in zip(range(1,(n+1)//2),b)}
def mat(n,S):return [[(j-i)%n in S for j in range(n)] for i in range(n)]
# Every subset, not pair-closure, provides a separate module/subgroup test.
census=[]
for n in (3,5,7,9):
 decom=0
 for S in skew(n):
  adj=mat(n,S)
  allmods=[set(M) for k in range(2,n) for M in combinations(range(n),k) if module(adj,set(M))]
  subs=[set(range(0,n,d)) for d in range(2,n) if n%d==0]
  need(bool(allmods)==any(module(adj,H) for H in subs),'subgroup criterion')
  decom+=bool(allmods)
 census.append([n,2**((n-1)//2),decom])
# Necessity of strong hypothesis: overlapping maximal proper modules of transitive T3.
adj=[[j>i for j in range(3)] for i in range(3)]
need(module(adj,{0,1}) and module(adj,{1,2}),'non-strong overlap hostile')
# Every residue of exact conjugacy, plus wrong clock identification hostile.
need(all((3*(4*x+1)+1)%341 == 4*(3*x+1)%341 for x in range(341)),'CRT conjugacy')
need(all((3*((114*(x-1))%341)+1)%341==x for x in range(341)),'conjugacy inverse')
need(81%341!=4 and pow(81,5,341)!=1 and pow(4,5,341)==1,'golden clock is not trunk clock')
# Explicit ordinary subset signed minimality, with no complement operation.
for u,v in combinations(range(341),2):
 need({(u+v-x)%341 for x in (u,v)}=={u,v},'pair reflection')
# Reconstruct full profiles independently as triple loops on CRT coordinates.
profiles=[]
for p,q in ((11,31),(31,11)):
 R={a*a%p for a in range(1,p)};Q={b*b%q for b in range(1,q)}
 pts=list(product(range(p),range(q)))
 for delta in (1,-1):
  def s(a,b):return (1 if b in Q else -1) if a==0 else (1 if a in R else -1)*(delta if b==0 else 1)
  adj=[[i!=j and s((c-a)%p,(d-b)%q)==1 for j,(c,d) in enumerate(pts)] for i,(a,b) in enumerate(pts)]
  vals=[sum(adj[0][k] and adj[j][k] for k in range(p*q)) for j in range(1,p*q)]
  need({pts[j] for j,v in enumerate(vals,1) if v==max(vals)}=={(0,b) for b in range(1,q)},'intrinsic fibres')
  from collections import Counter
  profile=Counter(vals[j-1] for j in range(1,p*q) if adj[0][j]);profiles.append([p,delta,dict(sorted(profile.items()))])
  # Changing just one block's b translation preserves its internal tournament but destroys matching.
  perm=[pts.index((a,(b+1)%q if a==0 else b)) for a,b in pts]
  isauto=all(adj[i][j]==adj[perm[i]][perm[j]] for i in range(p*q) for j in range(p*q))
  need(isauto==(delta==1),'independent block permutation hostile')
# Both factors decomposable: lexicographic C3[C3] and C5[C5], orders9 and25.
# Test every rectangular subgroup directly, without inferring from primality of factors.
m,n=9,25
A={x for x in range(1,m) if (x%3==1 if x%3 else x//3==1)}
B={x for x in range(1,n) if (x%5 in {1,2} if x%5 else x//5 in {1,2})}
pts=list(product(range(m),range(n)))
def edge(u,v):
 a,b=(v[0]-u[0])%m,(v[1]-u[1])%n
 return (b in B) if a==0 else ((a not in A) if b==0 else a in A)
adj=[[u!=v and edge(u,v) for v in pts] for u in pts]
rects=0
for d in (1,3,9):
 for e in (1,5,25):
  H={i for i,(a,b) in enumerate(pts) if a%d==0 and b%e==0}
  if 1<len(H)<m*n:
   need(not module(adj,H),'composite-factor rectangular subgroup hostile');rects+=1
print(json.dumps({'checks':checks,'all_subset_cyclic_census':census,'independent_CRT_profiles':profiles,'composite_9x25_proper_rectangles':rects},indent=2,sort_keys=True))
