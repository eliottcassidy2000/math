from itertools import permutations
from collections import Counter,defaultdict
from math import log,factorial,isclose

def code(A,p):return tuple(A[p[i]][p[j]] for i in range(len(A)) for j in range(i+1,len(A)))
def canon(A):return min(code(A,p) for p in permutations(range(len(A))))
def rev(A):return tuple(tuple(0 if i==j else 1-A[i][j] for j in range(len(A))) for i in range(len(A)))
def merged(A):return min(canon(A),canon(rev(A)))
def tilings(n):
 pairs=[(i,j) for i in range(n) for j in range(i+2,n)]
 for z in range(1<<len(pairs)):
  A=[[0]*n for _ in range(n)]
  for i in range(n-1):A[i][i+1]=1
  for k,(i,j) in enumerate(pairs):A[i][j]=(z>>k)&1;A[j][i]=1-A[i][j]
  yield tuple(map(tuple,A))
def ends(A):
 n=len(A);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for s in range(1,1<<n):
  for i in range(n):
   if not(s>>i&1):continue
   for j in range(n):
    if not(s>>j&1) and A[i][j]:D[s|1<<j][j]+=D[s][i]
 return D[-1]
def aut(A):return sum(code(A,p)==code(A,range(len(A))) for p in permutations(range(len(A))))
def delete(A,v):return tuple(tuple(A[i][j] for j in range(len(A)) if j!=v) for i in range(len(A)) if i!=v)
prev=None
for n in range(1,6):
 ts=list(tilings(n)); M=Counter(map(merged,ts)); reps={merged(A):A for A in ts}
 for u,A in reps.items():
  factor=1 if canon(A)==canon(rev(A)) else 2
  assert M[u]==factor*sum(ends(A))//aut(A)
 entropy=-sum((v/len(ts))*log(v/len(ts)) for v in M.values())
 assert isclose(entropy+sum(v*log(v)/len(ts) for v in M.values()),log(len(ts)),abs_tol=1e-12)
 if prev:
  E=Counter((merged(delete(A,n-1)),merged(A)) for A in ts)
  assert Counter({u:sum(v for (a,b),v in E.items() if a==u) for u in prev})==Counter({u:v*2**(n-2) for u,v in prev.items()})
  assert Counter({u:sum(v for (a,b),v in E.items() if b==u) for u in M})==M
  for d,A in reps.items():
   numerator=Counter()
   AA=[A] if canon(A)==canon(rev(A)) else [A,rev(A)]
   for B in AA:
    for v,e in enumerate(ends(B)):numerator[merged(delete(B,v))]+=e
   for c,v in numerator.items():assert v==aut(A)*E[c,d]
  # Actual extensions of individual tilings, not merely aggregated node rows.
  profiles=defaultdict(set)
  for B in tilings(n-1):
   prof=Counter()
   for mask in range(1<<(n-2)):
    A=[list(row)+[0] for row in B]+[[0]*n]
    A[n-2][n-1]=1
    for i in range(n-2):A[i][n-1]=mask>>i&1;A[n-1][i]=1-A[i][n-1]
    prof[merged(A)]+=1
   profiles[merged(B)].add(tuple(sorted(prof.items())))
  bad=sum(len(p)>1 for p in profiles.values())
  print('order',n,'tilings',len(ts),'merged nodes',len(M),'nonlumpable parent fibres',bad,'endpoint transport PASS')
 else:print('order',n,'tilings',len(ts),'merged nodes',len(M))
 prev=M
print('PASS all merged weights, endpoint deletion, row/column transport, entropy, through order 5')
# Explicit nonlumpability hostile, same strong class, different terminal presentation.
def children_hist(B):
 n=len(B)+1;hist=Counter()
 for mask in range(1<<(n-2)):
  A=[list(row)+[0] for row in B]+[[0]*n];A[n-2][n-1]=1
  for i in range(n-2):A[i][n-1]=mask>>i&1;A[n-1][i]=1-A[i][n-1]
  hist[sum(ends(A))]+=1
 return dict(sorted(hist.items()))
B0,B1=list(tilings(4))[:2]
assert canon(B0)==canon(B1)
assert children_hist(B0)=={5:1,9:4,11:2,13:1}
assert children_hist(B1)=={5:1,9:4,11:1,13:1,15:1}
print('same H5 source class, two last-path presentations:',children_hist(B0),children_hist(B1))
