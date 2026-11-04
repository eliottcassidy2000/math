#!/usr/bin/env python3
"""Independent literal ordered-cover and insertion audit; no author imports."""
import itertools,json

def check(p,label):
 if not p:raise RuntimeError(label)

def h(a):
 return sum(all(a[p[i]][p[i+1]] for i in range(len(p)-1)) for p in itertools.permutations(range(len(a))))

def covers(a):
 n=len(a);K=[[0]*n for _ in a]
 # Each permutation and split is one ordered two-path cover, ignoring split edge.
 for p in itertools.permutations(range(n)):
  for s in range(1,n):
   if all(a[p[i]][p[i+1]] for i in range(n-1) if i!=s-1):K[p[s-1]][p[s]]+=1
 return K

def ext(a,b):
 n=len(a);r=[row+[not (b>>i&1)] for i,row in enumerate(a)];r.append([bool(b>>i&1) for i in range(n)]+[False]);return r

tournaments=responses=submod=0
for n in range(1,5):
 pairs=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  a=[[False]*n for _ in range(n)]
  for q,(i,j) in enumerate(pairs):a[i][j]=bool(mask>>q&1);a[j][i]=not a[i][j]
  H=h(a);K=covers(a);W={(i,j):(K[i][j]+K[j][i])//2 for i,j in pairs}
  for i,j in pairs:
   check(K[i][j]+K[j][i]==2*W[i,j]>0,'positive integer')
   a[i][j],a[j][i]=a[j][i],a[i][j];other=covers(a);a[i][j],a[j][i]=a[j][i],a[i][j]
   check((K[i][j],K[j][i])==(other[i][j],other[j][i]),'mutual-edge independence')
  F=[h(ext(a,b)) for b in range(1<<n)]
  for b in range(1<<n):
   cut=sum(W[min(i,j),max(i,j)] for i in range(n) for j in range(n) if a[i][j] and b>>i&1 and not b>>j&1)
   check(F[b]==H+2*cut,'literal cut');responses+=1
  for A in range(1<<n):
   for B in range(1<<n):
    cross=sum(W[min(i,j),max(i,j)] for i in range(n) for j in range(n) if A>>i&1 and not B>>i&1 and B>>j&1 and not A>>j&1)
    check(F[A]+F[B]-F[A|B]-F[A&B]==2*cross,'strict submod');submod+=1
  tournaments+=1
print(json.dumps(dict(status='PASS',literal_tournaments=tournaments,literal_extensions=responses,submodularity_pairs=submod,independent_method='permutations and split positions; no imports from author code'),indent=2))
