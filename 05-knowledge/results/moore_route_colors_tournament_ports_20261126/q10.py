from math import isqrt
from collections import defaultdict
def sd(x):
 for r in range(isqrt(x),0,-1):
  if x%(r*r)==0:return x//(r*r),r
typ=defaultdict(list)
for n in range(1,400):
 d,r=sd(n*(400-n))
 typ[d].extend([(n,-r),(n,r)])
found=[]
for D,ss in sorted(typ.items()):
 vs=sorted(ss+[(400,0)])
 def val(x,y):
  n,r=x; m,s=y
  q,d=divmod(n*m+D*r*s,50*n)
  if d or not 0<=q<=min(57,m):return None
  return q
 adj=[[val(x,y) is not None and val(y,x) is not None for y in vs] for x in vs]
 def go(ch,ii,left,total,bal):
  if not left:
   if not total and not bal:found.append((D,tuple(vs[j] for j in ch)))
   return
  for ix in range(ii,len(vs)):
   n,r=vs[ix]
   if n*left>total:break
   if total-n>400*(left-1):continue
   if not all(adj[ix][j] for j in ch):continue
   go(ch+[ix],ix,left-1,total-n,bal+r)
 go([],0,10,3250,0)
print("candidate color profiles",len(found))
for D,rows in found[:100]:
 print(D,rows)

# Full matrix verification, not just the eigenvalue-trace condition.
invariants=set()
for D,r in found:
 n=[x[0] for x in r];t=[x[1] for x in r]
 b=[[0 if i==j else (n[i]*n[j]+D*t[i]*t[j])//(50*n[i])
     for j in range(10)] for i in range(10)]
 assert all(sum(row)==57 for row in b)
 for i in range(10):
  for j in range(10):
   assert n[i]*b[i][j]==n[j]*b[j][i]
   assert sum(b[i][h]*b[h][j] for h in range(10))+b[i][j]-(56 if i==j else 0)==n[j]
 # Normal form under global sign reversal (order within same-size classes irrelevant).
 invariants.add((D,min(r,tuple(sorted((a,-c) for a,c in r)))))
print("up to global sign reversal",len(invariants))
for D,r in sorted(invariants):
 print("types (size, sign label):", D,r)
