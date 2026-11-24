import numpy as np
from moore_grid_pipeline import moore, grid, target
from itertools import product,combinations

def two_rows(C,d):
 s=d-2;n=s*d
 if s<3:return
 M=target(s,d)-C@C-C
 if not ((M.diagonal()==2).all() and ((M-np.diag(M.diagonal())<=1)&(M>=0)).all()):return
 R=(M==1)
 T=sorted({tuple(sorted((int(u),int(v),*map(int,np.flatnonzero(R[u]&R[v])))))
           for u,v in zip(*np.where(np.triu(R)))})
 if len(T)!=2*d or any(len(t)!=s or len({x//d for x in t})!=s or len({x%d for x in t})!=s for t in T):return
 H=np.zeros((n,2*d),dtype=np.int64)
 for i,t in enumerate(T):H[list(t),i]=1
 if not np.array_equal(H@H.T,M):return
 L=H.T@H-s*np.eye(2*d,dtype=np.int64)
 if not ((L==0)|(L==1)).all() or not (L.sum(1)==s).all():return
 sides=np.full(2*d,-1);sides[0]=0;queue=[0]
 for v in queue:
  for w in np.flatnonzero(L[v]):
   if sides[w]==sides[v]:return
   if sides[w]<0:sides[w]=1-sides[v];queue.append(w)
 if list(np.bincount(sides+1))!=[0,d,d]:return
 opts=[tuple(sorted(set(range(d))-{x%d for x in t})) for t in T]
 G=[[] for _ in range(4*d)]
 def forbid(i,a,j,b):
  G[2*i+a].append(2*j+1-b);G[2*j+b].append(2*i+1-a)
 for i,j in combinations(range(2*d),2):
  if sides[i]==sides[j] or L[i,j]:
   for a,b in product(range(2),repeat=2):
    if opts[i][a]==opts[j][b]:forbid(i,a,j,b)
 for x in range(n):
  i,j=np.flatnonzero(H[x])
  miss=set(range(d))-{x%d,*[y%d for y in np.flatnonzero(C[x])]}
  for a,b in product(range(2),repeat=2):
   if {opts[i][a],opts[j][b]}!=miss:forbid(i,a,j,b)
 K=(C+np.eye(n,dtype=np.int64))@H
 for i in range(2*d):
  j,k=np.flatnonzero((sides!=sides[i])&(L[i]==0))
  for a in range(2):
   for b in range(2):
    if opts[i][a]==opts[j][b] and not np.array_equal(K[:,i]+H[:,k],np.arange(n)%d!=opts[i][a]):forbid(i,a,j,b)
    if opts[i][a]==opts[k][b] and not np.array_equal(K[:,i]+H[:,j],np.arange(n)%d!=opts[i][a]):forbid(i,a,k,b)
 # reachability-based constructive 2-SAT, portable replay
 assign=set()
 def closure(start):
  seen={start};q=[start]
  for u in q:
   for v in G[u]:
    if v not in seen:seen.add(v);q.append(v)
  return seen
 for i in range(2*d):
  for a in range(2):
   h=2*i+a
   if h in assign:break
   q=assign|closure(h)
   if not any((v^1) in q for v in q):assign=q;break
  else:return
 perm=sorted(range(2*d),key=lambda i:(sides[i],opts[i][2*i+1 in assign]))
 D=np.zeros((2*d,2*d),dtype=np.int64)
 for i,j in combinations(range(2*d),2):
  if sides[i]!=sides[j] and L[i,j]==0 and opts[i][2*i+1 in assign]!=opts[j][2*j+1 in assign]:D[i,j]=D[j,i]=1
 H=H[:,perm];D=D[np.ix_(perm,perm)]
 B=np.block([[C,H],[H.T,D]])
 assert np.array_equal(B@B+B,target(d,d))
 return B

if __name__=='__main__':
 A=moore(7); d=6;n=len(A);count=0
 for r in range(n):
  for a in np.flatnonzero(A[r]):
   B=grid(A,r,a)
   assert two_rows(B[:-2*d,:-2*d],d) is not None;count+=1
 print('Hoffman–Singleton two-row completions:',count)
