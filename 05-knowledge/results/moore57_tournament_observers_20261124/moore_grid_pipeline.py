import numpy as np

def moore(d):
 if d==2:
  V=range(5)
  return np.array([[(i-j)%5 in (1,4) for j in V] for i in V], dtype=np.int64)
 if d==3:
  from itertools import combinations
  V=list(map(set,combinations(range(5),2)))
  return np.array([[not i&j for j in V] for i in V],dtype=np.int64)
 assert d==7
 V=[(side,i,x) for side in range(2) for i in range(5) for x in range(5)]
 def edge(v,w):
  a,i,x=v;b,j,y=w
  if a==b:return i==j and (x-y)%5 in (1+a,4-a)
  if a:return edge(w,v)
  return (y-x-i*j)%5==0
 return np.array([[edge(v,w) for w in V] for v in V],dtype=np.int64)

def grid(A, root, anchor):
 k=int(A[root].sum());d=k-1
 arms=np.flatnonzero(A[root]);assert anchor in arms
 ids=[i for i in arms if i!=anchor]
 labels=[x for x in np.flatnonzero(A[anchor]) if x!=root]
 X=np.array([[next(x for x in np.flatnonzero(A[i]) if A[l,x]) for l in labels] for i in ids])
 B=A[np.ix_(X.ravel(),X.ravel())]
 assert np.array_equal(B@B+B,target(d,d))
 return B

def target(s,d):
 return (d-1)*np.eye(s*d,dtype=np.int64)+np.kron(1-np.eye(s,dtype=np.int64),1-np.eye(d,dtype=np.int64))

def inverse_last(C,d):
 s=d-1;n=s*d
 assert C.shape==(n,n)
 M=target(s,d)-C@C-C
 if not ((M.diagonal()==1).all() and ((M==0)|(M==1)).all()):return None
 label=np.full(n,-1)
 for v in range(n):
  if label[v]!=-1:continue
  clique=np.flatnonzero(M[v]);parts={x//d for x in clique};cols={x%d for x in clique}
  if len(clique)!=s or len(parts)!=s or len(cols)!=s:return None
  if not np.array_equal(M[clique],np.broadcast_to(M[v],(s,n))):return None
  miss,=set(range(d))-cols
  if miss in label:return None
  label[clique]=miss
 for v,l in enumerate(label):
  seen=[v%d,*[x%d for x in np.flatnonzero(C[v])]]
  if len(set(seen))!=s or int(l) in seen:return None
 H=np.eye(d,dtype=np.int64)[label]
 B=np.block([[C,H],[H.T,np.zeros((d,d),dtype=np.int64)]])
 assert np.array_equal(B@B+B,target(d,d))
 return B

if __name__=='__main__':
 for degree in [2,3,7]:
  A=moore(degree); n=len(A)
  assert np.array_equal(A@A+A,(degree-1)*np.eye(n,dtype=np.int64)+1)
  if degree==2: continue
  d=degree-1
  count=0
  for r in range(n):
   for a in np.flatnonzero(A[r]):
    B=grid(A,r,a)
    assert np.array_equal(inverse_last(B[:-d,:-d],d),B);count+=1
  print('degree',degree,'exact inversions',count,flush=True)
