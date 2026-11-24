import numpy as np
from itertools import combinations
from scipy.optimize import linprog
from moore_grid_pipeline import moore

def encode(A,U,markers):
 U=sorted(U);O=sorted(set(range(len(A)))-set(U));m=len(O);s=len(U)
 j,x=next((j,x) for j in O for x in U if A[j,x])
 O.remove(j);O.append(j);U.remove(x);U.append(x)
 labels=O+U
 W=list(reversed(range(m+s,2*m+s)))
 X=list(range(m));Y=list(range(m,m+s))
 n=2*m+s
 T=np.zeros((n+4*markers,n+4*markers),dtype=int)
 def edge(i,j,p=True):T[i,j]=int(p);T[j,i]=int(not p)
 for L in [X,Y,W[::-1]]:
  for i,j in combinations(L,2):edge(i,j)
 for i,v in enumerate(O):
  for j,z in enumerate(O):edge(X[i],W[j],v==z or A[v,z])
  for j,z in enumerate(U):
   edge(Y[j],W[i],A[z,v])
   edge(X[i],Y[s-1-j],A[z,v])
 spine=list(range(n))
 if markers:
  a,b,alpha,beta=range(n,n+4)
  for L,l,r in [(X,1,0),(Y,1,1),(W,0,1)]:
   for v in L:edge(a,v,l);edge(v,b,r);edge(alpha,v);edge(v,beta)
  for i,j in [(a,b),(a,alpha),(a,beta),(alpha,b),(beta,b),(beta,alpha)]:edge(i,j)
  spine=[a,alpha,*spine,beta,b]
 assert all(T[i,j] for i,j in zip(spine,spine[1:]))
 R=T[np.ix_(spine,spine)]
 assert np.array_equal(R,R[::-1,::-1].T)
 return T,labels

def decode(T,d):
 s=T.sum(1)
 alpha,=np.flatnonzero(s==len(T)-3)
 beta,=np.flatnonzero(s==2)
 a,=(i for i in np.flatnonzero(T[:,alpha]) if i!=beta)
 b,=(i for i in np.flatnonzero(T[beta]) if i!=alpha)
 kinds={(1,0):[],(1,1):[],(0,1):[]}
 for v in set(range(len(T)))-{a,b,alpha,beta}: kinds[T[a,v],T[v,b]].append(v)
 def ordered(L):return sorted(L,key=lambda v:-T[v,L].sum())
 X,Y,W=map(ordered,kinds.values());W=W[::-1]
 D=T[np.ix_(X,W)]-np.eye(len(X),dtype=int)
 B=T[np.ix_(Y,W)].T
 q=len(Y); expected=d-B.sum(0)
 P=list(combinations(range(q),2))
 # Polynomial rational feasibility decoder. Here we use HiGHS then check its
 # integral output exactly, justified by the independent uniqueness proof.
 E=np.zeros(((len(X)+1)*q,len(P)),dtype=int)
 for z,(i,j) in enumerate(P):
  E[0*q+i,z]=E[0*q+j,z]=1
  E[q+j::q,z]=B[:,i];E[q+i::q,z]=B[:,j]
 rhs=np.r_[expected,(1-D@B-B).ravel()]
 sol=linprog(np.zeros(len(P)),A_eq=E,b_eq=rhs,bounds=(0,None),method='highs')
 assert sol.success
 c=np.rint(sol.x).astype(int); assert (c>=0).all() and (E@c==rhs).all()
 C=np.zeros((q,q),dtype=int)
 for (i,j),w in zip(P,c):C[i,j]=C[j,i]=w
 return np.block([[D,B],[B.T,C]])

if __name__=='__main__':
 for d in [2,3,7]:
  A=moore(d)
  if d==7:
   U=list(range(25));T,labels=encode(A,U,True)
   np.random.default_rng(0).shuffle(perm:=np.arange(len(T)))
   assert np.array_equal(decode(T[np.ix_(perm,perm)],d),A[np.ix_(labels,labels)])
   print('unmarked blue HS:',len(T))
  else:
   U=[0,*np.flatnonzero(A[0])];T,labels=encode(A,U,False)
   m=len(labels)-len(U);B=T[m:len(labels),len(labels):][:,::-1].T
   r,=np.flatnonzero(d-B.sum(0)==d)
   C=np.zeros((len(U),len(U)),dtype=int);C[r,:]=C[:,r]=1;C[r,r]=0
   assert np.array_equal(C,A[np.ix_(labels[m:],labels[m:])])
   print('fixed-spine blue degree',d,':',len(T))
