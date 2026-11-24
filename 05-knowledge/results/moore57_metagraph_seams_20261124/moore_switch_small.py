from itertools import combinations
import numpy as np
for k in [2,3]:
 V=list(range(5)) if k==2 else list(combinations(range(5),2))
 n=len(V)
 A=np.array([[int((v-u)%5 in (1,4)) if k==2 else int(not set(u)&set(v)) for v in V]for u in V])
 hits={}
 for x in (2,4,6):
  for y in (2,4,6):
   if x>y or x+y>n:continue
   out=[]
   for X in combinations(range(n),x):
    for Y in combinations([i for i in range(n) if i not in X],y):
     if x==y and X>Y:continue
     B=A[np.ix_(X,Y)]
     if not ((B.sum(0)==x//2).all() and (B.sum(1)==y//2).all()): continue
     C=A.copy(); C[np.ix_(X,Y)]=1-B;C[np.ix_(Y,X)]=1-B.T
     if np.array_equal(C@C+C,(k-1)*np.eye(n,dtype=int)+1):out.append((X,Y))
   if out:hits[(x,y)]=(len(out),out[0])
 print(k,V,hits)
