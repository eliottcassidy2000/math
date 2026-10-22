import numpy as np
from itertools import combinations,product

def rootsA(n):
 r=[]
 for i in range(n+1):
  for j in range(n+1):
   if i!=j:
    x=np.zeros(n+1,dtype=int);x[i]=1;x[j]=-1;r.append(x)
 return np.array(r)
def rootsD(n):
 r=[]
 for i,j in combinations(range(n),2):
  for a,b in product([-1,1],repeat=2):
   x=np.zeros(n,dtype=int);x[i]=a;x[j]=b;r.append(x)
 return np.array(r)
def count(r,N=144):
 G=r@r.T;c=0;ext=0
 for i in range(len(r)):
  for j in np.where(G[i]==-1)[0]:
   for k in np.where((G[j]==-1)&(G[i]==0))[0]:
    c+=1;ext+=np.count_nonzero((G[i]==0)&(G[j]==0)&(G[k]==0))+N-len(r)
 return c,ext
for name,r in [('A5',rootsA(5)),('D4',rootsD(4))]:print(name,len(r),count(r))
a=count(rootsA(5))[1];d=count(rootsD(4))[1];print('D4^6',6*d,'A5^4D4',4*a+d,'difference',5*d-4*a)
