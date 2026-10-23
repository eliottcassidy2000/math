import itertools, collections, sympy as s
from fractions import Fraction as F

def classes(n,perms=None):
 E=list(itertools.combinations(range(n),2));idx={e:i for i,e in enumerate(E)}
 perms=list(itertools.permutations(range(n))) if perms is None else perms
 maps=[]
 for p in perms:
  mp=[idx[tuple(sorted((p[a],p[b])))] for a,b in E]
  flip=sum(1<<mp[k] for k,(a,b) in enumerate(E) if p[a]>p[b]);maps.append((mp,flip))
 G={};T={}
 for x in range(1<<len(E)):
  ys=[sum(1<<mp[k] for k in range(len(E)) if x>>k&1) for mp,_ in maps]
  G[x]=min(ys);T[x]=min(y^maps[j][1] for j,y in enumerate(ys))
 return E,G,T

def chain(A,G,T,paired=False):
 full=max(G)
 g=lambda x:min(G[x],G[full^x]) if paired else G[x]
 t=lambda x:min(T[x],T[full^x]) if paired else T[x]
 cells=collections.Counter((g(x),t(x)) for x in A);gs=sorted({i for i,j in cells});ts=sorted({j for i,j in cells})
 C=s.Matrix([[cells[i,j] for j in ts] for i in gs]);r=[sum(C.row(i)) for i in range(C.rows)];c=[sum(C.col(j)) for j in range(C.cols)]
 K=s.diag(*[1/s.Integer(v) for v in r])*C*s.diag(*[1/s.Integer(v) for v in c])*C.T
 return gs,ts,C,K
if __name__=='__main__':
 for n in range(2,6):
  E,G,T=classes(n);path=sum(1<<E.index((i,i+1)) for i in range(n-1));A=[x for x in G if x&path==0]
  for paired in [False,True]:
   gs,ts,C,K=chain(A,G,T,paired)
   print(n,paired,len(A),K.shape, 'poly',s.factor(K.charpoly().as_expr()))
   if n<=4:print('G',gs,'T',ts,'C',C,'K',K)
