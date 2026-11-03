def require(x):
 if not x: raise AssertionError("verification failed")
from itertools import combinations
from collections import Counter

def adj(n,c):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(c>>k)&1; a[j][i]=1-a[i][j]
 return a

def paths(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for m in range(1,1<<n):
  for i in range(n):
   if m>>i&1:
    d[m][i]+=sum(d[m^(1<<i)][j]*a[j][i] for j in range(n) if m>>j&1 and j!=i)
 return sum(d[-1])
def ear(a,s):
 n=len(a);return [row+[int(not(s>>i&1))] for i,row in enumerate(a)]+[[s>>i&1 for i in range(n)]+[0]]
for c in [8,759]:
 a=adj(5,c); F=[paths(ear(a,s)) for s in range(32)]; H=paths(a); A=[F[1<<i]-H for i in range(5)]; K={(i,j):F[1<<i]+F[1<<j]-F[(1<<i)|(1<<j)]-H for i,j in combinations(range(5),2)}
 print(c,H,Counter(F[1:-1]),A,K, 'constant',F[0],F[-1]);print('by size',*[dict(sorted(Counter(F[s] for s in range(32) if s.bit_count()==k).items())) for k in range(6)])
a=adj(5,8);F=[paths(ear(a,s)) for s in range(32)];print('pointwise complement',all(F[s]==F[31^s] for s in range(32)));print('edges',[(i,j) if a[i][j] else (j,i) for i,j in combinations(range(5),2)]);print('masks',[(s,F[s]) for s in range(1,31)])
def strong(a):
 n=len(a)
 for i in range(n):
  R={i}
  while 1:
   U=R|{j for x in R for j in range(n) if a[x][j]}
   if U==R:break
   R=U
  if len(R)<n:return False
 return True
families=Counter(); hs=Counter();
for c in range(1024):
 a=adj(5,c)
 if strong(a):
  H=paths(a);hs[H]+=1
  im=tuple(sorted({paths(ear(a,s)) for s in range(1,31)}));families[(H,im)]+=1
print('strong parent H census',hs)
print('families',len(families));print(*families.items(),sep='\n')
a=adj(5,8); f=[paths(ear(a,s)) for s in range(32)]
require(paths(a)==9)
require(Counter(f[1:31])=={15:2,17:4,19:2,23:6,25:2,27:2,29:2,33:4,37:4,41:2})
require(hs=={9:240,11:120,13:120,15:64} and len(families)==6)
require(set.union(*(set(im) for h,im in families))=={15,17,19,23,25,27,29,31,33,37,41,43,45})
require(all(f[s]==9+sum([8,6,14,24,14][i] for i in range(5) if s>>i&1)-sum(k for (i,j),k in {(0,1):6,(0,2):4,(0,3):4,(0,4):8,(1,2):10,(1,3):10,(1,4):4,(2,3):10,(2,4):4,(3,4):6}.items() if s>>i&1 and s>>j&1) for s in range(32)))
print('ALL ASSERTIONS PASS')
