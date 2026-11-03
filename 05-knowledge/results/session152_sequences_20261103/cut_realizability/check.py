def require(x):
 if not x: raise AssertionError("cut restriction check failed")
from itertools import combinations, permutations
from collections import Counter

def adjacency(n,code):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(code>>k)&1; a[j][i]=1-a[i][j]
 return a

def hp(a):
 n=len(a); d=[[0]*n for _ in range(1<<n)]
 for i in range(n): d[1<<i][i]=1
 for s in range(1,1<<n):
  for i in range(n):
   if s>>i&1:
    d[s][i]+=sum(d[s^(1<<i)][j]*a[j][i] for j in range(n) if j!=i and s>>j&1)
 return sum(d[-1])

def test(n,code):
 a=adjacency(n,code); perms=list(permutations(range(n)))
 good=[]; bad=[]; K=[[0]*n for _ in range(n)]
 for p in perms:
  breaks=[k for k in range(n-1) if not a[p[k]][p[k+1]]]
  if not breaks: good.append(p); marks=range(n-1)
  elif len(breaks)==1: bad.append((p,breaks[0])); marks=breaks
  else: continue
  for k in marks:
   i,j=p[k:k+2]; K[i][j]+=1;K[j][i]+=1
 H=len(good); require(H==hp(a)) 
 require(sum(K[i][j] for i,j in combinations(range(n),2))==(n-1)*H+len(bad)) 
 for m in range(1, (1<<n)-1):
  cross=lambda i,j: bool((m>>i&1)!=(m>>j&1))
  C=sum(sum(cross(p[k],p[k+1]) for k in range(n-1)) for p in good)
  B=sum(cross(p[k],p[k+1]) for p,k in bad)
  W=sum(K[i][j] for i,j in combinations(range(n),2) if cross(i,j))
  require(W==C+B and W>=H and W%2==0) 
  cyc=[p for p in good if a[p[-1]][p[0]]]
  if cyc: require(W>=2*((H+n-1)//2)) # ceil((H+n-2)/2)*2
 F=[]
 for m in range(1<<n):
  child=[row+[int(not(m>>i&1))] for i,row in enumerate(a)]+[[m>>i&1 for i in range(n)]+[0]]
  F.append(hp(child))
 for m in range(1,(1<<n)-1):
   W=sum(K[i][j] for i,j in combinations(range(n),2) if bool((m>>i&1)!=(m>>j&1)))
   require(F[m]+F[((1<<n)-1)^m]==2*H+W) 
 return H,len(bad),F

for n in range(2,5):
 for code in range(1<<(n*(n-1)//2)): test(n,code)
 print('exhaustive n=',n,'count=',1<<(n*(n-1)//2))
for code in (8,759):
 h,b,f=test(5,code); print('n=5 code=',code,'H=',h,'B=',b,'image=',sorted(set(f[1:-1])))
 print('target masks:',{v:[m for m in range(1,31) if f[m]==v] for v in (3*h,3*h+2,3*h+4)})
print('PASS')
