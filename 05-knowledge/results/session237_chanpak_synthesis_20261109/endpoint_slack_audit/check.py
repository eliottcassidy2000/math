"""Literal permutation check, template-isomorphism equality recognition."""
import itertools as it,json

def require(x,msg):
 if not x:raise RuntimeError(msg)
def graph(n,m):
 a=[[False]*n for _ in range(n)]
 for k,(i,j) in enumerate(it.combinations(range(n),2)):
  a[i][j]=bool(m&(1<<k));a[j][i]=not a[i][j]
 return a
def endpoints(a):
 E=[0]*len(a)
 for p in it.permutations(range(len(a))):
  if all(a[x][y] for x,y in zip(p,p[1:])):E[p[-1]]+=1
 return E
def templates(n):
 out=set()
 for prefix in range(n-2):
  a=[[i<j for j in range(n)] for i in range(n)]
  a[prefix][n-1]=False;a[n-1][prefix]=True
  for p in it.permutations(range(n)):
   code=sum(1<<k for k,(i,j) in enumerate(it.combinations(range(n),2)) if a[p[i]][p[j]])
   out.add((code,p.index(n-2)))
 return out
import sys
require('--fail' not in sys.argv,'deliberate checker failure')
counts=[];tested=0;zero=0
for n in range(1,6):
 expected=templates(n);eq=0
 for code in range(1<<(n*(n-1)//2)):
  a=graph(n,code);E=endpoints(a);h=sum(E)
  for u,e in enumerate(E):
   if not any(a[u]):continue
   if not e:
    if h==1:zero+=1
    continue
   tested+=1;same=h==2*e+1
   require(same==((code,u) in expected),f'classify {n},{code},{u}')
   require(h>=2*e+1,'bound')
   if same:eq+=1;require(e&(e-1)==0,'power2')
 counts.append([n,eq])
checks=0
for l in range(7):
 for mask in range(1<<l):
  n=l+2;a=[[i<j for j in range(n)] for i in range(n)]
  for i in range(l):a[-1][i]=bool(mask&(1<<i));a[i][-1]=not a[-1][i]
  E=endpoints(a);h=sum(E);e=E[l]
  S=[i for i in range(l) if mask&(1<<i)]
  m=sum(1<<(l-1-i) for i in S)
  q=sum(2 if j==i+1 else 1<<(j-i-1) for i,j in it.combinations(S,2))
  require((e,h,h-2*e)==(m-q,1+2*m-q,1+q),'apex')
  checks+=1
print(json.dumps({'status':'PASS','positive_nonsink_rooted_cases':tested,'equality_by_order':counts,'zero_endpoint_nonsink_equality_cases':zero,'literal_apex_masks_l_le_6':checks},indent=2))
