import itertools,json,hashlib
from collections import defaultdict

def canonical(a,aps):
 n=len(a)
 if sorted(map(sum,a))==list(range(n)):return ('transitive',n)
 codes=[]
 for p in aps:
  rest=[i for i in range(n) if i!=p]
  rest=sorted(rest,key=lambda i:-sum(a[i][j] for j in rest))
  codes.append(''.join(str(a[p][j]) for j in rest))
 return min(codes)

def inspect(a):
 n=len(a); aps=set(range(n)); cycles=0
 for i,j,k in itertools.combinations(range(n),3):
  if a[i][j]==a[j][k]==a[k][i]:
   cycles+=1;aps.intersection_update((i,j,k))
   if not aps:return None
 return aps

out=[]
for n in range(8,9):
 pairs=[(i,j) for i in range(n) for j in range(i+2,n)]
 counts=defaultdict(int); blues=defaultdict(int); conv={}; internal=0; loops=0
 for mask in range(1<<len(pairs)):
  a=[[0]*n for _ in range(n)]
  for i in range(n-1):a[i][i+1]=1
  for b,(i,j) in enumerate(pairs):a[i][j]=(mask>>b)&1;a[j][i]=1-a[i][j]
  aps=inspect(a)
  if aps is None:continue
  c=canonical(a,aps);counts[c]+=1
  ac=[[a[j][i] for j in range(n)] for i in range(n)]
  conv[c]=canonical(ac,aps)
  blue=all(a[i][j]==a[n-1-j][n-1-i] for i in range(n) for j in range(i+1,n))
  blues[c]+=blue
  kap=[r[:] for r in a]
  for i,j in pairs:kap[i][j]^=1;kap[j][i]^=1
  kas=inspect(kap)
  if kas is not None:
   internal+=1;kc=canonical(kap,kas)
   loops+=kc==c or kc==conv[c]
 r=n//2
 A=2**(n-1)-(n-1)*(n-2)//2-1
 C=r if n%2==0 else 2**r+r-1
 V=4**(n-2)-3*2**(n-2)-(n*n-7*n+4)//2
 B=2**r+r-2 if n%2==0 else 2**(2*r-1)+r-1
 result=dict(n=n,A=len(counts),C=sum(c==conv[c] for c in counts),V=sum(counts.values()),B=sum(blues.values()),only_blue_classes=[str(c) for c in counts if counts[c]==blues[c]],internal_tilings=internal,loop_tilings=loops)
 assert [result[k] for k in ['A','C','V','B']]==[A,C,V,B],result
 out.append(result);print(json.dumps(result),flush=True)
json.dump(out,open('/workspace/shared/session26-blue-review/enumeration_n8.json','w'),indent=2)
