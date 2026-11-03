import itertools
from collections import deque

def paths(n,mask):
 pairs=list(itertools.combinations(range(n),2)); edges={(i,j):(mask>>k)&1 for k,(i,j) in enumerate(pairs)}
 def beat(i,j): return edges[min(i,j),max(i,j)] if i<j else 1-edges[j,i]
 return sum(all(beat(a,b) for a,b in zip(p,p[1:])) for p in itertools.permutations(range(n)))
def calc(n,mask):
 pairs=list(itertools.combinations(range(n),2)); h=paths(n,mask); f=[]
 for s in range(1,2**n-1):
  ext=mask
  # extension pair bit indexes not contiguous: rebuild
  allp=list(itertools.combinations(range(n+1),2)); ext=0
  for k,(i,j) in enumerate(allp):
   if j<n: val=(mask>>pairs.index((i,j)))&1
   else: val=not ((s>>i)&1)
   ext|=int(val)<<k
  f.append((s,paths(n+1,ext)))
 return h, sorted(set(v for s,v in f))
def strong(n,m):
 pairs=list(itertools.combinations(range(n),2)); E={(i,j):((m>>k)&1) for k,(i,j) in enumerate(pairs)}
 def beat(i,j):return E[i,j] if i<j else 1-E[j,i]
 for v in range(n):
  seen={v}
  while True:
   more=seen|{j for i in seen for j in range(n) if j!=i and beat(i,j)}
   if more==seen: break
   seen=more
  if len(seen)<n:return False
 return True
for n in range(3,6):
 first=None; count=0
 for m in range(2**(n*(n-1)//2)):
  if not strong(n,m): continue
  count+=1
  h,vals=calc(n,m)
  holes=[v for v in range(min(vals),max(vals)+1,2) if v not in vals]
  if holes:
   first=(n,m,h,vals,holes); break
 if n<5 and first is not None: raise RuntimeError(first)
 if n==5 and first != (5,8,9,[15,17,19,23,25,27,29,33,37,41],[21,31,35,39]): raise RuntimeError(first)
 print('order',n,'strong checked',count,'first nonsolid',first)
n=6; pairs=list(itertools.combinations(range(n),2)); m=0
for k,(i,j) in enumerate(pairs):
 val=1 if i//3<j//3 else int((i%3,j%3)!=(0,2))
 m|=val<<k
h, image=calc(n,m)
if (h,image)!=(9,[9,15,31,41,45,61,75,99]): raise RuntimeError((h,image))
print('C3 -> C3', 'code',m,'H',h,'ear image',image)
