from itertools import permutations,combinations,product
from collections import defaultdict

def hp(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for i,x in enumerate(d[s]):
   if x:
    t=a[i]&~s
    while t:
     b=t&-t;t-=b;d[s|b][b.bit_length()-1]+=x
 return sum(d[-1])
def profile(a,v):
 rest=[i for i in range(len(a)) if i!=v];m=len(rest);out=[0]*(m+2)
 for p in permutations(rest):
  for mask in range(1,1<<(m+1)):
   ok=True
   for j in range(m+1):
    if mask>>j&1:
     if j and not (a[p[j-1]]>>v&1):ok=False;break
     if j<m and not (a[v]>>p[j]&1):ok=False;break
    elif 0<j<m and not (a[p[j-1]]>>p[j]&1):ok=False;break
   if ok:out[mask.bit_count()]+=1
 return tuple(out[1:])
def rows(n,mask):
 a=[0]*n
 for b,(i,j) in enumerate(combinations(range(n),2)):
  if mask>>b&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a
def strong(a):
 n=len(a)
 for i in range(n):
  s=1<<i
  while True:
   t=s
   for j in range(n):
    if s>>j&1:t|=a[j]
   if t==s:break
   s=t
  if s!=(1<<n)-1:return False
 return True
if __name__=='__main__':
 g=defaultdict(dict)
 for mask in range(1024):
  a=rows(5,mask)
  if strong(a):
   p=profile(a,0);g[p[1:]].setdefault(p[0],(a,p))
 for tail,d in g.items():
  for x in d:
   if x+2 in d:print(d[x],d[x+2])
def inflate(a,v,r):
 others=[x for x in range(len(a)) if x!=v];n=r+len(others);b=[0]*n
 for i in range(r):
  for j in range(i+1,r):b[i]|=1<<j
 for i,x in enumerate(others,r):
  for j,y in enumerate(others,r):
   if a[x]>>y&1:b[i]|=1<<j
  for j in range(r):
   if a[v]>>x&1:b[j]|=1<<i
   else:b[i]|=1<<j
 return b
if __name__=='__main__':
 a=[12,17,2,6,13]
 for r in range(1,8):
  b=inflate(a,0,r);target=r+1;vals=[hp(b)]
  for j in range(r):
   b[j]^=1<<target;b[target]^=1<<j;vals.append(hp(b))
  print(r,vals)
