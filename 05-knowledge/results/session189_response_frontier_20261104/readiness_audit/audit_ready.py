import itertools,json

def tournament(n,c):
 a=[0]*n
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if c>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a

def induced(a,vs):return [sum((1<<j) for j,v in enumerate(vs) if a[u]>>v&1) for u in vs]

def strong(a):
 n=len(a)
 for i in range(n):
  seen=1<<i
  while True:
   nxt=seen
   for j in range(n):
    if seen>>j&1:nxt|=a[j]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<n)-1:return False
 return True

def hpends(a,r):return sum(1 for p in itertools.permutations(range(len(a))) if p[-1]==r and all(a[p[i]]>>p[i+1]&1 for i in range(len(p)-1)))
def cycles(a):return sum(1 for p in itertools.permutations(range(1,len(a))) if all(a[x]>>y&1 for x,y in zip((0,)+p,p+(0,))))
def phi(a,v):
 b=a[:]
 for i,j in itertools.combinations([x for x in range(len(a)) if x!=v],2):b[i]^=1<<j;b[j]^=1<<i
 return b
res={}
for n in [4,5]:
 marked=cases=0
 for c in range(1<<(n*(n-1)//2)):
  a=tournament(n,c)
  if not strong(a):continue
  cy=cycles(a);ex=[v for v in range(n) if a[v].bit_count()==1];ad=[]
  for v in ex:
   vs=[u for u in range(n) if u!=v];core=induced(a,vs);r=vs.index((a[v]&-a[v]).bit_length()-1);digit=hpends(core,r)
   if digit and not strong(core):raise RuntimeError(('strong core',n,c,v,digit))
   if strong(core):ad.append(v)
   if digit:
    z=phi(a,v)
    if cycles(z)!=digit or not strong(z):raise RuntimeError(('involution',n,c,v))
    marked+=1
  if ex and cy in (1,2):
   if len(ad)!=1:raise RuntimeError(('canonical',n,c,cy,ad))
   cases+=1
 res[n]={'marked_positive':marked,'canonical_cases':cases}
print(json.dumps({'status':'PASS','by_order':res},sort_keys=True))
