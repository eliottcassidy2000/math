from itertools import permutations, combinations
from functools import lru_cache, reduce
from collections import Counter
import random,json,hashlib

def need(x,m):
 if not x:raise RuntimeError(m)
def orient(n,code):
 a=[[False]*n for _ in range(n)];seen=set();k=0
 for i in range(n):
  for j in range(i+1,n):
   if (i,j) in seen:continue
   orb={(i,j),(n-1-j,n-1-i)};seen|=orb
   for u,v in orb:a[u][v]=bool(code>>k&1);a[v][u]=not a[u][v]
   k+=1
 return a,k
def paths(a):
 n=len(a)
 @lru_cache(None)
 def ending(mask,v):
  rem=mask^(1<<v)
  return 1 if not rem else sum(ending(rem,u) for u in range(n) if rem>>u&1 and a[u][v])
 return ending

def oddcycles(a):
 n=len(a);cs=[]
 def visit(p,used):
  if len(p)>=3 and len(p)%2 and a[p[-1]][p[0]]:cs.append(tuple(p))
  for v in range(p[0]+1,n):
   if not used>>v&1 and a[p[-1]][v]:visit(p+[v],used|(1<<v))
 for root in range(n):visit([root],1<<root)
 return cs

def poly_for_items(items,n):
 # Vertex-pivot recurrence, not recursive enumeration of combinations.
 @lru_cache(None)
 def rec(mask):
  if not mask:return (1,)
  bit=mask&-mask;out=list(rec(mask^bit))
  for support in items:
   if support&bit and support&mask==support:
    p=rec(mask^support)
    if len(out)<len(p)+1:out += [0]*(len(p)+1-len(out))
    for i,x in enumerate(p):out[i+1]+=x
  return tuple(out)
 return rec((1<<n)-1)
def canonical(c):
 j=c.index(min(c));return c[j:]+c[:j]
def inspect(a):
 n=len(a);full=(1<<n)-1;dp=paths(a);h=sum(dp(full,v) for v in range(n));cs=oddcycles(a);idx={c:i for i,c in enumerate(cs)};m=[sum(1<<v for v in c) for c in cs]
 tau=[idx[canonical(tuple(n-1-v for v in c[::-1]))] for c in cs]
 alpha=poly_for_items(m,n);need(h==sum((1<<k)*x for k,x in enumerate(alpha)),'OCF')
 fixed=[i for i,t in enumerate(tau) if t==i];pairs=[(i,t) for i,t in enumerate(tau) if i<t];good=[(i,j) for i,j in pairs if not m[i]&m[j]];I=len(pairs)-len(good)
 need((h-1-2*len(fixed)-4*I)%8==0,'mod8')
 if n%2==0:need(h%4==1 and not fixed,'even')
 # Independently enumerate all feasible packing sets only for the involution check.
 fixed_poly=Counter({0:1})
 for k in range(1,n//3+1):
  for pack in combinations(range(len(cs)),k):
   masks=[m[i] for i in pack]
   if sum(x.bit_count() for x in masks)!=reduce(int.__or__,masks,0).bit_count():continue
   if set(pack)=={tau[i] for i in pack}:fixed_poly[k]+=1
 q=[m[i]|m[j] for i,j in good];expected=Counter()
 for k,x in enumerate(poly_for_items(q,n)):expected[2*k]+=x
 for f in fixed:
  for k,x in enumerate(poly_for_items([x for x in q if not x&m[f]],n)):expected[2*k+1]+=x
 need(+fixed_poly==+expected,'fixed polynomial')
 need(all((x-fixed_poly[k])%2==0 and x>=fixed_poly[k] for k,x in enumerate(alpha)),'moving orbits')
 return h,len(cs)

def port_and_seam(a):
 n=len(a);full=(1<<n)-1;dp=paths(a);h=sum(dp(full,v) for v in range(n));rev=paths([list(t) for t in zip(*a)])
 for v in range(n):
  S=rev(full,v);R=dp(full,v)
  for mask in range(1,full):
   if mask>>v&1:R+=dp(mask,v)*sum(dp(full^mask,u) for u in range(n) if (full^mask)>>u&1)
  need((R-S)%2==0,'port parity')
  b=[[False]*(2*n) for _ in range(2*n)]
  for i in range(n):
   for j in range(n):b[i][j]=a[i][j];b[n+i][n+j]=a[j][i];b[n+i][j]=True
  b[n+v][v]=False;b[v][n+v]=True
  hd=paths(b);hh=sum(hd((1<<(2*n))-1,u) for u in range(2*n))
  need(hh==h*h-S*S+R*R,'literal mirror seam')
 return n
def symmetry_controls(a):
 n=len(a);ident=tuple(range(n));cs=oddcycles(a)
 def comp(p,q):return tuple(p[q[i]] for i in range(n))
 def power(p,k):
  r=ident
  for _ in range(k):r=comp(r,p)
  return r
 autos=[];antis=[]
 for p in permutations(range(n)):
  if all(a[p[i]][p[j]]==a[i][j] for i in range(n) for j in range(n)):autos.append(p)
  if all(a[p[i]][p[j]]==a[j][i] for i in range(n) for j in range(n)) and comp(p,p)==ident:antis.append(p)
 need(len(autos)%2==1 and antis,'automorphism order and existence')
 stats=set()
 for p in antis:
  maps={c:canonical(tuple(p[v] for v in c[::-1])) for c in cs}
  F=sum(c==d for c,d in maps.items());I=sum(c<d and bool(set(c)&set(d)) for c,d in maps.items());stats.add((F,I))
 need(len(stats)==1,'choice-independent cycle orbit data')
 sigma=antis[0]
 for tau in antis:
  g=comp(tau,sigma);m=1
  while power(g,m)!=ident:m+=1
  need(m%2==1,'odd product order')
  h=power(g,(m+1)//2)
  need(h in autos and comp(h,sigma)==comp(tau,h),'explicit conjugator')
 return len(antis)
rng=random.Random(291037);cases=0;hist={};anti_checks=0
for n in range(1,8):
 _,d=orient(n,0);codes=list(range(1<<d)) if n<=5 else rng.sample(range(1<<d),40)
 counts=Counter()
 for code in codes:
  a=orient(n,code)[0];h,c=inspect(a);counts[h%8]+=1;cases+=1
  if n<=5:anti_checks+=symmetry_controls(a)
 hist[n]={'cases':len(codes),'mod8':dict(counts)}
ports=0
for n in range(1,5):
 es=list(combinations(range(n),2))
 for code in range(1<<len(es)):
  a=[[False]*n for _ in range(n)]
  for k,(i,j) in enumerate(es):a[i][j]=bool(code>>k&1);a[j][i]=not a[i][j]
  ports+=port_and_seam(a)
print(json.dumps({'status':'PASS','self_converse_cases':cases,'antiinvolutions_conjugacy_checked':anti_checks,'orders':hist,'literal_port_and_mirror_seam_cases':ports},sort_keys=True,indent=2))
