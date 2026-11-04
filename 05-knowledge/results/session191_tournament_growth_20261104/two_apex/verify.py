#!/usr/bin/env python3
"""Exact two-apex transfer and its minimal rational continuation dimension.
No third-party packages; all checks survive python -O.
"""
import hashlib,itertools,json,random
from fractions import Fraction as Q
from pathlib import Path

def require(ok,msg):
 if not ok:raise ValueError(msg)

# Symbol c records a->x in bit 0 and b->x in bit 1, with a->b.
TRANS=[[[] for _ in range(36)] for _ in range(4)]
ACCEPT=[]
for o,x,y,z in itertools.product(range(2),range(3),range(3),range(2)):
 s=o*18+x*6+y*2+z
 if x!=2 and (y==1 or (y==0 and o==0)):ACCEPT.append(s)
 for c in range(4):
  p=(c>>o)&1;q=(c>>(1-o))&1
  TRANS[c][s].append(o*18+(2 if p else 1)*6+y*2+z)
  if y or p:TRANS[c][s].append(o*18+x*6+(2 if q else 1)*2+z)
  if z or q:TRANS[c][s].append(o*18+x*6+y*2+1)
INIT=(1,)+(0,)*17+(1,)+(0,)*17

def step(v,c):
 out=[0]*36
 for x,h in enumerate(v):
  for y in TRANS[c][x]:out[y]+=h
 return tuple(out)
def state(w):
 v=INIT
 for c in w:v=step(v,c)
 return v
def value(w):return sum(state(w)[i] for i in ACCEPT)
def graph(w):
 n=len(w)+2;rows=[0]*n;rows[0]=2
 for j,c in enumerate(w,2):
  for a in (0,1):
   if c>>a&1:rows[a]|=1<<j
   else:rows[j]|=1<<a
  for i in range(2,j):rows[i]|=1<<j
 return rows

def dp(rows):
 n=len(rows);d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if not(s>>v&1):continue
   old=s^(1<<v)
   if old:d[s][v]=sum(d[old][u] for u in range(n) if rows[u]>>v&1)
 return sum(d[-1])

def strong_rows(rows):
 n=len(rows);allmask=(1<<n)-1
 for source in range(n):
  seen=1<<source
  while True:
   new=seen
   for u in range(n):
    if seen>>u&1:new|=rows[u]
   if new==seen:break
   seen=new
  if seen!=allmask:return False
 return True

def strong_word(w):
 return bool(w) and w[0]!=0 and w[-1]!=3 and (2 in w or any(a==3 and b==0 for i,a in enumerate(w) for b in w[i+1:]))

def literal_slots(w):
 # Independent direct enumeration of a physical Hamilton word, not status states.
 rows=graph(w);total=0
 for order in ((0,1),(1,0)):
  for alloc in itertools.product(range(3),repeat=len(w)):
   word=[]
   for gap in range(3):
    word.extend(j+2 for j,a in enumerate(alloc) if a==gap)
    if gap<2:word.append(order[gap])
   total+=all(rows[u]>>v&1 for u,v in zip(word,word[1:]))
 return total

class Span:
 def __init__(self):self.rows={}
 def add(self,v):
  v=list(map(Q,v))
  for p,b in sorted(self.rows.items()):
   if v[p]:
    a=v[p];v=[x-a*y for x,y in zip(v,b)]
  try:p=next(i for i,a in enumerate(v) if a)
  except StopIteration:return False
  a=v[p];self.rows[p]=[x/a for x in v];return True

def determinant(matrix):
 a=[list(map(Q,r)) for r in matrix];d=Q(1);n=len(a)
 for j in range(n):
  p=next((i for i in range(j,n) if a[i][j]),None)
  if p is None:return 0
  if p!=j:a[p],a[j]=a[j],a[p];d=-d
  z=a[j][j];d*=z
  for i in range(j+1,n):
   c=a[i][j]/z
   for k in range(j+1,n):a[i][k]-=c*a[j][k]
  # No normalization: the diagonal pivots form the determinant.
 return int(d)

def main():
 receipt={'alphabet':{'0':'x->a,b','1':'a->x->b','2':'b->x->a','3':'a,b->x'}}
 # Every two-apex presentation with at most five chain vertices: 1+4+...+4^5.
 checks=0
 for m in range(6):
  for w in itertools.product(range(4),repeat=m):
   h=value(w);require(h==dp(graph(w))==literal_slots(w),'independent path count')
   require(h%2==1,'parity');require(strong_word(w)==strong_rows(graph(w)),'strongness characterization');checks+=1
 receipt['complete_independent_presentations']=checks
 rng=random.Random(191)
 samples=[]
 for m in range(6,11):
  for _ in range(8):
   w=tuple(rng.randrange(4) for _ in range(m));h=value(w)
   require(h==dp(graph(w)),'larger independent sample')
   samples.append({'word':''.join(map(str,w)),'H':h})
 receipt['larger_independent_samples']=samples
 # Reachable span is invariant because every child of every accepted basis vector
 # is tested, not because a fixed word depth happened to stop growing.
 R=[INIT];W=[()];span=Span();span.add(INIT);i=0
 while i<len(R):
  for c in range(4):
   v=step(R[i],c)
   if span.add(v):R.append(v);W.append(W[i]+(c,))
  i+=1
 require(len(R)==28,'reachable dimension')
 # All observable rows are restricted to the invariant reachable span.
 V=[()];rows=[[sum(v[j] for j in ACCEPT) for v in R]];sp=Span();sp.add(rows[0]);i=0
 while i<len(V):
  for c in range(4):
   suffix=(c,)+V[i];row=[value(w+suffix) for w in W]
   if sp.add(row):V.append(suffix);rows.append(row)
  i+=1
 require(len(V)==17,'minimal dimension')
 # Independent columns choose a square Hankel minor.
 cp=Span();ix=[]
 for j in range(len(W)):
  if cp.add([row[j] for row in rows]):ix.append(j)
 prefixes=[W[j] for j in ix]
 H=[[value(p+s) for p in prefixes] for s in V]
 for s,row in zip(V,H):
  for p,h in zip(prefixes,row):require(dp(graph(p+s))==h,'independent Hankel entry')
 det=determinant(H);require(det==150994944,'nonzero Hankel determinant')
 receipt.update(reachable_dimension=28,minimal_rational_dimension=17,
  reachable_basis_words=[''.join(map(str,w)) for w in W],
  hankel_prefixes=[''.join(map(str,w)) for w in prefixes],
  hankel_suffixes=[''.join(map(str,w)) for w in V],hankel=H,hankel_determinant=det)
 # Minimal all-order one-apex hostile, and constant-run spectral identity.
 require(value((1,2,0))==13,'two-apex 13')
 for c in range(4):
  for v in R:
   v1=step(v,c);v2=step(v1,c);v3=step(v2,c);v4=step(v3,c)
   require(all(d-6*e+11*f-6*g==0 for d,e,f,g in zip(v4,v3,v2,v1)), 'clock annihilator')
 receipt['constant_run_annihilator']='M^4-6M^3+11M^2-6M=0 on the reachable span'
 # These three initials plus the proved recurrence certify the all-k family.
 family=[value((1,2)+(0,)*k) for k in (1,2,3)]
 require(family==[13,33,85],'explicit two-apex family initials')
 receipt['explicit_family']={'word':'12 0^k','k_min':1,'H':'-1+4*2^k+2*3^k','initials':family}

 raw=(json.dumps(receipt,sort_keys=True,indent=2)+'\n').encode()
 Path(__file__).with_name('receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','exhaustive_presentations':checks,'random_larger_presentations':len(samples),
 'hankel_entries_independently_checked':289,'minimal_dimension':17,'determinant':det,
 'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
