#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations,permutations,product
from functools import lru_cache
from fractions import Fraction as Q
import hashlib,json
D=Path(__file__).resolve().parent.parent

def check(t,msg):
 if not t: raise RuntimeError(msg)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def act(m,s):return tuple(dot(m[3*i:3*i+3],s) for i in range(3))
def mul(a,b):return tuple(sum(a[3*i+k]*b[3*k+j] for k in range(3)) for i in range(3) for j in range(3))
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
G=[tuple(map(int,l.split('|')[0].split())) for l in (D/'orbit/generators_base3_15.txt').read_text().splitlines() if not l.startswith('#')]
F=[tuple(map(int,l.split())) for l in (D/'last_letter/matrices_base3_order4.txt').read_text().splitlines()]
T=G[0]; v=(3,1,2); seed=(1,0,1)
# Reconstruct each transfer using three physically independent states, not its path-sum formula.
@lru_cache(None)
def state(adj,port):
 n=len(adj); full=(1<<n)-1; cycles=[]
 for k in range(3,n+1,2):
  for sub in combinations(range(n),k):
   for tail in permutations(sub[1:]):
    p=(sub[0],)+tail
    if all(adj[p[j]]>>p[(j+1)%k]&1 for j in range(k)):cycles.append(sum(1<<x for x in p))
 @lru_cache(None)
 def H(mask):
  if not mask:return 1
  bit=mask&-mask
  return H(mask^bit)+3*sum(H(mask^c) for c in cycles if c&bit and c&mask==c)
 vals=[0,0]
 def paths(p,mask,length):
  vals[length%2]+=H(full^mask)
  for q in range(n):
   if not mask>>q&1 and adj[p]>>q&1: paths(q,mask|1<<q,length+1)
 paths(port,1<<port,1)
 return (H(full),*vals)
def tournaments(n):
 pairs=list(combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  a=[0]*n
  for z,(i,j) in enumerate(pairs):
   if bits>>z&1:a[i]|=1<<j
   else:a[j]|=1<<i
  yield tuple(a)
inputs=[(0,),(2,0),(6,4,0)]
check([state(a,0) for a in inputs]==[(1,0,1),(1,1,1),(1,2,2)],'independent input basis')
def physical_matrix(b,w,up):
 ys=[]
 for a in inputs:
  n=len(a); k=len(b)
  adj=[x|(((1<<k)-1)<<n) for x in a]+[x<<n for x in b]
  adj[0]^=1<<(n+w);adj[n+w]|=1
  ys.append(state(tuple(adj),n+up))
 e=tuple(y-x for x,y in zip(ys[0],ys[1]));o=tuple(z-x-2*y for x,y,z in zip(ys[0],e,ys[2]));h=tuple(x-y for x,y in zip(ys[0],o))
 return tuple(t for row in zip(h,e,o) for t in row)
rebuilt=set();count=0
for k in range(1,5):
 for b in tournaments(k):
  for w,up in product(range(k),repeat=2):rebuilt.add(physical_matrix(b,w,up));count+=1
check(count==1105 and rebuilt==set(F),'physical alphabet reconstruction')
check(set(G)<=set(F),'15 subset of 57')
print('PASS: direct physical gluing of 1105 blocks on three basis states independently reconstructs exactly 57 matrices.',flush=True)
# Geometric certificates with independently calculated supporting rows.
full=[(1,0,1),(97,32,33),(25,17,17),(1,8,8)]
nonempty=[(193,33,129),(292,65,162),*full[1:]]
words=[(),(4,6),(4,2),(4,)]
qwords=[(4,6,1),(4,6,5),*words[1:]]
from math import gcd
facets=[]
for rays,ws in [(full,words),(nonempty,qwords)]:
 fs=[]
 for j,r in enumerate(rays):
  s=seed
  for i in ws[j]:s=act(G[i-1],s)
  check(s==r,'actual extreme word')
  f=cross(r,rays[(j+1)%len(rays)])
  if dot(f,rays[(j+2)%len(rays)])<0:f=tuple(-x for x in f)
  d=gcd(*f);f=tuple(x//d for x in f);fs.append(f)
  check(all(dot(f,z)>0 for z in rays if z!=r and z!=rays[(j+1)%len(rays)]),'strict convex polygon')
 check(all(dot(f,act(m,r))>=0 and act(m,r)[0]>0 for f in fs for m in G for r in rays),'invariant hull')
 check(all(dot(f,act(m,seed))>=0 for f in fs for m in G),'initial images')
 facets.append(fs)
 print('PASS: exact orbit cone',rays,'facets',fs,flush=True)
check(all(min(dot(f,(10,4,5)) for f in fs)>=0 for fs in facets),'hostile cone point')
check(min(Q(r[1],r[0]) for r in nonempty)==Q(33,193),'sharp E bound')
check(min(Q(act(m,r)[0],r[0]) for m in G for r in nonempty)==Q(292,193),'sharp growth')
check(act(T,nonempty[0])[0]==292,'growth equality')
for rs in (full,nonempty):
 check(max(Q(r[1],r[0]) for r in rs)==8 and max(Q(r[2],r[0]) for r in rs)==8 and max(Q(r[0],r[2]) for r in rs)==Q(97,33) and max(Q(r[1]+r[2],r[0]) for r in rs)==16,'sharp ratios')
print('PASS: all sharp cone ratios and growth 292/193, attained by T at (193,33,129).',flush=True)
# Independent kernel classification: nonzero cross product of the two annihilator rows spans preimage line.
powers=[];p=(1,0,0,0,1,0,0,0,1)
for k in range(4):p=mul(T,p);powers.append(p)
classes={'tower':[],'nonnegative_kernel':[],'excluded':[]}
for i,m in enumerate(F):
 check(min(m)>=0 and m[0]>0 and min(3*m[6+j]-m[j] for j in range(3))>=0,'invariance C')
 a=tuple(m[j]-3*m[3+j] for j in range(3));b=tuple(m[6+j]-2*m[3+j] for j in range(3));z=cross(a,b)
 # Rank-one singular annihilator systems need the elementary positive-row obstruction.
 if z==(0,0,0):
  check(any(r[0]>0 and min(r)>=0 for r in (a,b,tuple(-x for x in a),tuple(-x for x in b))),'rank-one obstruction')
  classes['excluded'].append(i);continue
 if z[0]<0:z=tuple(-x for x in z)
 if m in powers:
  check(cross(z,v)==(0,0,0),'tower preimage');classes['tower'].append(i)
 elif z[0]>0 and min(z)>=0:
  check(cross(z,(45,1,5))==(0,0,0) and z[0]>3*z[2],'exception');classes['nonnegative_kernel'].append(i)
 else:classes['excluded'].append(i)
check([len(classes[k]) for k in classes]==[4,1,52],'classification')
check(all(F[i] not in G for i in classes['nonnegative_kernel']),'15 stronger obstruction')
check(act(T,v)==tuple(2*x for x in v),'tower doubling')
I=(1,0,0,0,1,0,0,0,1);ell=(1,2,1)
check(powers[2]==tuple(I[3*i+j]+v[i]*ell[j] for i in range(3) for j in range(3)),'tower cubed')
print('PASS: independent cross-product kernel classification',classes,flush=True)
# Unsaturated exact closure, retaining the full state.
seen={seed};queue=[seed]
for s in queue:
 for m in G:
  t=act(m,s)
  if t[0]<=3000 and t not in seen:seen.add(t);queue.append(t)
gaps=sorted(set(range(1,3001,3))-{s[0] for s in seen})
check(len(seen)==53732 and gaps==[10,31,61,106,217,274,337,346,607,865],'finite closure')
print('PASS: unsaturated closure',len(seen),'states, gaps',gaps)
# Check the distributed certificate artifact, independently of its generator.
cert=json.loads((D/'last_letter/CHECK.json').read_text())
check(cert['alphabet_sha256']==hashlib.sha256((D/'last_letter/matrices_base3_order4.txt').read_bytes()).hexdigest(),'frozen certificate hash')
for rec in cert['records']:
 m=F[rec['index']]
 if rec['kind']=='separator':
  a=tuple(m[j]-3*m[3+j] for j in range(3));b=tuple(m[6+j]-2*m[3+j] for j in range(3))
  f=tuple(rec['p']*x+rec['q']*y for x,y in zip(a,b))
  check(tuple(rec['positive_form'])==f and f[0]>0 and min(f)>=0,'stored separator')
# Cramer determinants check every printed rational inverse witness independently.
def det(m):return dot(m[:3],cross(m[3:6],m[6:9]))
expected={1:('3/2','1/2','1'),7:('0','0','1/2'),8:('1','-1/3','1/3'),9:('-33/10','5/2','-7/10'),10:('57/104','7/104','-5/13'),11:('1','-1/3','-1/6'),12:('0','1','-1'),13:('0','0','1/3'),14:('51/131','10/131','-19/131'),15:('3/71','14/71','-5/71')}
for idx,m in enumerate(G,1):
 d=det(m)
 if d:
  x=[]
  for j in range(3):
   a=list(m)
   for i in range(3):a[3*i+j]=v[i]
   x.append(Q(det(a),d))
  check(tuple(x)==tuple(map(Q,expected[idx])),'printed inverse table')
 else:
  r=(0,1,-1) if idx in (2,3,4) else (1,-2,-1)
  check(all(sum(r[i]*m[3*i+j] for i in range(3))==0 for j in range(3)) and dot(r,v)!=0,'singular separator')
P=(1,1,0,0,0,1,0,1,1);Pinv=(1,1,-1,0,-1,1,0,1,0)
check(mul(P,Pinv)==I,'unimodular inverse')
for m in G:
 check(min(mul(mul(Pinv,m),P))>=0,'elementary cone')
 check(m[0]>=1 and m[1]>=3,'growth first row')
check(act(F[54],(45,1,5))==(393,131,262),'exceptional hostile')
source=D.parents[0]/'kps_s18_ternary_chains_digraph_weights_20261109/out/generators_base3_15.txt'
check(source.read_bytes()==(D/'orbit/generators_base3_15.txt').read_bytes(),'inherited bytes unchanged')
print('PASS: stored CHECK.json certificates, Cramer inverse table, unimodular cone and inherited generator bytes.')
print('ALL INDEPENDENT CHECKS PASS')
