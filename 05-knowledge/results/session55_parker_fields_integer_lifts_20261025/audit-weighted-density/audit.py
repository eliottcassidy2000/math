#!/usr/bin/env python3
"""Independent closed controls; no producer imports; exceptions survive -O."""
import collections, hashlib, itertools, json, math
from pathlib import Path
from fractions import Fraction as Q
D=Path(__file__).resolve().parent; R=D.parent

def ck(x,s):
 if not x: raise RuntimeError(s)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def grid(c,u,v,n):return tuple((c+s*u+t*v)%n for s,t in [(1,0),(-1,-1),(0,1),(-1,1),(0,0),(1,-1),(0,-1),(1,1),(-1,0)])
def signature(g):
 names={};return ''.join(str(names.setdefault(v,len(names))) for v in g)
def histogram(n):
 sq={x*x%n for x in range(n)};out=collections.Counter()
 for c,u,v in itertools.product(range(n),repeat=3):
  g=grid(c,u,v,n)
  if all(x in sq for x in g):out[signature(g)]+=1
 return out
def source_path(path):
 mapping={'session55-crt-collision-wiring':'crt-wiring','session55-finite-field-geometry':'finite-fields','session55-prime-power-counts':'prime-power-counts'}
 if path.startswith('math/05-knowledge/results/'):
  return R.parent/path.removeprefix('math/05-knowledge/results/')
 first,rest=path.split('/',1)
 return R/mapping.get(first,first)/rest
pins=[]
for path,h in json.loads((R/'weighted-crt/SOURCE_PINS.json').read_text()).items():
 ck(digest(source_path(path))==h,'weighted source pin '+path);pins.append(path)
for item in json.loads((R/'density-bridge/source-pins.json').read_text()):
 ck(digest(source_path(item['source']))==item['sha256'],'density original source pin');ck(digest(R/'density-bridge'/item['copy'])==item['sha256'],'density copy pin');pins.append(item['source'])
cert=json.loads((R/'weighted-crt/certificate.json').read_text());hist={}
for n in [1,2,3,4,5,7,8,9,11,16,25,27,29,35,55]:
 h=histogram(n);hist[n]=h
 row=cert['census'][n-1];ck(sum(h.values())==row['total_arrays'] and h['012345678']==row['distinct_labeled'],'independent full plane '+str(n))
for n in [35,55]:
 h=hist[n];target={r['partition_rgs']:r['count'] for r in cert['examples'][str(n)]['full_histogram']};ck(h==target,'full histogram')
 ck(h['012345678']=={35:24,55:40}[n],'explicit distinct count')
ck(hist[3]=={'000000000':2},'mod3 no root weights')
# All 3^9 raw arrays, no parameterization, directly test eight line equations.
lines=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6));raw=collections.Counter()
for g in itertools.product(range(3),repeat=9):
 if set(g)<={0,1} and len({sum(g[i] for i in L)%3 for L in lines})==1:raw[signature(g)]+=1
ck(raw==hist[3],'raw nine-cell mod3')
# Independent partition generation by inserting the next position into blocks.
def parts(blocks,i):
 if i==9:yield blocks;return
 for j in range(len(blocks)):yield from parts(blocks[:j]+[blocks[j]+[i]]+blocks[j+1:],i+1)
 yield from parts(blocks+[[i]],i+1)
ps=list(parts([],0));ck(len(ps)==21147,'Bell9')
for n in [35,55]:
 a,b=cert['examples'][str(n)]['factors'];total=0
 for blocks in ps:
  def z(h):return sum(w for g,w in h.items() if all(len({g[i] for i in B})==1 for B in blocks))
  zn=z(hist[n]);ck(zn==z(hist[a])*z(hist[b]),'zeta multiplicative')
  mu=(-1)**(9-len(blocks))*math.prod(math.factorial(len(B)-1) for B in blocks);total+=mu*zn
 ck(total==hist[n]['012345678'],'Mobius direction and sign')
# Standalone tuple polynomial fields, explicit group-law validation.
class Field:
 def __init__(self,p,f):
  self.p=p;self.f=f;self.n=len(f)-1;self.E=list(itertools.product(range(p),repeat=self.n));self.z=(0,)*self.n;self.one=(1,)+(0,)*(self.n-1);self.q=len(self.E)
  self.sq={self.mul(x,x) for x in self.E};self.ch={x:0 if x==self.z else 1 if x in self.sq else -1 for x in self.E}
  ck(len(self.sq)==(self.q+1)//2,'field square count');ck(all(self.pow(x,self.q-1)==self.one for x in self.E if x!=self.z),'field multiplicative orders')
 def add(self,a,b):return tuple((x+y)%self.p for x,y in zip(a,b))
 def scale(self,k,a):return tuple(k*x%self.p for x in a)
 def mul(self,a,b):
  c=[0]*(2*self.n-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):c[i+j]+=x*y
  for i in range(len(c)-1,self.n-1,-1):
   for j in range(self.n):c[i-self.n+j]-=c[i]*self.f[j]
  return tuple(x%self.p for x in c[:self.n])
 def pow(self,a,n):
  b=self.one
  while n:
   if n%2:b=self.mul(b,a)
   a=self.mul(a,a);n//=2
  return b
fieldrows=[];density=json.loads((R/'density-bridge/certificate.json').read_text());known={x['q']:x for x in density['field_counts']}
for p,f in [(3,[0,1]),(5,[0,1]),(29,[0,1]),(59,[0,1]),(3,[1,0,1]),(3,[1,2,0,1]),(7,[1,0,1]),(3,[2,1,0,0,1]),(5,[1,1,0,1])]:
 F=Field(p,f);q=F.q;one=F.one;E={F.scale(i,one) for i in [0,1,-1,2,-2]};A={v for v in F.E if all(F.ch[F.add(one,F.scale(s,v))]==1 for s in [-1,1])};G=A-E
 ck(Q(len(A))==Q(q-F.ch[F.scale(-1,one)],4)-Q(1+F.ch[F.scale(2,one)],2),'exact conic identity')
 counts={};good=0;terms=0
 for c in [F.z,one]:
  C=0
  for v in F.E:
   N=0
   for u in F.E:
    g=[F.add(c,F.add(F.scale(s,u),F.scale(t,v))) for s,t in itertools.product([-1,0,1],repeat=2)]
    valid=len(set(g))==9 and all(x in F.sq for x in g)
    C+=valid
    if c==one and v in G and valid and F.z not in g:N+=1
   if c==one and v in G:
    ck(max(0,abs(64*N-q)-655)**2<=114**2*q,'two-sided slice');good+=N;terms+=1
  counts[c]=C
 ck(counts[one]==known[q]['center_one'] and counts[F.z]==known[q]['center_zero'],'independent full normalized field planes')
 ck(0<=counts[one]-good<=13*q,'bad-v and zero-cell remainder')
 total=(q-1)//2*counts[one]+counts[F.z];ck(total==known[q]['total'],'field total')
 fieldrows.append({'q':q,'C':counts[one],'Z':counts[F.z],'total':total,'good_slices':terms})
ck(Q(1,512)/(1-Q(1,64))==Q(1,504),'dyadic density');ck(Q(1,27)/(1-Q(1,729))==Q(27,728),'ternary density');ck(Q(27,728)*512==Q(1728,91),'ternary field ratio')
out={'status':'PASS','source_pins':pins,'weighted_controls':{str(n):dict(h) for n,h in hist.items()},'partition_count':len(ps),'independent_fields':fieldrows,'normal_optimized_safe':True}
s=json.dumps(out,sort_keys=True,indent=2)+'\n'
if (D/'independent-results.json').exists():ck((D/'independent-results.json').read_text()==s,'repeat results')
else:(D/'independent-results.json').write_text(s)
print('PASS: independent weighted, raw mod3, partition inversion, tuple-field, source pins and density controls')
