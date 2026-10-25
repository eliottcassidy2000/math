#!/usr/bin/env python3
import itertools as I, json, hashlib, pathlib, math
from fractions import Fraction
SRC=pathlib.Path(__file__).resolve().parent.parent/'tournament-hierarchy'
def require(ok, msg):
 if not ok: raise RuntimeError(msg)
def mat(n,m):
 a=[[0]*n for _ in range(n)]
 for k,(u,v) in enumerate(I.combinations(range(n),2)): a[u][v]=(m>>k)&1;a[v][u]=1-a[u][v]
 return a
def mods(a):
 n=len(a); full=(1<<n)-1; out=[]
 for s in range(1,full+1):
  verts=[i for i in range(n) if s>>i&1]
  if all(all(a[x][y]==a[x][verts[0]] for y in verts) for x in range(n) if not s>>x&1):out.append(s)
 return out
def strongtree(a):
 mm=mods(a); st=[s for s in mm if all(not s&t or (s&t)==s or (s&t)==t for t in mm)]
 root=(1<<len(a))-1
 parent={s:min((t for t in st if s!=t and s&t==s),key=int.bit_count) for s in st if s!=root}
 ch={s:[] for s in st}
 for s,p in parent.items():ch[p].append(s)
 return st,parent,ch
class Index:
 def __init__(self,ch,root):
  self.ch=ch;self.parent={}; self.depth={};self.euler=[];self.first={};self.up={}
  def visit(v,p,d):
   self.parent[v]=p;self.depth[v]=d;self.up[v]=[p]
   for j in range(1,(len(ch)+1).bit_length()):self.up[v].append(self.up[self.up[v][j-1]][j-1] if p!=v else v)
   self.first[v]=len(self.euler);self.euler.append(v)
   for c in ch[v]:visit(c,v,d+1);self.euler.append(v)
  visit(root,root,0);self.rmq=[self.euler[:]]
  j=1
  while 1<<j<=len(self.euler):
   prev=self.rmq[-1]; span=1<<(j-1)
   self.rmq.append([min((prev[i],prev[i+span]),key=self.depth.get) for i in range(len(self.euler)-(1<<j)+1)]);j+=1
 def lca(self,u,v):
  l,r=sorted((self.first[u],self.first[v])); k=(r-l+1).bit_length()-1
  return min((self.rmq[k][l],self.rmq[k][r-(1<<k)+1]),key=self.depth.get)
 def port(self,v,m):
  delta=self.depth[v]-self.depth[m]-1
  require(delta>=0,'port ancestor')
  for j in range(delta.bit_length()):
   if delta>>j&1:v=self.up[v][j]
  return v

def h(a):return sum(all(a[p[k]][p[k+1]] for k in range(len(a)-1)) for p in I.permutations(range(len(a))))
def flux(a):return tuple(sum(1 if a[i][j] else -1 for i in range(k) for j in range(k,len(a))) for k in range(len(a)+1))
def connected(a):
 n=len(a)
 for s in range(n):
  seen={s}; old=-1
  while old!=len(seen):old=len(seen);seen|={v for u in list(seen) for v in range(n) if a[u][v]}
  if len(seen)!=n:return False
 return True
report={'status':'PASS','exhaustive':[],'source_pins':{}}
for manifest in ['SHA256SUMS','source-hashes.txt']:
 for line in (SRC/manifest).read_text().splitlines():
  digest,name=line.split();actual=hashlib.sha256((SRC/name).read_bytes()).hexdigest();require(actual==digest,'hash '+name); report['source_pins'][name]=actual
cert=json.loads((SRC/'certificate.json').read_text())
counts=[];strongcounts=[];primecounts=[]
for n in range(1,6):
 N=1<<(n*(n-1)//2);sc=pc=0
 for mask in range(N):
  a=mat(n,mask);st,pa,ch=strongtree(a);root=(1<<n)-1;ix=Index(ch,root)
  require(len(st)<=2*n-1,'node bound'); require(sum(len(c)-1 for c in ch.values() if c)==n-1,'arity identity')
  bits=pb=0;w=0
  for s,cc in ch.items():
   if not cc:continue
   require(len(cc)>=2 and sum(cc)==s,'child partition')
   reps=[(c&-c).bit_length()-1 for c in cc]; q=[[a[u][v] for v in reps] for u in reps];k=len(q)
   trans=not any(q[u][v] and q[v][z] and q[z][u] for u,v,z in I.permutations(range(k),3))
   prime=len(mods(q))==k+1
   require(trans or prime,'quotient classification');bits+=k*(k-1)//2
   if not trans:pb+=k*(k-1)//2;w=max(w,k)
  require(bits<=n*(n-1)//2 and 2*pb<=w*(n-1),'payload bound')
  for u,v in I.permutations(range(n),2):
   m=ix.lca(1<<u,1<<v);oracle=min((s for s in st if s>>u&1 and s>>v&1),key=int.bit_count)
   require(m==oracle,'LCA');cu=ix.port(1<<u,m);cv=ix.port(1<<v,m)
   ru=(cu&-cu).bit_length()-1;rv=(cv&-cv).bit_length()-1
   require(cu!=cv and cu in ch[m] and cv in ch[m] and a[u][v]==a[ru][rv],'arc reconstruction')
  f=flux(a);scores=list(map(sum,a));require(f==tuple(sum(2*scores[i]-(n-1) for i in range(k)) for k in range(n+1)),'cut identity')
  sc+=connected(a);pc+=n>=3 and len(mods(a))==n+1
 counts.append(N);strongcounts.append(sc);primecounts.append(pc);report['exhaustive'].append({'n':n,'tournaments':N,'strongly_connected':sc,'prime':pc})
require(sum(counts)==1099,'total')
# Labelled SEQ recurrence: n! coefficient convolution, not ordinary counts.
A=[1]+counts;S=[0]+strongcounts
for n in range(1,6):require(A[n]==sum(math.comb(n,k)*S[k]*A[n-k] for k in range(1,n+1)),'labelled SCC SEQ')
# Labelled prime substitution via ordered block-size compositions.
P=[0,0,0]+primecounts[2:]
def compositions(n,k):
 if k==0:
  if n==0:yield ()
 else:
  for i in range(1,n-k+2):
   for rest in compositions(n-i,k-1):yield (i,)+rest
for n in range(1,6):
 rhs=int(n==1)
 for k in range(3,n+1):
  for sizes in compositions(n,k):rhs+=Fraction(P[k]*math.factorial(n)*math.prod(A[i] for i in sizes), math.factorial(k)*math.prod(math.factorial(i) for i in sizes))
 require(rhs==S[n],'labelled prime species')
report['prime_star_H']=[]
for ex in cert['prime_star_H_examples']:
 a=ex['matrix'];require(a==mat(5,ex['mask']),'prime mask');require(len(mods(a))==6,'prime');require(h(a)==ex['H'],'prime H');report['prime_star_H'].append(h(a))
require(sorted(report['prime_star_H'])==[9,13,15],'prime examples')
report['fixed_Hamilton_order']=[]
for n in range(1,6):
 free=[(u,v) for u,v in I.combinations(range(n),2) if v!=u+1];seen={};collisions=0
 for mask in range(1<<len(free)):
  a=[[int(u<v) for v in range(n)] for u in range(n)]
  for k,(u,v) in enumerate(free):a[u][v]=(mask>>k)&1;a[v][u]=1-a[u][v]
  f=flux(a);hh=h(a)
  if f in seen and hh not in seen[f]:collisions+=1
  seen.setdefault(f,set()).add(hh)
 require((collisions==0)==(n<=4),'minimality')
 report['fixed_Hamilton_order'].append({'n':n,'universe':1<<len(free),'H_separating_collisions':collisions})
witness=cert['cut_flux_search'][-1]['collision'];aa=witness['first']['matrix'];bb=witness['second']['matrix']
require((h(aa),h(bb))==(11,13),'witness H');require(list(map(sum,aa))==list(map(sum,bb))==[3,2,1,2,2],'scores')
require(flux(aa)==flux(bb)==(0,2,2,0,0,0),'witness flux')
cycle=[0,4,1,3];edges={tuple(sorted((cycle[i],cycle[(i+1)%4]))) for i in range(4)}
require(all(aa[cycle[i]][cycle[(i+1)%4]] for i in range(4)),'directed cycle')
require({(u,v) for u,v in I.combinations(range(5),2) if aa[u][v]!=bb[u][v]}==edges,'cycle flip')
T=[1,1,1,2,4,12,56,456];ss=[0]*8
for n in range(1,8):ss[n]=T[n]-sum(ss[k]*T[n-k] for k in range(1,n))
require(ss==[0,1,0,1,1,6,35,353],'unlabelled inverse')
report['notes']=['No asserts: same validations under python -O.','Independent sparse-table RMQ and binary-lifting child navigation verified; producer uses linear scans.','Source hashes verified; no all-size enumeration or novel complexity claim.']
print(json.dumps(report,indent=2))
