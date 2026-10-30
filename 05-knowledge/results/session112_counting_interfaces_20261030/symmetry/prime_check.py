"""Prime fixed-bank convolution versus literal signed ideal-cut intersections."""
from math import gcd
from itertools import product
from collections import Counter
import json
CHECKS=0
def need(ok,why):
 global CHECKS
 CHECKS+=1
 if not ok: raise RuntimeError(why)

def scalar_counts(p,m,a):
 h=(p-1)//2
 A=[1]+[1<<(h*sum(gcd(a,p**((l-1)//2)) for l in range(1,j+1))) for j in range(1,m+1)]
 P=[0]*(m+1)
 for j in range(1,m+1):P[j]=A[j]-sum(P[r]*A[j-r] for r in range(1,j))
 return A,P

def prime_fixed(p,k,a,b):
 m=2*k;n=p**(k-1);c=n//p;h=(p-1)//2
 a%=n; shifted=(a+c*b)%n
 C=(p-1)*sum(gcd(a,p**j) for j in range(k-1))+h*(gcd(a,n)+gcd(shifted,n))
 A,P=scalar_counts(p,m,a)
 Ash,_=scalar_counts(p,m,shifted)
 return (1<<C)-P[1]*Ash[m-1]-sum(P[r]*A[m-r] for r in range(2,m))

def grouped_sum(p,k):
 n=p**(k-1);c=n//p
 total=0;rows=[]
 for v in range(k-2):
  mult=(p-1)*p**(k-1-v)
  val=prime_fixed(p,k,p**v,0)
  rows.append({'type':'valuation','v':v,'multiplicity':mult,'prime_fixed':val})
  total+=mult*val
 for az,bz in product((True,False),repeat=2):
  a=0 if az else c
  shifted=0 if bz else c
  b=((shifted-a)//c)%p
  mult=(1 if az else p-1)*(1 if bz else p-1)
  val=prime_fixed(p,k,a,b)
  rows.append({'type':'top_pair','a_zero':az,'shifted_zero':bz,'multiplicity':mult,'prime_fixed':val})
  total+=mult*val
 need(sum(r['multiplicity'] for r in rows)==p**k,'grouped mass')
 need(total%p**k==0,'prime orbit integrality')
 return total//p**k,rows

class DSU:
 def __init__(self,n):self.pa=list(range(n));self.par=[0]*n;self.n=n;self.ok=True
 def find(self,x):
  if self.pa[x]!=x:
   parent=self.pa[x];self.pa[x]=self.find(parent);self.par[x]^=self.par[parent]
  return self.pa[x]
 def add(self,a,b,e):
  ra,rb=self.find(a),self.find(b);v=self.par[a]^self.par[b]^e
  if ra==rb:
   if v:self.ok=False
   return False
  self.pa[ra]=rb;self.par[ra]=v;self.n-=1;return True

class Ring:
 def __init__(self,p,k,d=1):
  self.p,self.k,self.q,self.D=p,k,p**k,p*d
  q=self.q
  self.points=list(product(range(q),repeat=2))
  inv=pow(1-self.D,-1,q)
  g=((1+self.D)*inv%q,2*inv%q)
  owner={};orbits=[]
  for z in self.points[1:]:
   if z in owner:continue
   orb=[];w=z
   while w not in owner:owner[w]=len(orbits);orb.append(w);w=self.mul(g,w)
   need(w==z,'clock orbit closes');orbits.append(orb)
  self.lookup={};self.reps=[]
  for orb in orbits:
   z=orb[0]
   if z in self.lookup:continue
   i=len(self.reps);self.reps.append(z)
   for w in orb:
    self.lookup[w]=(i,0);self.lookup[(-w[0]%q,-w[1]%q)]=(i,1)
  self.bits=len(self.reps)
  need(self.bits==p**k-1,'literal bank dimension')
  self.cuts=[]
  for r in range(1,2*k):
   ux,vy=p**((r+1)//2),p**(r//2)
   generators=[(ux%q,0),(0,vy%q)]
   eqs=[];ds=DSU(self.bits)
   for z in self.points[1:]:
    if z[0]%ux==0 and z[1]%vy==0:continue
    i,e=self.lookup[z]
    for a,b in generators:
     w=((z[0]+a)%q,(z[1]+b)%q)
     j,f=self.lookup[w]
     if ds.add(i,j,e^f):eqs.append((i,j,e^f))
   need(ds.ok,'consistent ideal cut')
   self.cuts.append(eqs)
 def mul(self,z,w):
  q=self.q
  return ((z[0]*w[0]+self.D*z[1]*w[1])%q,(z[0]*w[1]+z[1]*w[0])%q)
 def action(self,a,b):
  q=self.q;p=self.p;k=self.k
  scale=pow(1+p,a,q);twist=pow(1+p**(k-1),b,q)
  return [self.lookup[(scale*x%q,scale*twist*y%q)] for x,y in self.reps]
 def count(self,action,mask):
  ds=DSU(self.bits)
  for i,(j,e) in enumerate(action):ds.add(i,j,e)
  for r,eqs in enumerate(self.cuts):
   if mask>>r&1:
    for i,j,e in eqs:ds.add(i,j,e)
  return 1<<ds.n if ds.ok else 0

reports=[]
for p,k,d in ((5,2,1),(5,2,2),(7,2,1),(5,3,1)):
 ring=Ring(p,k,d);n=p**(k-1);hist=Counter();total=0
 for a,b in product(range(n),range(p)):
  action=ring.action(a,b)
  vals=[ring.count(action,mask) for mask in range(1<<(2*k-1))]
  prime=sum((-1 if mask.bit_count()%2 else 1)*val for mask,val in enumerate(vals))
  need(prime==prime_fixed(p,k,a,b),'fixed prime convolution versus literal cut inclusion-exclusion')
  hist[prime]+=1;total+=prime
 predicted,rows=grouped_sum(p,k)
 need(total//p**k==predicted and total%p**k==0,'literal and grouped orbit totals')
 reports.append({'p':p,'k':k,'D':p*d,'bank_bits':ring.bits,'prime_orbits':str(predicted),
                 'fixed_prime_histogram':{str(a):b for a,b in sorted(hist.items())},'grouped_rows':[{**r,'prime_fixed':str(r['prime_fixed'])} for r in rows]})
need(reports[0]['prime_orbits']=='665856','625 prime O2 orbit control')
A,P=scalar_counts(5,4,0)
need(65536-sum(P[r]*A[4-r] for r in range(1,4))==-65280,'omitted inner twist is a negative hostile')
# Independent rational-tail recurrence identities for fixed finite valuation.
for p in (5,7):
 for v in range(1,4):
  A,P=scalar_counts(p,2*v+8,p**v)
  c=1<<(((p-1)//2)*p**v)
  Q=[1]+[A[i]-c*A[i-1] for i in range(1,2*v+1)]
  need(Q[-1]!=0,'scalar recurrence exact degree')
  for m in range(2*v+1,len(P)):
   need(sum(Q[i]*P[m-i] for i in range(len(Q)))==0,'scalar rational-tail recurrence')
_,Pv1=scalar_counts(5,6,5)
need(Pv1[5]==17045913600 and Pv1[6]==17454749184000,'old depth-five hostile extends to new depth-six row')
large=[]
for p in (5,7,13):
 for k in range(2, (8 if p==5 else 6 if p==7 else 4)+1):
  val,rows=grouped_sum(p,k)
  need(all(int(r['prime_fixed'])>=0 for r in rows),'nonnegative fixed prime strata')
  large.append({'p':p,'k':k,'prime_orbit_bits':val.bit_length(),'prime_orbits':str(val) if k<=3 else None})
print(json.dumps({'status':'PASS','checks':CHECKS,'literal_families':reports,'formula_families':large},indent=2))
