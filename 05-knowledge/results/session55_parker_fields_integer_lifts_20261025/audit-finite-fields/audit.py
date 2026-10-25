"""Independent tuple-polynomial field arithmetic; no producer imports."""
import itertools,json,math
from pathlib import Path
class K:
 def __init__(self,p,f):
  self.p=p; self.f=f; self.n=len(f)-1; self.E=list(itertools.product(range(p),repeat=self.n)); self.z=(0,)*self.n;self.one=(1,)+(0,)*(self.n-1)
  self.sq={self.mul(x,x) for x in self.E}; self.ch={x:0 if x==self.z else 1 if x in self.sq else -1 for x in self.E}
  assert len(self.sq)==(len(self.E)+1)//2
  assert all(self.pow(x,len(self.E)-1)==self.one for x in self.E if x!=self.z)
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
 def factors(self,u,v):return [self.add(self.one,self.add(self.scale(s,u),self.scale(t,v))) for s,t in [(1,0),(-1,0),(1,1),(1,-1),(-1,1),(-1,-1)]]
rows=[]
# Distinct independent irreducible presentations, verified via x^(q-1)=1.
for p,f in [(3,[0,1]),(5,[0,1]),(7,[0,1]),(11,[0,1]),(29,[0,1]),(3,[1,0,1]),(5,[2,0,1]),(3,[1,2,0,1]),(7,[1,0,1]),(3,[2,1,0,0,1]),(5,[1,1,0,1]),(3,[1,2,0,0,0,1])]:
 F=K(p,f);q=len(F.E);exceptions={F.scale(a,F.one) for a in [0,1,-1,2,-2]};admiss=[]
 for v in F.E:
  # Root criterion, including all excluded and small-characteristic collapses.
  roots=[F.scale(s,F.add(F.one,F.scale(t,v))) for s,t in [(1,0),(-1,0),(1,1),(-1,1),(1,-1),(-1,-1)]]
  assert (len(set(roots))==6)==(v not in exceptions),(p,f,v)
  if v!=F.z:
   bad={F.scale(a,v) for a in [0,1,-1,2,-2,pow(2,-1,p),-pow(2,-1,p)]}
   for u in F.E:
    offsets={F.add(F.scale(s,u),F.scale(t,v)) for s in [-1,0,1] for t in [-1,0,1]}
    assert (len(offsets)==9)==(u not in bad),(p,f,u,v)
  if v not in exceptions and all(F.ch[F.add(F.one,F.scale(s,v))]==1 for s in [-1,1]):admiss.append(v)
 assert q<=25 or admiss
 checks=[]
 for v in admiss[:1]+admiss[-1:]:
  vals=[F.factors(u,v) for u in F.E]; sums=[]
  for mask in range(1,64):
   d=mask.bit_count(); S=0
   for L in vals:
    a=F.one
    for i in range(6):
     if mask>>i&1:a=F.mul(a,L[i])
    S+=F.ch[a]
   assert S*S<=(d-1)**2*q
   if d==1:assert S==0
   if d==2:assert abs(S)==1
   sums.append(S)
  W=sum(math.prod(1+F.ch[x] for x in L) for L in vals)
  assert W==q+sum(sums)
  assert W>=q-15-114*math.sqrt(q)
  zero_weight=sum(math.prod(1+F.ch[x] for x in L) for L in vals if F.z in L)
  assert zero_weight<=192
  checks.append({'v':v,'weighted_sum_numerator':W,'root_weight_numerator':zero_weight})
 rows.append({'p':p,'q':q,'modulus':f,'admissible_v_count':len(admiss),'checks':checks})
 print('PASS',q,flush=True)
Path(__file__).with_name('audit-results.json').write_text(json.dumps({'status':'PASS','tests':rows},indent=2)+'\n')
