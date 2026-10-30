#!/usr/bin/env python3
"""Exact CRT/Gauss spectral census; standard library, explicit failures."""
import itertools,json,math
from fractions import Fraction
from pathlib import Path

def require(x,msg):
 if not x: raise RuntimeError(msg)
def factor(n):
 out=[];p=2
 while p*p<=n:
  a=0
  while n%p==0:n//=p;a+=1
  if a:out.append((p,a))
  p+=1
 if n>1:out.append((n,1))
 return out
def divisors(n):return [d for d in range(1,n+1) if n%d==0]
def chi(x,p):return 0 if x%p==0 else (1 if pow(x%p,(p-1)//2,p)==1 else -1)
def responses(r):
 reps=[(1,)+s for s in itertools.product((-1,1),repeat=r-1)]
 for vals in itertools.product((-1,1),repeat=len(reps)):
  table={}
  for s,v in zip(reps,vals):table[s]=v;table[tuple(-x for x in s)]=-v
  coeff={}
  for mask in range(1<<r):
   c=Fraction(sum(v*math.prod(s[i] for i in range(r) if mask>>i&1) for s,v in table.items()),1<<r)
   if c:coeff[mask]=c
  yield table,coeff

def contribution(N,m,coeff,k,primes):
 fs=factor(m);out=[Fraction(0)]*(1<<len(primes))
 for mask,c in coeff.items():
  v=c;globalmask=0
  for j,(p,a) in enumerate(fs):
   pa=p**a;h=k*pow(m//pa,-1,pa)%pa
   if mask>>j&1:
    if h%(p**(a-1)) or h==0:v=0;break
    v*=p**(a-1)*chi(h//p**(a-1),p);globalmask|=1<<primes.index(p)
   else:
    v*= (p**(a-1)*(p-1) if h==0 else -p**(a-1) if h%p**(a-1)==0 else 0)
  out[globalmask]+=v
 return out

def multiply(a,b,ps):
 out=[0]*len(a)
 for i,x in enumerate(a):
  if not x:continue
  for j,y in enumerate(b):
   if y:out[i^j]+=x*y*math.prod(-p for h,p in enumerate(ps) if (i&j)>>h&1)
 return out

def symbol(N,ds,choices):
 out=[0]*N
 for m,(table,_) in zip(ds,choices):
  ps=[p for p,a in factor(m)]
  for t in range(1,m):
   if math.gcd(t,m)==1:out[N//m*t]=table[tuple(chi(t,p) for p in ps)]
 return out

def direct4(s):
 n=len(s)
 return n*sum(sum(s[x]*s[(x+h)%n] for x in range(n))**2 for h in range(n))

def moment(N,ds,choices):
 ps=[p for p,a in factor(N)];total=0
 for d in ds:
  eig=[Fraction(0)]*(1<<len(ps))
  for m,(_,cs) in zip(ds,choices):
   eig=[a+b for a,b in zip(eig,contribution(N,m,cs,N//d,ps))]
  sq=multiply(eig,eig,ps)
  total+=sum(math.gcd(u,d)==1 for u in range(d))*multiply(sq,sq,ps)[0]
 return total

def census231():
 N=231;ps=[3,7,11];ds=divisors(N)[1:];opts=[list(responses(len(factor(m)))) for m in ds]
 weights=[sum(math.gcd(u,d)==1 for u in range(d)) for d in ds]
 pre=[]
 for m,os in zip(ds,opts):
  row=[]
  for table,cs in os:
   vals=[[2*x for x in contribution(N,m,cs,N//d,ps)] for d in ds]
   require(all(x.denominator==1 for v in vals for x in v),'scale2 integrality')
   row.append([[int(x) for x in v] for v in vals])
  pre.append(row)
 # Independent full census: construct actual positive arc masks, rotate bits,
 # and count intersections. No Gauss sums, Walsh coefficients, or radicals.
 bitopts=[]
 for m,os in zip(ds,opts):
  row=[]
  for table,cs in os:
   pp=[p for p,a in factor(m)]
   row.append(sum(1<<(N//m*t) for t in range(1,m) if math.gcd(t,m)==1 and table[tuple(chi(t,p) for p in pp)]==1))
  bitopts.append(row)
 hist={'character':{},'nonlinear':{}};witness={};checked=[]
 for indices in itertools.product(*(range(len(o)) for o in opts)):
  total=0
  for di,w in enumerate(weights):
   v=[sum(pre[j][ix][di][b] for j,ix in enumerate(indices)) for b in range(8)]
   a,b,c,d=v[1],v[2],v[4],v[7]
   C=-(3*a*a+7*b*b+11*c*c+231*d*d)
   total+=w*(C*C+21*(2*a*b-22*c*d)**2+33*(2*a*c-14*b*d)**2+77*(2*b*c-6*a*d)**2)
  require(total%16==0,'moment integrality');total//=16
  bits=0
  for j,ix in enumerate(indices):bits|=bitopts[j][ix]
  require(bits.bit_count()==115,'positive degree')
  autoc=(N-1)**2;allbits=(1<<N)-1
  for h in range(1,N):
   rotated=((bits<<h)|(bits>>(N-h)))&allbits
   corr=2-N+4*(bits&rotated).bit_count()
   autoc+=corr*corr
  require(N*autoc==total,'independent full census autocorrelation')
  typ='character' if len(opts[-1][indices[-1]][1])==1 else 'nonlinear'
  h=hist[typ];h[total]=h.get(total,0)+1
  if typ not in witness or total<witness[typ][0]:witness[typ]=(total,indices)
 for typ,(value,inds) in witness.items():
  s=symbol(N,ds,[opts[j][ix] for j,ix in enumerate(inds)])
  require(direct4(s)==value,'231 independent min check');checked.append({'kind':typ,'indices':inds,'moment4':value})
 return {'N':N,'family_size':8192,'independent_bit_autocorrelation_checks':8192,'summaries':{t:{'count':sum(h.values()),'min':min(h),'max':max(h),'min_multiplicity':h[min(h)],'distinct_moments':len(h),'histogram':{str(k):v for k,v in sorted(h.items())}} for t,h in hist.items()},'direct_min_checks':checked}

def main():
 controls=[]
 for N in (3,7,9,21,27,33,49):
  ds=divisors(N)[1:];opts=[list(responses(len(factor(m)))) for m in ds];count=0;vals={}
  for choices in itertools.product(*opts):
   s=symbol(N,ds,choices);m=moment(N,ds,choices)
   require(m==direct4(s),f'direct moment N={N}')
   require(all(s[-i]==-s[i] for i in range(1,N)),'oddness')
   count+=1;vals[int(m)]=vals.get(int(m),0)+1
  controls.append({'N':N,'symbols_checked':count,'moments':vals})
 receipt={'status':'PASS','controls':controls,'census':census231()}
 print(json.dumps(receipt,sort_keys=True,indent=2))
if __name__=='__main__':main()
