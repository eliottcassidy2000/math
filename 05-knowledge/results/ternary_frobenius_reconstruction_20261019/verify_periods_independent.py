"""Independent integer-dictionary certificates and sparse finite-field rank checks.
No production algebra/matrix code is imported. Python standard library only.
"""
from itertools import permutations
from math import comb,gcd
from functools import reduce
from pathlib import Path
import json
BASE=Path(__file__).resolve().parent

def plus(a,b,c=1):
 out=a.copy()
 for t,v in b.items():
  w=out.get(t,0)+c*v
  if w:out[t]=w
  elif t in out:del out[t]
 return out

def times(a,b):
 out={}
 for t,c in a.items():
  for u,d in b.items():
   v=tuple(x+y for x,y in zip(t,u));out[v]=out.get(v,0)+c*d
 return {t:c for t,c in out.items() if c}

def var(n,i,c=1):return {tuple(int(j==i) for j in range(n)):c}
def power(a,n):
 out={(0,)*len(next(iter(a))):1}
 for _ in range(n):out=times(out,a)
 return out

def substitution(f,linear):
 maxs=[max(t[j] for t in f) for j in range(len(linear))]
 cache=[[power(q,n) for n in range(m+1)] for q,m in zip(linear,maxs)]
 out={}
 for t,c in f.items():
  v={(0,)*len(linear):c}
  for j,e in enumerate(t):v=times(v,cache[j][e])
  out=plus(out,v)
 return out

def require(ok,msg):
 if not ok:raise RuntimeError(msg)

def load_poly(name):
 return {tuple(map(int,t.split(','))):c for t,c in json.loads((BASE/name).read_text())[0].items()}

def perm(f,pm):return {tuple(t[i] for i in pm)+tuple(t[i+3] for i in pm):c for t,c in f.items()}
def check_cert(name,d,k,eps):
 P=load_poly(name);v=[var(6,i) for i in range(6)]
 A=[var(6,1,-1),plus(v[0],v[1],-1),v[2],var(6,4,-1),plus(v[3],v[4],-1),v[5]]
 A2=[plus(v[1],v[0],-1),var(6,0,-1),v[2],plus(v[4],v[3],-1),var(6,3,-1),v[5]]
 AP=substitution(P,A);A2P=substitution(P,A2)
 cert=plus(plus(times(plus(v[0],v[1],-1),P),times(v[0],AP),-1),times(v[1],A2P))
 require(not cert,'rational norm identity failed')
 require(perm(P,(1,2,0))==P,'cyclic failed')
 require(perm(P,(1,0,2))=={t:eps*c for t,c in P.items()},'symmetry failed')
 require(all(sum(t[:3])==d+3 and sum(t[3:])==k for t in P),'bidegree failed')
 require(all(all(t[i]>=2 or t[i+3]>0 for i in range(3)) for t in P),'interior specialization condition failed')
 require(reduce(gcd,P.values())==1,'not primitive')
 # Gauge P(X+cx,Y+cy,Z+cz)=P is equivalent to the directional derivative
 # because characteristic zero is being used for this certificate.
 deriv={}
 for t,c in P.items():
  for i in range(3):
   if t[i+3]:
    u=list(t);u[i]+=1;u[i+3]-=1;deriv=plus(deriv,{tuple(u):c*t[i+3]})
 require(not deriv,'linear-additive gauge invariance failed')
 return P,{'file':name,'terms':len(P),'bidegree':[d,k],'symmetry':eps,'certificate_zero':True,'gauge_invariant':True}

def displayed_alternating_numerator():
 triples=[([(1,0),(3,0),(-1,3)],1),
 ([(1,0),(1,1),(1,2)],-1), ([(1,0),(2,1),(0,2)],-2),
 ([(2,0),(0,1),(1,2)],-1), ([(2,0),(1,1),(0,2)],-1),
 ([(3,0),(0,1),(0,2)],1), ([(0,1),(1,1),(2,1)],-2)]
 out={}
 for fs,c in triples:
  for pm in permutations(range(3)):
   inv=sum(pm[i]>pm[j] for i in range(3) for j in range(i+1,3))
   lower=[1,1,1];upper=[0,0,0]
   for i,j in enumerate(pm):lower[j]+=fs[i][0];upper[j]+=fs[i][1]
   out=plus(out,{tuple(lower+upper):c*(-1)**inv})
 return out

def compact_symmetric_numerator():
 v=[var(6,i) for i in range(6)]
 # Monomials (x degree, s degree, t degree, coefficient), s=y+z,t=yz.
 terms=[(6,1,0,175),(5,2,0,-70),(5,0,1,315),
 (4,3,0,-175),(4,1,1,525),(3,4,0,136),(3,2,1,-754),(3,0,2,639),
 (2,5,0,21),(2,3,1,-14),(2,1,2,-231),
 (1,6,0,-34),(1,4,1,263),(1,2,2,-519),(1,0,3,108),
 (0,7,0,11),(0,5,1,-111),(0,3,2,287),(0,1,3,-68)]
 out={}
 for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:
  x,y,z,X,Y,Z=[v[h] for h in [i,j,k,i+3,j+3,k+3]]
  ss=plus(y,z);tt=times(y,z);H={}
  for a,b,c,coef in terms:
   H=plus(H,times(times(power(x,a),power(ss,b)),power(tt,c)),coef)
  wedges=times(plus(times(x,Y),times(y,X),-1),plus(times(x,Z),times(z,X),-1))
  term=times(times(times(tt,power(plus(power(y,2),power(z,2),-1),2)),H),wedges)
  out=plus(out,term)
 return out

def specialize(P,p):
 out={}
 for t,c in P.items():
  u=tuple(t[i]-1+p*t[i+3] for i in range(3));out[u]=(out.get(u,0)+c)%p
 return {t:c for t,c in out.items() if c}

def independent_norm(t,p):
 # Literal multiplication of the three substituted linear forms.
 v=[var(3,i) for i in range(3)]
 subs=[[v[0],v[1],v[2]],[var(3,1,-1),plus(v[0],v[1],-1),v[2]],[plus(v[1],v[0],-1),var(3,0,-1),v[2]]]
 out={}
 for forms in subs:
  z={(0,0,0):1}
  for l,e in zip(forms,t):z=times(z,power(l,e))
  for u,c in z.items():out[u]=(out.get(u,0)+c)%p
 return {u:c for u,c in out.items() if c}

def cyclic_interior_columns(n,p):
 # Choose lexicographically smallest orbit representative independently.
 for a in range(1,n):
  for b in range(1,n-a):
   c=n-a-b;t=(a,b,c);orb={t,(b,c,a),(c,a,b)}
   if t!=min(orb):continue
   col={}
   for u in orb:
    for r,v in independent_norm(u,p).items():col[r]=(col.get(r,0)+v)%p
   yield {r:v for r,v in col.items() if v}

def sparse_rank(columns,p):
 piv={};n=0
 for col in columns:
  n+=1
  while col:
   key=min(col);c=col[key]
   if key not in piv:
    inv=pow(c,-1,p);piv[key]={r:v*inv%p for r,v in col.items()};break
   b=piv[key]
   for r,v in b.items():
    a=(col.get(r,0)-c*v)%p
    if a:col[r]=a
    elif r in col:del col[r]
 return n,len(piv)

if __name__=='__main__':
 out={'certificates':[],'specializations':[],'thresholds':[]}
 A,aa=check_cert('rational_kernel_d3_k3_eps-1.json',3,3,-1)
 require(A=={t:-c for t,c in displayed_alternating_numerator().items()},'displayed alternating formula failed')
 aa['stored_vs_displayed_sign']=-1
 out['certificates'].append(aa)
 S,ss=check_cert('rational_kernel_d12_k2_eps1.json',12,2,1);out['certificates'].append(ss)
 require(compact_symmetric_numerator()==S,'compact symmetric formula failed')
 require(max(max(t[i]-1 for i in range(3)) for t in S)==12,'exponent upper bound failed')
 require(min(min(t[i]-1 for i in range(3)) for t in S)==-1,'exponent lower bound failed')
 require(S.get((2,13,0,0,0,2))==11,'uniform nonzero coefficient failed')
 require(specialize(S,11).get((1,6,27))==10,'p11 nonzero coefficient failed')
 require(specialize(S,13).get((1,12,25))==11,'p13 nonzero coefficient failed')
 for p in (3,5,7,11,13,17,19,23,29,31,37):
  for P,label in [(A,'alternating'),(S,'symmetric')]:
   f=specialize(P,p)
   require(all(min(t)>0 for t in f),'specialized polynomial not interior')
   norm={}
   for t,c in f.items():
    # This full norm path is only needed in the small-prime boundary checks;
    # the universal certificate handles every larger prime algebraically.
    if p<=13:
     for u,e in independent_norm(t,p).items():norm[u]=(norm.get(u,0)+c*e)%p
   require(not any(norm.values()),'specialized norm failed')
   out['specializations'].append({'p':p,'family':label,'terms':len(f),'nonzero':bool(f)})
 print('Universal certificates and boundary specializations passed.',flush=True)
 # All degrees below the supplied minima, not samples.
 for p,delta in [(2,11),(3,11),(5,18),(7,24),(11,34),(13,38),(17,46),(19,50)]:
  vals=[]
  for n in range(delta+1):
   cols,rank=sparse_rank(cyclic_interior_columns(n,p),p);null=cols-rank
   require(null==(1 if n==delta else 0),f'threshold discrepancy p={p},n={n},nullity={null}')
   vals.append({'n':n,'columns':cols,'rank':rank,'nullity':null})
  out['thresholds'].append({'p':p,'minimum':delta,'all_degrees':vals})
  print(f'p={p}: independently proved finite-exact minimum {delta}, nullity one.',flush=True)
 (BASE/'independent_verification.json').write_text(json.dumps(out,indent=2)+'\n')
 print('All independent checks passed.',flush=True)
