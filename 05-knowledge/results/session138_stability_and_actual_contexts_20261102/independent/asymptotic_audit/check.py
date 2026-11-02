"""Exact rational independent checks for triangle + inherited-variance bounds."""
from fractions import Fraction as F
from itertools import combinations,product
import random,json,hashlib,functools

def require(p,message):
 if not p:raise RuntimeError(message)
def matrix(n,p):
 A=[[F(0) for _ in range(n)] for _ in range(n)]
 for (i,j),x in zip(combinations(range(n),2),p):A[i][j]=x;A[j][i]=1-x
 return A
def det(A):
 B=[list(map(F,r)) for r in A];out=F(1)
 for k in range(len(B)):
  r=next((r for r in range(k,len(B)) if B[r][k]),None)
  if r is None:return F(0)
  if r!=k:B[k],B[r]=B[r],B[k];out=-out
  x=B[k][k];out*=x
  for j in range(k,len(B)):B[k][j]/=x
  for i in range(k+1,len(B)):
   x=B[i][k]
   for j in range(k,len(B)):B[i][j]-=x*B[k][j]
 return out
def uv(A):
 n=len(A);return det([[A[i][j]+(i==j) for j in range(n)] for i in range(n)]),(-1)**n*det(A)
def all_h(A):
 n=len(A);d=[[F(0)]*n for _ in range(2**n)];h=[F(0)]*2**n;h[0]=1
 for i in range(n):d[1<<i][i]=1
 for mask in range(1,2**n):
  h[mask]=sum(d[mask])
  for i in range(n):
   if d[mask][i]:
    for j in range(n):
     if not mask>>j&1:d[mask|(1<<j)][j]+=d[mask][i]*A[i][j]
 return h
def triangles(A):
 return {(i,j,k):A[i][j]*A[j][k]*A[k][i]+A[i][k]*A[k][j]*A[j][i] for i,j,k in combinations(range(len(A)),3)}
def pow2(e):return F(2**e) if e>=0 else F(1,2**(-e))
def c(n):return (24-pow2(3-n)*(n*n+2*n+3))/(n*(n-1))
def beta(n):return (24-pow2(4-n)*(n*n-n+3))/(n*(n-1))

def gamma(n,k):return 4*(n-k+2-F(n+k+2,2**k))/F(n*(n-1))
def gap(A):
 u,v=uv(A);return all_h(A)[-1]-u-abs(v)
@functools.lru_cache(None)
def binary_gap(n,code):return gap(matrix(n,[F(code>>i&1) for i in range(n*(n-1)//2)]))
def expectation(A):
 n=len(A);edges=list(combinations(range(n),2));fixed=0;uncertain=[]
 for e,(i,j) in enumerate(edges):
  p=A[i][j]
  if p==1:fixed|=1<<e
  elif p:uncertain.append((e,p))
 value=F(0)
 for bits in range(2**len(uncertain)):
  code=fixed;prob=F(1)
  for k,(e,p) in enumerate(uncertain):
   if bits>>k&1:code|=1<<e;prob*=p
   else:prob*=1-p
  value+=prob*binary_gap(n,code)
 return value
def uncertainty_cover(A):
 n=len(A);unc=[(i,j) for i,j in combinations(range(n),2) if A[i][j] not in (0,1)]
 return min(C.bit_count() for C in range(1<<n) if all((C>>i&1) or (C>>j&1) for i,j in unc))

counts={'matrices':0,'margin_product_checks':0,'quadratic_edge_checks':0,'rational_gap_checks':0,'matching_graphs':0,'matching_groups':0,'baseline_square_checks':0,'near_candidate_checks':0};digest=hashlib.sha256()
def margin_pair(A):
 u,v=uv(A);return [u+v,u-v]
def quadratic_check(A,i,j):
 vals=[]
 for p in [F(0),F(1,2),F(1)]:
  B=[r[:] for r in A];B[i][j]=p;B[j][i]=1-p;vals.append(margin_pair(B))
 vs=[k for k in range(len(A)) if k not in [i,j]];R=[[A[k][l] for l in vs] for k in vs];mr=margin_pair(R)
 f=margin_pair(A);p=A[i][j];q=p*(1-p)
 for sign in range(2):
  c=vals[1][sign];a=(vals[2][sign]+vals[0][sign])/2-c;b=(vals[2][sign]-vals[0][sign])/2
  require(a==mr[sign]/4,'quadratic coefficient pair-deletion normalization')
  if a:require(4*a*(c-a)>=b*b,'completed-square floor')
  else:require(b==0,'zero quadratic coefficient has no linear part')
  L=(1-p)*vals[0][sign]+p*vals[2][sign]
  require(L==f[sign]+4*q*a,'endpoint interpolation exact identity')
  require(f[sign]*(1+4*q)>=L,'multiplicative one-edge lower bound')
 counts['quadratic_edge_checks']+=1
def audit(A,all_edges=True):
 n=len(A);margins=margin_pair(A);qs=[A[i][j]*A[j][i] for i,j in combinations(range(n),2)];P=F(1)
 for q in qs:P/=1+4*q
 require(min(margins)>=P,'multiplicative full-margin bound');counts['margin_product_checks']+=1
 es=list(combinations(range(n),2))
 for i,j in (es if all_edges else es[:3]):quadratic_check(A,i,j)
 if n>=3:
  D=all_h(A)[-1]-uv(A)[0]-abs(uv(A)[1]);E=expectation(A)
  for order in [sorted(qs),sorted(qs,reverse=True)]:
   correction=F(0);later_product=F(1)
   for q in reversed(order):correction+=q*later_product;later_product/=1+4*q
   require(D>=E+correction,'exact rational telescoping gap bound');counts['rational_gap_checks']+=1
  if n>=4:
   r=n//2;kappa=F(4,r)*(1-F(3,4)**r-(F(1,4)**r if n%2==0 else 0))
   tau=sum(triangles(A).values());qdict={e:q for e,q in zip(es,qs) if q};V=sum(qs)
   degrees=[sum(q for e,q in qdict.items() if i in e) for i in range(n)]
   W=sum(d*d for d in degrees)/2-sum(q*q for q in qs)
   require(E>=tau>=W,'binary baseline triangle/wedge reserve')
   pending=qdict.copy();M=F(0)
   while pending:
    e=max(pending,key=lambda e:(pending[e],e));M+=pending[e]
    pending={f:q for f,q in pending.items() if not(set(e)&set(f))}
   require(D>=E+kappa*M,'baseline-preserving matching support bound')
   require(D>=kappa*V-F(n,2)*kappa*kappa,'baseline square-completion bound')
   counts['baseline_square_checks']+=1
   if D<=kappa*V and V:
    require(E<=2*n*kappa*kappa,'near-candidate binary baseline')
    require(V-M<=2*n*kappa,'near-candidate matching mass loss')
    counts['near_candidate_checks']+=1
  digest.update(str((n,min(margins),P,D,E)).encode())
 counts['matrices']+=1
for n in range(1,5):
 for p in product([F(0),F(1,2),F(1)],repeat=n*(n-1)//2):audit(matrix(n,p))
rng=random.Random(138608)
for n in range(5,9):
 for trial in range(8):
  es=list(combinations(range(n),2));p=[F(rng.randrange(2)) for e in es]
  for k in rng.sample(range(len(es)),6):
   den=rng.choice([3,5,7,11]);p[k]=F(rng.randrange(1,den),den)
  audit(matrix(n,p),False)
# Fair matching controls attain the candidate constant exactly.
for n in range(5,9):
 es=list(combinations(range(n),2));pairs={(i,i+1) for i in range(0,n-1,2)}
 audit(matrix(n,[F(1,2) if e in pairs else F(1) for e in es]),False)
# Empty determinant margins and n2 gap are explicit excluded boundaries.
require(margin_pair([])==[2,0],'empty-boundary coordinates')
require(min(margin_pair([]))<1,'empty margin product theorem must fail')
A=matrix(2,[F(1,2)]);require(gap(A)==0,'n2 zero gap')
require(margin_pair(A)==[F(1,2),F(1)],'n2 nonempty product sharp control')
for n in range(2,26):
 for trial in range(16):
  es=list(combinations(range(n),2));q={e:F(rng.randrange(9),rng.choice([3,7,11])) for e in es};q={e:w for e,w in q.items() if w}
  V=sum(q.values(),F(0));d=[sum(w for e,w in q.items() if i in e) for i in range(n)]
  W=sum(x*x for x in d)/2-sum(w*w for w in q.values());remaining=q.copy();M=F(0);groups=[]
  while remaining:
   e=max(remaining,key=lambda e:(remaining[e],e));w=remaining[e];M+=w
   assigned={f:z for f,z in remaining.items() if f!=e and(set(f)&set(e))};S=sum(assigned.values(),F(0))
   aa=sum(z for f,z in assigned.items() if e[0] in f);bb=S-aa
   local=w*S+(aa*aa+bb*bb-sum(z*z for z in assigned.values()))/2
   require(all(z<=w for z in assigned.values()),'greedy dominance')
   require(local>=S*S/4+w*S/2,'local wedge lower bound');groups.append(S)
   remaining={f:z for f,z in remaining.items() if not(set(f)&set(e))}
  require(W>=sum(S*S for S in groups)/4,'distinct edge-group wedge pairs')
  require(2*n*W>=(V-M)**2,'global greedy matching capture')
  counts['matching_graphs']+=1;counts['matching_groups']+=len(groups)
  digest.update(str((n,V,W,M)).encode())
print(json.dumps({'status':'PASS','scope':'all ternary n1..4;8 six-fractional-edge rational cases and one sharp fair-matching control each n5..8;16 nonnegative rational graph weights per order n2..25','counts':counts,'explicit_boundaries':['r1 margins1','r2 fair edge margins1/2,1','empty margins2,0 excludes empty product claim','n2 gap0 excludes nonlinear positive correction'],'semantic_sha256':digest.hexdigest()},indent=2,sort_keys=True))
