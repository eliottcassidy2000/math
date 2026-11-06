"""Independent rank-pencil audit, no author-code imports or assert statements."""
import hashlib,itertools,json,math
from pathlib import Path
from fractions import Fraction

T=(1,0,2,1);B=(3,1,3,1);AUX=(1,0,2,0)
M={
 'O':((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1)),
 'Z':((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1)),
 'J':((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))}
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
P=(T,(1,0,4,1),B,(5,2,6,3))
def need(p,m):
 if not p:raise RuntimeError(m)
def act(c,v):return tuple(sum(x*y for x,y in zip(r,v))for r in M[c])
def word(s,w):
 for c in w:s=act(c,s)
 return s
def dagger(w):return w[::-1].translate(str.maketrans('OZJ','ZOJ'))
def cov(v):return tuple(sum(v[i]*Q[i][j]for i in range(4))for j in range(4))
def dot(v,w):return sum(x*y for x,y in zip(v,w))
def K(v,w):return dot(cov(v),w)
def positive(v):
 d,e,u,z=v;return d-e-z,u-d+e-z,e,z-e

def rank2(w,q):
 # Exact rational elimination, rather than the author's minor guard.
 i=next(j for j,x in enumerate(w)if x)
 r=Fraction(q[i],w[i])
 return any(Fraction(y)!=r*x for x,y in zip(w,q))
def determinant(A):
 A=[[Fraction(x)for x in row]for row in A];sgn=1
 for j in range(len(A)):
  i=next((i for i in range(j,len(A))if A[i][j]),None)
  if i is None:return 0
  if i!=j:A[i],A[j]=A[j],A[i];sgn=-sgn
  for k in range(j+1,len(A)):
   s=A[k][j]/A[j][j]
   A[k]=[a-s*b for a,b in zip(A[k],A[j])]
 return sgn*math.prod(A[i][i]for i in range(len(A)))

def bank(H,paired=False,prune=False):
 initial=(B,AUX)if paired else T
 height=lambda v:v[0][2]if paired else v[2]
 out={initial:(0,'')};buckets=[set()for _ in range(H+1)];buckets[height(initial)].add(initial)
 for u in range(H+1):
  for v in sorted(buckets[u]):
   cost,w=out[v]
   for c in 'JZO':
    z=tuple(act(c,x)for x in v)if paired else act(c,v);h=height(z)
    need(h>u,'strict monotone physical cover')
    if h>H or(prune and not rank2(*z)):continue
    recipe=(cost+int(c=='J'),w+c)
    if z not in out:buckets[h].add(z);out[z]=recipe
    elif recipe<out[z]:out[z]=recipe
 return out

def prime(n):return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))

def put(D,k,c):D[k]=min(D.get(k,10**9),c)
def images(S,R,N):
 f={};g={};tests=0
 for v,(a,_)in S.items():
  vq=cov(v)
  for(w,q),(b,_)in R.items():
   m=dot(vq,w);h=dot(vq,q)
   need(3<=h<m<2*h and m%2==h%2==1,'joint physicality and parity')
   if m>N:continue
   tests+=1
   if prime(m):put(f,m,a+b);put(g,m,a+b)
   power=2
   while m+2*(power-1)*h<=N:
    n=m+2*(power-1)*h
    if prime(n):put(f,n,a+b)
    power*=2
 return f,g,tests

def graph(w):
 arcs={(0,1)};n=2;r=0
 for c in w:
  for v in range(n):arcs.add((n,v)if c!='O'and v==r else(v,n))
  if c=='J':r=n
  n+=1
 return n,arcs

def literal(w):
 n,arcs=graph(w);hp=0;hc=0;steps=0
 def visit(start,last,used):
  nonlocal hp,hc,steps
  steps+=1
  if len(used)==n:
   hp+=1
   if start==0 and(last,0)in arcs:hc+=1
   return
  for x in range(n):
   if x not in used and(last,x)in arcs:visit(start,x,used|{x})
 for s in range(n):visit(s,s,{s})
 return dict(word=w,vertices=n,paths=hp,cycles=hc,partial_directed_paths=steps)

N=5001;H=math.isqrt(4*(N-1))
need(determinant(P)==-4 and determinant(Q)==-1,'exact basis and bilinear determinants')
PQ=tuple(cov(v)for v in P)
need(PQ==((1,-1,1,0),(1,-1,3,-2),(1,0,2,0),(3,-1,3,0))and determinant(PQ)==4,'invertible four-probe readout')
# The positive-cone generators in raw coordinates.
R=((1,0,1,0),(0,0,1,0),(2,1,2,1),(1,0,2,1))
positive_matrix=[[K(a,b)for b in R]for a in R]
need(positive_matrix==[[1,1,2,2],[1,1,1,1],[2,1,4,3],[2,1,3,3]],'positive Gram matrix')
need(all(K(a,b)>=3*a[3]*b[3]for a in R for b in R),'cone factor bound on generating rays')
S=bank(H);Rall=bank(H,True);Rp=bank(H,True,True);large=bank(N-1)
need({p:r for p,r in Rall.items()if rank2(*p)}==Rp,'rank-pruned full states and costs')
rank_one=0;factor_checks=0;recipe_checks=0;max_graph=0
for(w,q),(cost,beta)in Rall.items():
 need((word(B,beta),word(AUX,beta))==(w,q),'shared boundary word')
 ms=[K(v,w)for v in P];hs=[K(v,q)for v in P]
 deltas=[ms[0]*h-m*hs[0]for m,h in zip(ms[1:],hs[1:])]
 need(rank2(w,q)==any(deltas),'four physical probes exactly detect rank')
 need(all(x%2==0 for x in deltas),'even determinants')
 if rank2(w,q):
  need(any(Fraction(abs(hs[i]*ms[0]-ms[i]*hs[0]),ms[0]*ms[i])>=Fraction(1,6*w[2]**2)for i in range(1,4)),'exact rational slope gap')
  prefix=''
  for c in beta:
   prefix+=c;need(rank2(word(B,prefix),word(AUX,prefix)),'rank-two witness has no collapsed prefix')
 else:
  rank_one+=1;g=math.gcd(*w);p=tuple(x//g for x in w)
  i=next(i for i,x in enumerate(p)if x);a=Fraction(q[i],p[i])
  need(a.denominator==1 and q==tuple(a*x for x in p),'primitive integral multiplier')
  a=int(a)
  need(g%2==1 and g>=3 and g<2*a<2*g and min(positive(p))>=0 and p[3]>=1,'factor and cone bounds')
  need(all(not rank2(act(c,w),act(c,q))for c in 'OZJ'),'collapse is absorbing')
  for v in S:
   z=K(v,p)
   need(z>=3*v[3]*p[3]>=3,'primitive cone factor positive without physicality')
   for k in range(3):
    L=2*((1<<k)-1)
    need(K(v,w)+L*K(v,q)==(g+L*a)*z,'nontrivial exact target product')
   factor_checks+=1
 # Common-word physicality independently reconstructs every state on probes.
 for v,alpha in zip(P,('','O','Z','ZJ')):
  actual=word(T,alpha+dagger(beta))
  need((K(v,w),K(v,q))==(actual[0]+2*actual[2],2*actual[2]-actual[3]),'literal common-word joint reconstruction')
  recipe_checks+=1
f,g,tests=images(S,Rall,N);fp,gp,pruned_tests=images(S,Rp,N)
direct_f={};direct_g={}
for v,(cost,w)in large.items():
 F=v[0]-v[1]+v[2];G=v[0]+2*v[2]
 if 5<=F<=N and prime(F):put(direct_f,F,cost)
 if 5<=G<=N and prime(G):put(direct_g,G,cost)
need(f==fp==direct_f and g==gp==direct_g,'complete prime images and exact minimum J costs')
# Full forward cutoff is justified by F>=U+1 and G>=2U+1.
need(all(v[0]-v[1]>=v[3]>=1 for v in large),'complete direct numeric cutoff')
controls=[]
for alpha,beta,k,n,rk in [('', 'OJOOJ',1,171,1),('', 'ZOJ',2,139,2),('O','',2,51,2),('','',1,11,2),('','OJOJ',0,45,1)]:
 w,q=word(B,beta),word(AUX,beta);v=word(T,alpha)
 need((2 if rank2(w,q)else 1)==rk,'graph rank control')
 m,h=K(v,w),K(v,q)
 need(m+2*((1<<k)-1)*h==n,'graph arithmetic control')
 recipe=alpha+dagger(beta)+'O'+'Z'*k+'J'
 j=literal(recipe);need((j['paths'],j['cycles'])==(n,1),'literal DFS target and C1')
 j.update(rank=rk,m=m,h=h,k=k);controls.append(j)
x=(word(B,'OJOJ'),word(AUX,'OJOJ'));y=(word(B,'ZZOJ'),word(AUX,'ZZOJ'))
need(x[0]==y[0]and not rank2(*x)and rank2(*y),'same-w rank hostile')
print(json.dumps(dict(status='PASS',target_cap=N,forward_height_cap=N-1,forward_profiles=len(large),paired_height_cap=H,
 physical_probe_determinant=int(determinant(P)),measurement_determinant=int(determinant(PQ)),
 physical_profiles=len(S),synchronized_states=len(Rall),rank_one_states=rank_one,pruned_states=len(Rp),
 primitive_factor_checks=factor_checks,physical_recipe_checks=recipe_checks,
 full_pair_tests=tests,pruned_pair_tests=pruned_tests,prime_C1_targets=len(f),prime_balanced_targets=len(g),
 minimum_J_costs_preserved=True,literal_controls=controls,
 scope='Complete finite height-bounded comparison and literal controls; all-order assertions rely on reviewed proofs; no prime coverage theorem.',
 verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
