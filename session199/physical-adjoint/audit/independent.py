"""Independent positive-coordinate audit of the physical adjoint packet."""
from pathlib import Path
from itertools import product, permutations
from math import isqrt
import hashlib, heapq, json, random
ROOT=Path(__file__).resolve().parent.parent
O=((1,0,1,0),(1,2,1,2),(0,0,0,0),(0,0,1,1))
Z=((2,1,1,1),(0,1,0,0),(0,0,1,1),(0,0,0,0))
J=((0,1,0,1),(0,0,0,0),(1,0,1,1),(0,0,1,0))
A={'O':O,'Z':Z,'J':J};sig={'O':'Z','Z':'O','J':'J'}
t=(0,0,0,1);b=(1,0,1,0)
R=((1,1,2,2),(1,1,1,1),(2,1,4,3),(2,1,3,3))
u=(1,1,2,2);d=(0,0,1,1);fr=(2,1,3,3);gr=(3,2,6,5)
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
C=((1,0,2,1),(0,0,1,0),(1,1,2,2),(0,0,1,1))

def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def dot(v,w):return sum(a*b for a,b in zip(v,w))
def mv(M,v):return tuple(dot(row,v)for row in M)
def tr(M):return tuple(zip(*M))
def mm(M,N):return tuple(tuple(dot(row,col)for col in tr(N))for row in M)
def U(v):return dot(u,v)
def K(v,w):return dot(v,mv(R,w))
def state(word,seed=t):
 for c in word:seed=mv(A[c],seed)
 return seed
def dag(word):return ''.join(sig[c]for c in word[::-1])
def bank(cap,seed=t):
 if U(seed)>cap:return {}
 out={seed:(0,'')};heap=[(U(seed),seed)]
 while heap:
  _,v=heapq.heappop(heap);cost,word=out[v]
  for c in A:
   w=mv(A[c],v);h=U(w)
   require(U(v)<h<=2*U(v),'U-step')
   if h>cap:continue
   nc=cost+(c=='J')
   if w not in out:
    out[w]=(nc,word+c);heapq.heappush(heap,(h,w))
   elif nc<out[w][0]:out[w]=(nc,word+c)
 return out

def graph(word):
 out=[{1},set()];root=0
 for c in word:
  n=len(out);out.append(set())
  for i in range(n):
   if c!='O' and i==root:out[n].add(i)
   else:out[i].add(n)
  if c=='J':root=n
 return out

def literal(word):
 g=graph(word);H=C1=0
 for p in permutations(range(len(g))):
  if all(p[i+1]in g[p[i]]for i in range(len(p)-1)):
   H+=1
   if p[0]==0 and 0 in g[p[-1]]:C1+=1
 return H,C1

require(mm(mm(tr(C),Q),C)==R,'raw-positive congruence')
for c in A:require(mm(tr(A[c]),R)==mm(R,A[sig[c]]),'independent adjoint '+c)
require(mv(R,t)==fr and mv(R,b)==gr,'boundary rows')
require(mm(R,J)==mm(mm(tr(C),((1,0,0,0),(0,-1,0,0),(0,0,1,0),(0,0,0,0))),C),'graft factor')
for i in range(4):
 for j in range(4):
  require(u[i]*u[j]>=R[i][j],'upper matrix')
  require(2*R[i][j]>=u[i]*u[j]+2*d[i]*d[j],'lower matrix')
# Bare integer determinant through independent Leibniz expansion.
det=0
for p in permutations(range(4)):
 x=(-1)**sum(p[i]>p[j]for i in range(4)for j in range(i+1,4))
 for i in range(4):x*=Q[i][p[i]]
 det+=x
require(det==-1,'determinant')
checks=0
words=[''.join(w)for k in range(8)for w in product(A,repeat=k)]
rng=random.Random(199)
words += [''.join(rng.choice(tuple(A))for _ in range(50))for _ in range(100)]
for w in words:
 for z,readout in ((t,fr),(b,gr)):
  n=dot(readout,state(w));cap=isqrt(4*(n-1));found=False
  for cut in range(len(w)+1):
   l,r=state(w[:cut]),state(dag(w[cut:]),z)
   require(K(l,r)==n,'every split exact')
   require(n<=U(l)*U(r)<=2*(n-1),'every split product bounds')
   require(w[:cut].count('J')+dag(w[cut:]).count('J')==w.count('J'),'exact J budget')
   found |= max(U(l),U(r))<=cap
  require(found,'two-sided cutoff')
  checks+=1
N=2001;cap=isqrt(4*(N-1));big=bank(N-1);left=bank(cap);right=bank(cap,b)
full={dot(fr,v)for v in big if dot(fr,v)<=N}
bal=[{3}|{dot(gr,v)for v,(cost,_)in big.items()if cost<=r and dot(gr,v)<=N}for r in range(4)]
balall={3}|{dot(gr,v)for v in big if dot(gr,v)<=N}
pairfull=set();pairbal={3};pairbud=[{3}for _ in range(4)];pairs=0
for v,(cv,va)in left.items():
 for w,(cw,wa)in left.items():
  n=K(v,w);pairs+=1
  if n<=N:pairfull.add(n)
 for w,(cw,wa)in right.items():
  n=K(v,w);pairs+=1
  if n<=N:
   pairbal.add(n)
   for r in range(cv+cw,4):pairbud[r].add(n)
require(full==pairfull and balall==pairbal and bal==pairbud,'independent complete finite inverse')
literals=0
for alpha in ('','O','Z','J'):
 for beta in ('','O','Z','J'):
  require(literal(beta+dag(alpha)+'J')==(K(state(alpha),state(beta)),1),'literal full C1')
  require(literal(alpha+dag(beta)+'OJ')==(K(state(alpha),state(beta,b)),1),'literal balanced C1')
  literals+=2
report={'status':'PASS','proof_sha256':hashlib.sha256((ROOT/'PROOF.md').read_bytes()).hexdigest(),
 'candidate_verifier_sha256':hashlib.sha256((ROOT/'verify.py').read_bytes()).hexdigest(),
 'independent_verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'midpoint_word_tests':checks,'finite_target_cap':N,'source_height_cap':cap,
 'forward_profile_count':len(big),'ordinary_bounded_profiles':len(left),'triangle_bounded_profiles':len(right),
 'independent_pair_tests':pairs,'C1_support':len(full),'balanced_support':len(balall),
 'budget_support':[len(s)for s in bal],'literal_permutation_graphs':literals,
 'scope':'Physical adjoint, two-sided cutoff and exact inverse PASS; no universal coverage inference.'}
Path(__file__).with_name('INDEPENDENT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True))
