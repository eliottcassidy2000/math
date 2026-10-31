"""Independent audit: no producer imports; explicit partition enumeration + DP."""
from itertools import product
from fractions import Fraction
from pathlib import Path
import json, hashlib
import sympy as s
ROOT=Path(__file__).resolve().parents[1]/'code'
states=list(product(range(2),range(3),range(3),range(2)))
index={t:i for i,t in enumerate(states)}
# Local matrices from append-to-run specification, then tensor sum.
matrices=[]
for letter in range(4):
 blocks=[]
 for a,b in ((0,1),(1,0)):
  aa=(letter>>a)&1;bb=(letter>>b)&1
  p=s.zeros(3);q=s.zeros(3);r=s.zeros(2)
  for j in range(3):
   p[2 if aa else 1,j]=1
   if j or aa:q[2 if bb else 1,j]=1
  for j in range(2):
   if j or bb:r[1,j]=1
  blocks.append(s.kronecker_product(p,s.eye(3),s.eye(2))+s.kronecker_product(s.eye(3),q,s.eye(2))+s.kronecker_product(s.eye(3),s.eye(3),r))
 matrices.append(s.diag(*blocks))
output=s.Matrix([[int(p!=2 and (q==1 or (o==0 and q==0))) for o,p,q,r in states]])
initial=[int(p==q==r==0) for o,p,q,r in states]
sparse=[[(i,j,int(M[i,j])) for i in range(36) for j in range(36) if M[i,j]] for M in matrices]
def step(v,c):
 out=[0]*36
 for i,j,z in sparse[c]:out[i]+=z*v[j]
 return out
def transfer(word):
 v=initial
 for c in word:v=step(v,c)
 return v
def count(v):return sum(int(output[i])*v[i] for i in range(36))
def partitions(word):
 # Enumerate complete assignments, not transfer transitions.
 v=[0]*36
 for o,(a,b) in enumerate(((0,1),(1,0))):
  for colors in product(range(3),repeat=len(word)):
   runs=[[i for i,c in enumerate(colors) if c==k] for k in range(3)]
   P,Q,R=runs
   if Q and not (word[Q[0]]>>a)&1:continue
   if R and not (word[R[0]]>>b)&1:continue
   p=0 if not P else 2 if (word[P[-1]]>>a)&1 else 1
   q=0 if not Q else 2 if (word[Q[-1]]>>b)&1 else 1
   v[index[o,p,q,int(bool(R))]]+=1
 return v
def adjacency_dp(word):
 n=len(word)+2
 def edge(i,j):
  if i<2 and j<2:return i<j
  if i>=2 and j>=2:return i<j
  if i<2:return bool((word[j-2]>>i)&1)
  return not bool((word[i-2]>>j)&1)
 D={(1<<i,i):1 for i in range(n)}
 for mask in range(1,1<<n):
  for end in range(n):
   v=D.get((mask,end),0)
   if v:
    for nxt in range(n):
     if not mask>>nxt&1 and edge(end,nxt):
      key=(mask|(1<<nxt),nxt);D[key]=D.get(key,0)+v
 return sum(D.get(((1<<n)-1,i),0) for i in range(n))
checks=0
for length in range(6):
 for word in product(range(4),repeat=length):
  v=transfer(word)
  if v!=partitions(word):raise RuntimeError(('partition',word))
  if count(v)!=adjacency_dp(word):raise RuntimeError(('adjacency',word))
  for o in range(2):
   if sum(v[index[o,p,0,0]] for p in range(3))!=1:raise RuntimeError('unit')
  checks+=1
# Sole-outneighbor decomposition: enumerate one-ear paths as TWO increasing runs.
sole_checks=0
for n in range(5):
 for word in product((0,2),repeat=n):
  S=t=V=0
  for mask in range(1<<n):
   P=[i for i in range(n) if not mask>>i&1]
   R=[i for i in range(n) if mask>>i&1]
   right=not R or word[R[0]]==2
   left=not P or word[P[-1]]==0
   V+=right
   if left and right:
    S+=1
    t+=not R
  if adjacency_dp(word)!=S-t+V:raise RuntimeError('sole decomposition')
  sole_checks+=1
  for k in range(1,4):
   if adjacency_dp(word+(2,)*k)!=S+2**k*V:raise RuntimeError('sole repeat')
   sole_checks+=1
  # A final 2,0 produces a directed triangle in the one-ear subtournament.
  ww=word+(2,0);ss=0
  for mask in range(1<<len(ww)):
   P=[i for i in range(len(ww)) if not mask>>i&1]
   R=[i for i in range(len(ww)) if mask>>i&1]
   ss+=(not P or ww[P[-1]]==0) and (not R or ww[R[0]]==2)
  if ss<3:raise RuntimeError('triangle')
  sole_checks+=1
# Independent certificate multiplication; coefficients and target are external data only.
def word_matrix(word):
 A=s.eye(36)
 for c in word:A=matrices[c]*A
 return A
certificates=json.loads((ROOT/'ternary_certificates.json').read_text())
cert_results=[]
for cert in certificates:
 bank=cert['bank'];B=cert['radix']
 expected=[[0],[1],[2]] if B==3 else [[0,0],[0,1],[0,2],[0,3],[1,0],[1,2],[2,0],[2,1],[3,0]]
 if bank!=expected:raise RuntimeError('bank')
 L=len(bank[0])
 for term in cert['terms']:
  w=term['word']
  if len(w)%L or any(w[i:i+L] not in bank for i in range(0,len(w),L)):raise RuntimeError('illegal descendant')
 defect=output*(sum((word_matrix(w) for w in bank),s.zeros(36))-B*B*s.eye(36))
 result=s.zeros(1,36)
 for term in cert['terms']:
  result+=s.Rational(term['coefficient'])*defect*word_matrix(term['word'])
 target=s.Matrix([[int(o==0 and q==r==0) for o,p,q,r in states]])
 if result!=target or list(target)!=cert['target']:raise RuntimeError('certificate')
 if sum(Fraction(t['coefficient']) for t in cert['terms']):raise RuntimeError('sum')
 rows=defect;dims=[1]
 while True:
  expanded=rows.col_join(s.Matrix.vstack(*(rows*word_matrix(w) for w in bank)))
  ids=expanded.T.rref()[1];new=expanded[list(ids),:];dims.append(new.rows)
  if new.rows==rows.rows:break
  rows=new
 if dims!=cert['closure_dimensions']:raise RuntimeError(('rank',dims))
 if rows.col_join(output).rank()!=rows.rows:raise RuntimeError('output outside defect closure')
 second=s.Matrix([[int(o==1 and q==r==0) for o,p,q,r in states]])
 if rows.col_join(second).rank()!=rows.rows:raise RuntimeError('second unit outside defect closure')
 cert_results.append(dict(radix=B,coordinate_identities=36,closure_dimensions=dims,coefficient_sum=0))
# Common basis, all four generators and all 16 two-letter products.
P=s.Matrix([[1,0,0],[-1,1,1],[0,-1,0]])
Q=s.Matrix([[0,0,1],[1,1,0],[-1,0,0]])
R=s.Matrix([[0,1],[1,0]])
C=s.kronecker_product(P,Q,R);C=s.diag(C,C);Ci=C.inv()
diagonal=None
for M in matrices:
 T=Ci*M*C
 if any(T[i,j] for i in range(36) for j in range(i)):raise RuntimeError('flag')
 d=list(T.diagonal())
 if diagonal is None:diagonal=d
 if diagonal!=d:raise RuntimeError('diagonal')
if [diagonal.count(k) for k in range(4)]!=[8,16,10,2]:raise RuntimeError('multiplicity')
for w in product(range(4),repeat=2):
 T=Ci*word_matrix(w)*C
 if list(T.diagonal())!=[x*x for x in diagonal]:raise RuntimeError('product')
# Exact annihilator checks are stronger than numerical sequence fits.
I=s.eye(36);A=word_matrix((0,1));v=s.Matrix(initial)
if (A-4*I)*(A-I)**2*A*v!=s.zeros(36,1):raise RuntimeError('Jordan control')
for k in range(1,4):
 if (output*A**k*v)[0]!=s.Rational(4**k+6*k-1,9):raise RuntimeError('Jordan initials')
A=matrices[0]
if A*(A-I)*(A-2*I)*(A-3*I)!=s.zeros(36):raise RuntimeError('minimal polynomial')
v=s.Matrix(transfer((2,3,3)))
for k in range(1,4):
 if (output*A**k*v)[0]!=14*3**k-17:raise RuntimeError('cancellation')
# Recompute future-observable row closures, including exact invariant kernels.
obs=[]
for bank in ((0,1),(0,2),(1,2)):
 rows=output
 while True:
  X=rows.col_join(rows*matrices[bank[0]]).col_join(rows*matrices[bank[1]])
  ids=X.T.rref()[1];new=X[list(ids),:]
  if new.rows==rows.rows:break
  rows=new
 obs.append(dict(bank=bank,dimension=rows.rows))
if [v['dimension'] for v in obs]!=[10,10,18]:raise RuntimeError('observable')
print(json.dumps(dict(partition_and_adjacency_words=checks,sole_outneighbor_checks=sole_checks,max_length=5,certificates=cert_results,triangular_generators=4,triangular_products=16,diagonal_multiplicities=[8,16,10,2],observable=obs,hostile_spectral_controls=2,certificate_sha256=hashlib.sha256((ROOT/'ternary_certificates.json').read_bytes()).hexdigest()),sort_keys=True))
