"""Exact checks for the periodic obstruction. Python 3; no third-party packages."""
import itertools,json,sys,random
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'ribbon'))
from transfer import automaton,advance,graph,held_karp,count

D=json.loads((HERE/'representation.json').read_text());B=D['basis'];Rs=D['matrices'];f=D['output']
N=len(B);R=len(f);states,mat=automaton(3)
assert D['w']==3 and N==119 and R==62 and len(Rs)==8

def row_times(v,A):return [sum(x*A[i][j] for i,x in enumerate(v) if x) for j in range(len(A[0]))]
def ident():return [[int(i==j) for j in range(R)]for i in range(R)]
def product(word,p=None):
 A=ident()
 for a in word:A=mul_letter(A,a,p)
 return A
SP=[[[ (j,x) for j,x in enumerate(row) if x]for row in A]for A in Rs]
def mul_letter(A,a,p=None):
 C=[[0]*R for _ in range(R)]
 for i,row in enumerate(A):
  out=C[i]
  for k,x in enumerate(row):
   if x:
    for j,y in SP[a][k]:out[j]+=x*y
  if p is not None:
   for j in range(R):out[j]%=p
 return C

def determinant_minus(A,c,p):
 A=[r[:]for r in A]
 for i in range(R):A[i][i]=(A[i][i]-c)%p
 d=1
 for k in range(R):
  r=next((r for r in range(k,R) if A[r][k]),None)
  if r is None:return 0
  if r!=k:A[k],A[r]=A[r],A[k];d=-d
  x=A[k][k];iv=pow(x,-1,p);d=d*x%p
  for i in range(k+1,R):
   if A[i][k]:
    z=A[i][k]*iv%p;A[i][k]=0
    for j in range(k+1,R):A[i][j]=(A[i][j]-z*A[k][j])%p
 return d%p

# Intertwining is sufficient for every future output, without a minimality claim.
for i,st in enumerate(states):
 assert sum(B[i][j]*f[j]for j in range(R))==int(len(st)==1)
 for a in range(8):
  lhs=[sum(c*B[k][j]for k,c in mat[i][a])for j in range(R)]
  assert lhs==row_times(B[i],Rs[a]),(i,a)

stats={str(l):{'classes':0,'first_nonzero_modulus':{str(p):0 for p in [3,5,7,11,13]}}for l in range(1,4)}
def walk(word,A):
 if word and word==min(word[k:]+word[:k]for k in range(len(word))):
  z=stats[str(len(word))];z['classes']+=1
  for p in [3,5,7,11,13]:
   P=A if p==3 else product(word,p)
   if determinant_minus(P,2,p):z['first_nonzero_modulus'][str(p)]+=1;break
  else:raise AssertionError(('eigenvalue 2 unresolved',word))
 if len(word)<3:
  for a in range(8):walk(word+(a,),mul_letter(A,a,3))
walk((),ident())
assert [stats[str(l)]['classes']for l in range(1,4)]==[8,36,176]

# The forward-sink reset on every nonempty frontier basis state.
q=None
for i,st in enumerate(states):
 if not st:continue
 v=tuple(int(j==i)for j in range(N))
 for _ in range(3):v=advance(v,mat,0)
 if len(st)==1:
  if q is None:q=v
  assert v==q
 else:assert not any(v)

# Counterexample: boundary-deletion gcd is not a future divisibility invariant.
word=list(map(int,'0001500'));a=graph(word,3);deletions=[]
for mask in range(8):
 keep=[i for i in range(7)if i<4 or not(mask>>(i-4)&1)]
 deletions.append(held_karp([[a[i][j]for j in keep]for i in keep]))
assert deletions==[15,5,15,5,15,5,15,5]
assert held_karp(graph(word+[1],3))==91

# The seductive stronger no-eigenvalue-2 assertion is false.
E=json.loads((HERE/'eigenvalue2-counterexample.json').read_text());A=product(E['word']);z=E['right_eigenvector']
assert any(z) and [sum(row[j]*z[j]for j in range(R))for row in A]==[2*x for x in z]
v=(1,)+(0,)*(N-1);ys=[];hs=[]
for k in range(3):
 for a in E['word']:v=advance(v,mat,a)
 y=row_times(v,B);ys.append(y);hs.append(sum(y[j]*f[j]for j in range(R)))
assert hs==[9,543,40461]
assert all(c-77*b+150*a==0 for a,b,c in zip(*ys))
# Consequently H(W^k)=(7*75^k+66*2^k)/73 for every k>=1.

# Every possible overlapping-triangle geometry forces the R_(2m+1) graph.
overlap_checks=0
for w in range(1,8):
 for d1 in range(1,w+1):
  for d2 in range(1,w+1):
   L=d1+d2
   if L<=w:continue
   for m in range(1,7):
    vertices=sorted({j*L+x for j in range(m)for x in [0,d1,L]})
    back=set()
    for j in range(m):back|={(j*L+d1,j*L),(j*L+L,j*L+d1)}
    assert len(vertices)==2*m+1
    for j in range(len(vertices)):
     for i in range(j):
      if j==i+1:assert (vertices[j],vertices[i])in back and vertices[j]-vertices[i]<=w
      else:assert vertices[j]-vertices[i]>=L>w
    assert count([1]*len(vertices),1)>=3**m
    overlap_checks+=1

# Eventual transitive suffixes, including hostile prefix orientations.
rng=random.Random(197);tail_checks=0
for L in range(1,5):
 words=itertools.product(range(8),repeat=L)
 for word in words:
  suffix=list(word)*((L+6+L-1)//L);a=graph(suffix,3)
  cyclic=False
  for t in itertools.combinations(range(len(a)),3):
   if all(sum(a[i][j]for j in t)==1 for i in t):cyclic=True;break
  if cyclic:continue
  prefix=[rng.randrange(8)for _ in range(6)]
  k0=(6+L-1)//L;vals=[count(prefix+list(word)*k,3)for k in range(k0,k0+5)]
  assert len(set(vals))==1,(word,vals)
  tail_checks+=1

result={'short_block_certificates':stats,'intertwining_rows':119*8,'reset_rows':118,
 'boundary_deletion_counterexample':{'word':'0001500','deletion_counts':deletions,'append_1_count':91},
 'eigenvalue2_counterexample':{'word':'001476','initial_counts':hs,'formula':'(7*75^k+66*2^k)/73'},
 'overlap_geometry_checks':overlap_checks,'transitive_tail_prefix_checks':tail_checks}
(HERE/'VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
