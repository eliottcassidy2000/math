"""Independent Session116 arithmetic audit; standard library, exception gates only."""
from itertools import product, combinations
from collections import defaultdict, Counter
from math import gcd
import json
checks=0
def gate(p,s):
 global checks
 checks+=1
 if not p: raise RuntimeError(s)
def mul(a,b):
 n=len(a)
 return tuple(sum(a[i]*b[i^k] for i in range(n))%2 for k in range(n))
def norm_defect(S,x,y):
 n=len(S)
 z=[sum(S[i][i^k]*x[i]*y[i^k] for i in range(n)) for k in range(n)]
 return sum(t*t for t in z)-sum(t*t for t in x)*sum(t*t for t in y)
def expanded(S):
 # Multiply out each output coordinate squared directly, without rectangle formula.
 n=len(S);D=defaultdict(int)
 for k in range(n):
  for i,j in product(range(n), repeat=2):
   D[(tuple(sorted((i,j))),tuple(sorted((i^k,j^k))))]+=S[i][i^k]*S[j][j^k]
 for i,j in product(range(n),repeat=2):D[((i,i),(j,j))]-=1
 return {m:c for m,c in D.items() if c}
counts=Counter()
for bits in product((-1,1),repeat=9):
 S=[[1]*4 for _ in range(4)]
 for (i,j),s in zip(product(range(1,4),repeat=2),bits):S[i][j]=s
 D=expanded(S)
 gate(set(D.values())<= {-4,4},'expanded coefficient set')
 g=0
 for (rows,cols),c in D.items():
  gate(rows[0]!=rows[1] and cols[0]!=cols[1],'square-free mixed monomial')
  x=tuple(int(i in rows) for i in range(4));y=tuple(int(i in cols) for i in range(4))
  gate(norm_defect(S,x,y)==c,'every sparse coefficient witness')
  g=gcd(g,c)
 gate(g in (0,4),'coefficient-value ideal');counts[g]+=1
 gate(not D or all(norm_defect(S,tuple(int(i in r) for i in range(4)),tuple(int(i in c) for i in range(4)))%4==0 for r,c in product(combinations(range(4),2),repeat=2)),'direct norm residues')
gate(counts=={0:2,4:510},'normalized dimension4 census')
# Independently expand all eight-dimensional line tables and the fixed sigma bank.
from math import prod, isqrt
lines=sorted({tuple(sorted((i,j,i^j))) for i in range(1,8) for j in range(i+1,8)})
bank_good=0
for bits in product((0,1),repeat=7):
 S=[[1]*8 for _ in range(8)]
 for i in range(1,8):S[i][i]=-1
 for (i,j,k),b in zip(lines,bits):
  for u,v in ((i,j),(j,k),(k,i)):S[u][v]=(-1)**b;S[v][u]=-S[u][v]
 D=expanded(S);sigma=prod(p+1 for p,b in zip((19,2,59,29,5,89,7),bits) if b)
 square=sigma%30==0 and isqrt(sigma//30)**2==sigma//30
 gate((not D)==square,'fixed arithmetic bank polynomial classification')
 gate(set(D.values())<= {-4,4},'eight-dimensional polynomial coefficients')
 bank_good+=not D
gate(bank_good==16,'bank coherent census')
# Session112 forced-line convention: derive each of the five equations directly.
equations=[]
for t,a,b in ((2,5,4),(8,1,7),(11,4,1),(12,6,5),(13,4,3)):
 variables=set();constant=0
 for i,j in ((a,b),(a^t,b^t),(a,b^t),(a^t,b)):
  line=tuple(sorted((i,j,i^j)));u,v,w=line
  constant^=int((i,j) not in ((u,v),(v,w),(w,u)))
  if line in variables:variables.remove(line)
  else:variables.add(line)
 equations.append((variables,1^constant))
all_variables=set();rhs=0
for variables,r in equations:all_variables.symmetric_difference_update(variables);rhs^=r
gate(not all_variables and rhs==1,'five-equation exact contradiction')
# All pairs including nonunits, independently represented as coefficient tuples.
G={};pair_counts={}
for d in range(4):
 n=2**d;vectors=list(product((0,1),repeat=n));g=0
 for a in vectors:
  gate(mul(a,a)==(sum(a)%2,)+(0,)*(n-1),'unit square formula')
  for b in vectors:
   z=mul(a,b);delta=sum(z)-sum(a)*sum(b)
   gate(delta%4==0,'universal weight mod4')
   if sum(a)%2 and sum(b)%2:g=gcd(g,delta)
 G[d]=g;pair_counts[d]=len(vectors)**2
 gate(g=={0:0,1:0,2:8,3:4}[d],'unit defect ideal')
def phase(d,v,r,c):
 return tuple(c^((r&(v^x)).bit_count()%2)^int(x==v) for x in range(2**d))
subgroups={}
for d in (3,4,5):
 n=2**d;P=list(product(range(n),range(n),range(2)));F={p:phase(d,*p) for p in P}
 gate(len(set(F.values()))==2*n*n,'unique subgroup parameters')
 # Exhaustively certify multiplication by generators; associativity of convolution
 # then certifies entire parameter-additive map, independently of pair enumeration.
 generators=[(1<<i,0,0) for i in range(d)]+[(0,1<<i,0) for i in range(d)]+[(0,0,1)]
 for (v,r,c),a in F.items():
  gate(sum(a)%4==(1 if c==0 else 3),'weight character')
  for w,s,e in generators:
   gate(mul(a,F[w,s,e])==F[v^w,r^s,c^e],'generator multiplication')
  if d in (3,4):
   peak=max(abs(sum((-1)**(a[x]^((u&x).bit_count()%2)) for x in range(n))) for u in range(n))
   gate(peak==n-2,'Walsh extremality')
 subgroups[d]={'size':len(F),'index':2**(n-2*d-2),'generators':len(generators)}
 gate(sum(1 for p in F if p[2]==0)==n*n,'character fiber')
# Higher-d hostile is an embedded 3-coordinate support; no global enumeration.
for d in range(3,8):
 n=2**d;a=tuple(int(x in (0,1,2)) for x in range(n));b=tuple(int(x in (0,1,4)) for x in range(n))
 z=mul(a,b)
 gate({i for i,t in enumerate(z) if t}=={2,3,4,5,6},'exact hostile support')
 gate(sum(z)-sum(a)*sum(b)==-4,'sharp hostile defect')
# Sharp d2 failure of the affine square-zero assertion.
a=tuple(x&1 for x in range(4));b=tuple((x>>1)&1 for x in range(4))
gate(mul(a,b)==(1,1,1,1),'affine product boundary')
print(json.dumps({'status':'PASS','checks':checks,'dim4_ideal_census':dict(counts),'all_pair_counts':pair_counts,'unit_defect_gcd':G,'subgroups':subgroups},sort_keys=True,indent=2))
