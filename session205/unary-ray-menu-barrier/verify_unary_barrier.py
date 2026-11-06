"""Exact finite matrix certificate for the unary-ray menu obstruction."""
from collections import Counter
import itertools,json
import sympy as S

def patterns(labels,b):
 if b<0:return set()
 if not labels:return {tuple([(-1,-1)]*q) for q in range(b+1)}
 u,*rest=labels;out=set()
 for p in patterns(tuple(rest),b):
  for q in [p,p+((u,u),),p+((u,-1),)]:out.add(tuple(sorted(q)))
 for p in patterns(tuple(rest),b-1):out.add(tuple(sorted(p+((-1,u),))))
 for v in rest:
  for p in patterns(tuple(x for x in rest if x!=v),b):
   out.add(tuple(sorted(p+((u,v),))));out.add(tuple(sorted(p+((v,u),))))
 return out

def intro(parts,mask):
 out=Counter()
 left=[-1]+[i for i,(_,e) in enumerate(parts) if e==-1 or not(mask>>e&1)]
 right=[-1]+[i for i,(s,_) in enumerate(parts) if s>=0 and mask>>s&1]
 for i in left:
  for j in right:
   if i==j and i!=-1:continue
   q=[p for k,p in enumerate(parts) if k not in (i,j)]
   q.append((parts[i][0] if i!=-1 else -1,parts[j][1] if j!=-1 else -1))
   if sum(a==-1 for a,b in q)<=1:out[tuple(sorted(q))]+=1
 return out

def require(x,msg):
 if not x:raise RuntimeError(msg)
P=sorted(patterns((0,1),1));I={p:i for i,p in enumerate(P)}
require(len(P)==28,'state count')
x=S.zeros(28,1);x[I[((0,1),)]]=1;x[I[((0,0),(1,1))]]=1
ell=S.Matrix([[int(len(p)==1) for p in P]])
A=[];report=[]
t=S.Symbol('lambda')
for a in range(4):
 m=S.zeros(28)
 for j,p in enumerate(P):
  for q,c in intro(p,a).items():m[I[q],j]+=c
 A.append(m)
 cp=m.charpoly(t).as_expr();r=m.rank();r2=(m*m).rank()
 require(S.expand(cp-t**10*(t-1)**12*(t-2)**5*(t-3))==0,'characteristic polynomial')
 require(r==r2==18,'semisimple zero')
 report.append({'letter':a,'rank':r,'rank_squared':r2,'characteristic_polynomial':str(S.factor(cp))})
# One common triangular order certifies spectra for ALL nonempty word products.
order=[0,21,23,7,14,17,20,4,16,22,6,5,13,19,3,24,8,18,2,26,10,12,27,11,15,25,9,1]
require(sorted(order)==list(range(28)),'permutation')
pos={v:i for i,v in enumerate(order)}
for m in A:
 for i in range(28):
  require(m[i,i]==A[0][i,i],'common diagonal')
  for j in range(28):
   if m[i,j] and i!=j:require(pos[j]<pos[i],'common triangular order')
# The unary count follows from the order-two recurrence after its first vertex.
m=A[1];require((m-S.eye(28))*(m-2*S.eye(28))*m*x==S.zeros(28,1),'unary vector recurrence')
require((ell*m*x)[0]==3 and (ell*m*m*x)[0]==5,'unary initial values')
# Independently enumerate literal paths for short words.
def literal(w):
 n=len(w)+2;edge=[[False]*n for _ in range(n)];edge[0][1]=True
 for j,a in enumerate(w,2):
  for i in range(2):edge[j][i]=bool(a>>i&1);edge[i][j]=not edge[j][i]
  for i in range(2,j):edge[i][j]=True
 return sum(all(edge[p[i]][p[i+1]] for i in range(n-1)) for p in itertools.permutations(range(n)))
checks=0
for k in range(5):
 for w in itertools.product(range(4),repeat=k):
  z=x
  for a in w:z=A[a]*z
  require((ell*z)[0]==literal(w),'literal crosscheck');checks+=1
print(json.dumps({'status':'PASS','dimension':28,'matrices':report,'literal_checks':checks,'common_triangular_order':order,'product_spectrum':'0,1,2^k,3^k for every image word length k>=1','unary_formula':'H(1^n)=2^n+1 for n>=1'},indent=2))
