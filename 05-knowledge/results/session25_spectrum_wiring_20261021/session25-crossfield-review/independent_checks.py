"""Independent controls; no producer module imported; exceptions survive python -O."""
import itertools as it,json,math
import sympy as S

def require(ok,label):
 if not ok: raise RuntimeError(label)
B=S.Matrix([[0,-1],[1,1]]);J=S.Matrix([[1,1],[0,-1]]);I=S.eye(2)
require(B**3==-I and B**6==I and J**2==I and J*B*J==B**-1,'group relations')
group={tuple(B**k*J**e) for k in range(6) for e in range(2)}
require(len(group)==12,'order')
require(sum(S.Matrix(2,2,g)*B==B*S.Matrix(2,2,g) for g in group)==6,'centralizer')
x,y=S.symbols('x y',real=True);v=S.Matrix([x,y]);P=S.Matrix.hstack(v,B*v)
require(S.expand(P.det())==x*x+x*y+y*y and P*B==B*P,'transport')
lines=[y,x+2*y,x+y,2*x+y,x,x-y]
for k,line in enumerate(lines):
 require((J-B**k).rank()==1,'six rank-one equations')
 entries=list((J-B**k)*v)
 require(all(S.rem(e,line,x) == 0 if S.degree(line,x)==1 else S.rem(e,line,y)==0 for e in entries),'six lines')
v0=S.Matrix([1,2]);orbit=[B**k*v0 for k in range(6)]
require(J*v0 not in orbit,'fixed reflection hostile')
# Direct literal ten-family rule transcription, then separate boundary simulation.
def old(i,a,b,c):
 if b==0 and c==0:return None
 if b==0 and c==1:return (1,2,2,a+1) if i==0 else (0,1,0,a+3)
 if b==0:return (i,1+i,a+3,c-2)
 if c==0:return (1-i,a+1+i,0,b-1)
 return (i,a+3,b-1,c-1)
def encode(z):
 if z is None:return None
 i,a,b,c=z;return 2*a+i,b,c
def native(z):
 A,b,c=z
 return encode(old(A%2,A//2,b,c))
def first_return(A,c):
 z=(A,0,c);n=0
 while True:
  z=native(z);n+=1
  if z is None or z[1]==0:return z,n

def formula(A,c):
 k,e=divmod(A,2)
 if c==0:return None,1
 if c==1:
  if e:return (2,0,k+3),1
  if k==0:return (14,0,0),3
  return (17,0,k-1),3
 if c<=k+4:return (6*c-7+3*e,0,k-c+4),c
 return (6*k+20+3*e,0,c-k-5),k+4
returns=0
for A,c in it.product(range(211),repeat=2):
 require(first_return(A,c)==formula(A,c),('first return',A,c));returns+=1
require(encode(old(0,0,0,1))==(5,2,1) and encode(old(1,0,0,1))==(2,0,3),'phase hostile')
# Exhaust all labelled simple graphs <=4, exact width from elimination orders.
def treewidth(n,E):
 best=n-1
 for order in it.permutations(range(n)):
  adj=[set() for _ in range(n)]
  for a,b in E:adj[a].add(b);adj[b].add(a)
  width=0
  for a in order:
   ns=list(adj[a]);width=max(width,len(ns))
   for b in ns:
    adj[b].discard(a)
    adj[b].update(c for c in ns if c!=b)
   adj[a].clear()
  best=min(best,width)
 return best

def support(n,E,k):
 ans=set()
 for c in it.product(range(k),repeat=n):
  if all(c[a]!=c[b] for a,b in E):ans.add((c[0],c[1]))
 return ans
widthcases=0
for n in range(2,5):
 pairs=list(it.combinations(range(n),2))
 for bits in it.product(range(2),repeat=len(pairs)):
  E=[e for e,b in zip(pairs,bits) if b];t=treewidth(n,E)
  require(treewidth(n,set(E)|{(0,1)})<=t+1,'edge addition width')
  for k in range(2,5):
   if t<=k-2:
    require(any(a!=b for a,b in support(n,E,k)),'width separator');widthcases+=1
# Glue triangle-edge leaves: absent, unequal edge, flexible length-two path.
# At k>=3 a two-edge path supports all terminal color pairs.
gluings=0
for k in (3,4):
 for types in it.product(range(3),repeat=3):
  E=[];D=[];n=3
  for typ,(a,b) in zip(types,[(0,1),(0,2),(1,2)]):
   if typ==1:E.append((a,b));D.append((a,b))
   elif typ==2:E.extend([(a,n),(n,b)]);n+=1
  require(support(n,E,k)==support(3,D,k),'support reduction');gluings+=1
clocks=[]
for j in range(101):
 m=2+12*j;n=(2**m-1)//3
 require(n%105==1 and 3*n+1==2**m and n%2==1,'clock family')
 clocks.append(m+1)
print(json.dumps({'status':'PASS','symbolic_group_transport':True,'fixed_reflection_hostile':[list(z) for z in orbit],'A9_first_return_cases':returns,'bounded_width_separation_cases':widthcases,'explicit_fragment_support_cases':gluings,'Collatz_clock_cases':len(clocks),'scope':'Algebraic proofs audited separately; finite checks do not prove global halting or graph coloring claims.'},default=int,indent=2))
