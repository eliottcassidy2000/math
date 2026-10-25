#!/usr/bin/env python3
"""Independent finite controls for the Session56 niche atlas."""
import itertools,json
from fractions import Fraction

def chk(p,msg):
 if not p:raise RuntimeError(msg)
def tree(a,l,r):
 k=min(range(l,r+1),key=lambda i:(a[i],i));p={k:None}
 for lo,hi in [(l,k-1),(k+1,r)]:
  if lo<=hi:
   q=tree(a,lo,hi);q[next(i for i in q if q[i] is None)]=k;p.update(q)
 return p
def lca(p,x,y):
 s=set()
 while x is not None:s.add(x);x=p[x]
 while y not in s:y=p[y]
 return y
na=nq=0
for n in range(1,8):
 for a in itertools.product(range(3),repeat=n):
  p=tree(a,0,n-1);na+=1
  for l in range(n):
   for r in range(l,n):
    chk(lca(p,l,r)==min(range(l,r+1),key=lambda i:(a[i],i)),'stable RMQ');nq+=1
chk(tree([1,1],0,1)==tree([1,2],0,1),'tie-loss hostile')
chk(tree([1,1],0,1)==tree([5,5],0,1) and min([1,1])!=min([5,5]),'value-loss hostile')
def meet(a):
 z=3
 for x in a:z&=x
 return z
nm=0
for n in range(1,6):
 for a in itertools.product(range(4),repeat=n):
  for l in range(n):
   for r in range(l,n):
    k=1<<((r-l+1).bit_length()-1)
    chk(meet(a[l:r+1])==meet(a[l:l+k])&meet(a[r-k+1:r+1]),'mask interval');nm+=1
chk(sum([1,2,3])!=sum([1,2])+sum([2,3]),'overlap count hostile')
def rel(blocks):return frozenset((i,j) for s in blocks for i in s for j in s if i<j)
A=[{0,4,8},{2,3,7},{1,5,6}];B=[{2,4,6},{0,5,7},{1,3,8}]
chk(not(rel(A)&rel(B)) and rel(A)&rel(A)==rel(A),'CRT alignment hostile')
p12=rel([{0,1},{2}]);p13=rel([{0,2},{1}]);bot=frozenset()
chk(bot<=p12 and bot<=p13 and not p12<=p13 and not p13<=p12,'partition non-tree hostile')
def count(edges,k):return sum(all(c[i]!=c[j] for i,j in edges) for c in itertools.product(range(k),repeat=3))
ports={}
for k in range(2,6):
 x=count([(0,1)],k);y=count([(0,1),(0,2)],k)
 chk((x,y)==(k*k*(k-1),k*(k-1)**2),'port counts');ports[k]=[x,y]
for r in range(3,9):
 chk([min(abs(t)*(r-1)+s*abs(1-t)+(r-s)*abs(t) for t in range(-12,13)) for s in range(r+1)]==list(range(r))+[r-1],'higher arity')
for x in range(64):
 for y in range(x+1,64):chk(((x^y)&-(x^y))==((y-x)&-(y-x)),'low-bit trie valuation')
# Direct Fibonacci charge construction: atoms (1,0),(0,1), then sum.
f=[1,2];v=[(1,0),(0,1)]
while f[-1]<=100:f.append(f[-1]+f[-2]);v.append(tuple(v[-1][j]+v[-2][j] for j in range(2)))
def q(n):
 a=b=0
 for z,(x,y) in reversed(list(zip(f,v))):
  if z<=n:n-=z;a+=x;b+=y
 return a,b
def color(n):return tuple(t%2 for t in q(n))
def xor(a,b):return tuple(a[i]^b[i] for i in range(2))
chk(xor(color(1),color(4))!=xor(color(2),color(3)),'Zeckendorf overlap')
carrycases=0
for a in range(50):
 for b in range(50):
  qa,qb,qc=q(a),q(b),q(a+b);delta=qa[1]+qb[1]-qc[1]
  chk(delta in {-1,0,1} and tuple(qa[j]+qb[j]-qc[j] for j in range(2))==(-2*delta,delta),'carry');carrycases+=1
steps=0
for x in range(1,101):
 z=x;a=Fraction(1);b=Fraction(0)
 for j in range(12):
  if z%2:a*=3;b=3*b+1;z=3*z+1
  else:a/=2;b/=2;z//=2
  chk(a*x+b==z,'Collatz affine itinerary');steps+=1
chk(((3*1+1)&-(3*1+1)).bit_length()-1==2 and ((3*5+1)&-(3*5+1)).bit_length()-1==4,'valuation hostile')
def mm(a,b):return [[min(a[i][k]+b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
mat=[[[a,b],[c,d]] for a,b,c,d in itertools.product(range(2),repeat=4)]
for a,b,c in itertools.product(mat,repeat=3):chk(mm(mm(a,b),c)==mm(a,mm(b,c)),'min-plus associativity')
chk(mm([[1,1],[1,1]],[[1,1],[1,1]])==[[2,2],[2,2]],'min-plus nonidempotence')
print(json.dumps({'status':'FINITE-EXACT PASS','cartesian_arrays':na,'cartesian_interval_queries':nq,'mask_range_queries':nm,'named_port_counts':ports,'higher_arity_controls':{'r':[3,8],'tested_integer_coefficients':[-12,12]},'low_bit_trie_pairs':2016,'Zeckendorf_carry_pairs':carrycases,'Collatz_affine_steps':steps,'minplus_associativity_triples':len(mat)**3,'limits':['No native tournament realization','No persistent structure implemented or benchmarked','No arithmetic global closure']},sort_keys=True,indent=2))
