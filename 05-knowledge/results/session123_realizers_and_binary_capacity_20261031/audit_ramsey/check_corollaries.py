"""Independent controls for scalar residuation, hostile saturation and kernel corollaries."""
from itertools import product
from math import inf,comb
import json

def need(v,label):
 if not v:raise RuntimeError(label)
def add(a,b):
 if a==-inf or b==-inf:return -inf
 if a==inf or b==inf:return inf
 return a+b
def resid(d,e):
 if d==-inf or e==inf:return inf
 if e==-inf or d==inf:return -inf
 return e-d
scalars=[-inf,-7,-1,0,2,9,inf]
for h,d,e in product(scalars,repeat=3):need((add(h,d)<=e)==(h<=resid(d,e)),('adjunction',h,d,e))
def coins(cs,n):
 f=[0]+[-inf]*n
 for b in range(1,n+1):f[b]=max(f[b-l]+p for l,p in cs if l<=b)
 return f
f=coins([(1,0),(2,19),(5,50)],110)
g=coins([(1,0),(4,39),(5,50)],110)
fg=[max(f[i]+g[n-i] for i in range(n+1)) for n in range(111)]
bar=lambda v,n:max(v[n+5*t]-50*t for t in range(21))
need([bar(f,0),bar(g,0),bar(f,1),bar(g,1),bar(fg,1)]==[0,0,7,6,8],'arrival hostile')
def N(s,n):return sum((i&s)==s for i in range(n))
def K(a,b):return sum(1<<(i&j).bit_count() for i in range(a) for j in range(b))
maskchecks=0
for s in range(40):
 p=[N(s,n) for n in range(81)]
 for x in range(41):
  for y in range(41):need(p[x+y]>=p[x]+p[y],('mask',s,x,y));maskchecks+=1
for a in range(1,41):
 l=1<<(a-1).bit_length()
 vals=[K(a,b) for b in range(2*l+1)]
 for b in range(1,2*l+1):
  irr=all(vals[b]>vals[x]+vals[b-x] for x in range(1,b))
  need(irr==(b<=l),('generator',a,b))
# Direct Möbius double sum checks inverse without matrix library/elimination.
for n in range(1,21):
 z=lambda i,s:int((i&s)==s)
 zi=lambda i,s:(-1)**(i.bit_count()-s.bit_count()) if (i&s)==s else 0
 for i,j in product(range(n),repeat=2):
  need(sum(z(i,k)*zi(k,j) for k in range(n))==int(i==j),('inverse',n,i,j))
  need(z(i,j)==(comb(i,j)%2 if j<=i else 0),('Pascal',n,i,j))
  need(sum(N(s,i+1)*N(s,j+1) for s in range(n))==K(i+1,j+1),('Gram',n,i,j))
for depth in range(7):
 for i,j in product(range(1<<depth),repeat=2):
  tensor=1
  for bit in range(depth):tensor*=[[1,1],[1,2]][(i>>bit)&1][(j>>bit)&1]
  need(tensor==1<<(i&j).bit_count(),('tensor',depth,i,j))
print(json.dumps({'status':'PASS','scalar_adjunction_cases':len(scalars)**3,'hostile_values':[7,8],'mask_cases':maskchecks,'generator_A_max':40,'gram_inverse_sizes':20,'tensor_depths':7},sort_keys=True))
