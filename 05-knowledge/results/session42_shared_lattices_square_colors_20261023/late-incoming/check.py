import itertools,json

def need(x):
 if not x:raise RuntimeError('failed control')
H={(x,y,z) for x,y,z in itertools.product(range(3),repeat=3) if (x-y-z)%3==0}
need(len(H)==9)
u=(1,1,0);v=(1,0,1)
for g in (u,v):
 for k in range(3):
  x,y,z=(k*t%3 for t in g);need((x,y,z) in H and x-y-z==0)
need(tuple((2*u[i]+v[i])%3 for i in range(3))==(0,2,1))
need((0,2,1) in H and 0-2-1==-3)
for a,b,c in itertools.product(range(3),repeat=3):
 carry=lambda x,y:(x+y)//3
 need(carry(a,b)+carry((a+b)%3,c)==carry(b,c)+carry(a,(b+c)%3))
 need((a+b-(a+b)%3)%9==3*carry(a,b))
need(all(next(k for k in range(1,10) if k*x%9==0)==9 for x in (1,4,7)))
P=list(itertools.permutations(range(3)));e=(0,1,2)
def mul(a,b):return tuple(a[b[i]] for i in range(3))
def inv(a):return tuple(a.index(i) for i in range(3))
def parity(a):return sum(a[i]>a[j] for i in range(3) for j in range(i+1,3))%2
def times(x,y):return (mul(x[0],y[0]),mul(x[1],y[1]))
def inverse(x):return (inv(x[0]),inv(x[1]))
H={(x,y) for x in P for y in P if parity(x)==parity(y)};A=[x for x in P if parity(x)==0]
need(len(H)==18)
products={times((x,x),(e,a)) for x in P for a in A}
need(products==H and len(products)==len(P)*len(A))
for x,y in H:
 a=mul(inv(x),y);need(a in A and times((x,x),(e,a))==(x,y))
r=(e,(1,0,2));C={times(r,h) for h in H}
for x,y,z in itertools.product(C,repeat=3):need(times(times(x,inverse(y)),z) in C)
print(json.dumps({'status':'PASS','generator_hostile_group_size':9,'nonabelian_prefix_solutions':18,'nonabelian_heap_triples':18**3,'scope':'small independent algebra controls only; no incoming native census replay'},indent=2))
