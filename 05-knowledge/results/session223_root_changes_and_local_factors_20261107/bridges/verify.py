import itertools,json,math

def step(s,c):
 D,E,U,d=s
 return {'O':(D,0,2*U,d),'Z':(D-E+U,d,2*U-d,d),'J':(D-E+U,D-E,D+U,D)}[c]
def run(w,s=(1,0,2,1)):
 for c in w:s=step(s,c)
 return s
def F(s):D,E,U,d=s;return D-E+U
def require(x):
 if not x:raise RuntimeError('check failed')
states=0
for n in range(8):
 for word in itertools.product('OZJ',repeat=n):
  s=run(word);states+=1
  require(run('JZ',s)==run('OJ',s))
  for r in range(6):require(run('J'+'Z'*r,s)==run('O'*r+'J',s))
for a in range(1,41):
 for b in range(1,41):
  s=run('Z'*a+'O'*b);t=run('Z'*b+'O'*a)
  require(s==(2**a+1,0,2**b*(2**a+1),1))
  require(F(s)==F(t)==(2**a+1)*(2**b+1))
  require(F(step(s,'O'))-F(step(t,'O'))==2**b-2**a)
edges=[]
for a in range(1,201):
 for b in range(a+1,201):
  n=2**a+2**b+2
  if math.isqrt(n)**2==n:edges.append([a,b])
require(edges==[[1,5]])
# Check the reflection formula independently on arbitrary root lists and start labels.
checks=0
for k in range(1,7):
 for roots in itertools.product(range(1,5),repeat=k):
  for x in [-3,0,7]:
   y=x
   for r in roots:y=r*r-y
   formula=(-1)**k*x+sum((-1)**(k-1-i)*r*r for i,r in enumerate(roots))
   require(y==formula);checks+=1
# A closed reflection walk need not be simple.
x=1;walk=[x]
for r in [2,2,2,2]:x=r*r-x;walk.append(x)
require(walk==[1,3,1,3,1])
print(json.dumps({'reachable_state_checks':states,'two_run_pairs':1600,'dyadic_adjacency_exponents_through':200,'dyadic_edges':edges,'reflection_checks':checks,'closure_not_simple':walk},sort_keys=True,indent=2))
