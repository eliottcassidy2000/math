from itertools import product
from fractions import Fraction as F
from collections import defaultdict
def require(ok):
 if not ok: raise RuntimeError('audit check failed')
S=list(product(range(3),repeat=4)); idx={s:i for i,s in enumerate(S)}
def step(s,c):
 d,e,u,h=s
 return tuple(v%3 for v in ((d,0,2*u,h) if c==0 else (d-e+u,h,2*u-h,h) if c==1 else (d-e+u,d-e,d+u,d)))
def apply(dist,c):
 z=defaultdict(F)
 for s,p in dist.items():z[step(s,c)]+=p
 return dict(z)
def plus(*terms):
 z=defaultdict(F)
 for a,dist in terms:
  for s,p in dist.items():z[s]+=a*p
 return {s:p for s,p in z.items() if p}
def Pi(s):
 d,e,u,h=s
 if h:return {(x,y*h%3,z,h):F(1,18) for x,y,z in product(range(3),range(2),range(3))}
 if u:return {(x,0,y,0):F(1,6) for x in range(3) for y in [1,2]}
 return plus((F(1,2),{(d,0,0,0):F(1)}),(F(1,2),{((d-e)%3,0,0,0):F(1)}))
def Q(s):
 d,e,u,h=s
 if h or not u:return {}
 return {(x,0,y,0):F(1 if y==u else -1,6) for x in range(3) for y in [1,2]}
def kernel(dist,K):return plus(*[(p,K(s)) for s,p in dist.items()])
def P(dist):return plus((F(1,2),apply(dist,0)),(F(1,2),apply(dist,1)))
for s in S:
 p,q=Pi(s),Q(s)
 require(P(p)==p)
 require(kernel(P({s:F(1)}),Pi)==p)
 require(kernel(p,Pi)==p)
 require(P(q)==plus((-1,q)))
 require(kernel(P({s:F(1)}),Q)==plus((-1,q)))
 require(kernel(q,Q)==q)
 require(not kernel(p,Q) and not kernel(q,Pi))
print('PASS exact Pi/Q projector and P commutation identities on all81 states')
dist=Pi((1,0,2,1))
for k in range(3):
 out=[sum(p for (d,e,u,h),p in dist.items() if (d-e+u)%3==r) for r in range(3)]
 print(k,out)
 dist=kernel(apply(dist,2),Pi)
# Identify closed SCCs independently by transitive closure.
reach={s:{s,step(s,0),step(s,1)} for s in S}
for k in S:
 for s in S:
  if k in reach[s]:reach[s]|=reach[k]
classes=[];done=set()
for s in S:
 if s in done:continue
 c={v for v in reach[s] if s in reach[v]};done|=c
 if all(step(v,i) in c for v in c for i in [0,1]):classes.append(c)
print('Closed classes:',[(len(c),sorted({s[3] for s in c}),sorted({s[2] for s in c})) for c in classes])
require(sorted(map(len,classes))==[1,1,1,6,18,18])
print('PASS recurrent classification')
