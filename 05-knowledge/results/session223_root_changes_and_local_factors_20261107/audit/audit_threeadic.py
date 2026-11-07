from itertools import product
from fractions import Fraction as Q
from collections import defaultdict,Counter
import json
S=list(product(range(3),repeat=4)); idx={s:i for i,s in enumerate(S)}
def move(s,c):
 D,E,U,d=s
 if c=='O':v=(D,0,2*U,d)
 elif c=='Z':v=(D-E+U,d,2*U-d,d)
 else:v=(D-E+U,D-E,D+U,D)
 return tuple(x%3 for x in v)
def require(x):
 if not x:raise RuntimeError('audit check failed')
adj={s:{move(s,c) for c in 'OZ'} for s in S}
# Reachability, independently determine communicating classes and recurrence.
reach={}
for s in S:
 seen={s};todo=[s]
 while todo:
  for v in adj[todo.pop()]:
   if v not in seen:seen.add(v);todo.append(v)
 reach[s]=seen
remaining=set(S);classes=[]
while remaining:
 s=min(remaining);cl={u for u in reach[s] if s in reach[u]};classes.append(cl);remaining-=cl
closed=[cl for cl in classes if all(adj[s]<=cl for s in cl)]
which={s:i for i,cl in enumerate(closed) for s in cl}
# Uniform stationary mass verified from incoming edge multiplicities, not assumed.
for cl in closed:
 incoming=Counter(move(s,c) for s in cl for c in 'OZ')
 require(all(incoming[s]==2 for s in cl))
for s in S:require(all(move(s,c) in which for c in 'OZ'))
# Direct Cesaro projection from computed recurrent classes and one-step entry.
def project(v):
 out=defaultdict(Q)
 for s,p in v.items():
  for c in 'OZ':
   cl=closed[which[move(s,c)]]
   for u in cl:out[u]+=p/(2*len(cl))
 return dict(out)
def J(v):
 out=defaultdict(Q)
 for s,p in v.items():out[move(s,'J')]+=p
 return out
def residues(v):
 out=[Q(0)]*3
 for (D,E,U,d),p in v.items():out[(D-E+U)%3]+=p
 return out
t=(1,0,2,1);v=project({t:Q(1)});laws=[]
for j in range(4):
 laws.append([str(x) for x in residues(v)])
 v=project(J(v))
require(laws[0]==['1/3']*3)
require(laws[1]==['41/108','8/27','35/108'])
require(laws[2]==['34/81','47/162','47/162'])
# Integer DP independently counts full words, with J-count as an extra coordinate.
dp={(t,0):1};samples={}
for n in range(1,201):
 nxt=defaultdict(int)
 for (s,k),count in dp.items():
  for c in 'OZJ':
   kk=k+(c=='J')
   if kk<=2:nxt[(move(s,c),kk)]+=count
 dp=nxt
 if n in [2,3,5,10,25,50,100,200]:
  out=[0]*3
  for ((D,E,U,d),k),count in dp.items():
   if k==2:out[(D-E+U)%3]+=count
  total=n*(n-1)//2*2**(n-2)
  require(sum(out)==total)
  target=[Q(34,81),Q(47,162),Q(47,162)]
  samples[str(n)]={'probabilities':[float(Q(x,total)) for x in out],'L_scaled_errors':[float(n*(Q(x,total)-y)) for x,y in zip(out,target)]}
report={'closed_class_sizes':sorted(map(len,closed)),'transient_states':81-len(which),'projected_laws_J0_to_J3':laws,'fixed_length_checks':samples}
def apply(v,c):
 out=defaultdict(Q)
 for s,p in v.items():out[move(s,c)]+=p
 return dict(out)
v0=project({t:Q(1)});v1=project(J(v0))
conditional={str(k):{c:[str(x) for x in residues(project(J(apply(v,c))))] for c in 'OZ'} for k,v in [(1,v0),(2,v1)]}
report['conditional_last_bit']=conditional
# Literal physical check of the ratio partition through source length 9.
def physical(s,c):
 D,E,U,d=s
 if c=='O':return (D,0,2*U,d)
 if c=='Z':return (D-E+U,d,2*U-d,d)
 return (D-E+U,D-E,D+U,D)
ratio_checks=0
for n in range(3,10):
 for w in product('OZJ',repeat=n):
  if w.count('J')!=2:continue
  j2=max(i for i,c in enumerate(w) if c=='J')
  if w[j2-1]=='J':continue
  state=(1,0,2,1)
  for c in w[:j2]:state=physical(state,c)
  D,E,U,d=state
  require((2*D<=U)==(w[j2-1]=='O'))
  ratio_checks+=1
report['literal_ratio_partition_checks']=ratio_checks
print(json.dumps(report,indent=2,sort_keys=True))
