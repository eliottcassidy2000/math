from collections import defaultdict,Counter
from fractions import Fraction as Q
from itertools import product

def move(s,c,q):
 D,E,U,d=s
 return tuple(x%q for x in ((D,0,2*U,d) if c=='O' else (D-E+U,d,2*U-d,d) if c=='Z' else (D-E+U,D-E,D+U,D)))
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def compute(q):
 # After one O/Z letter E is 0 or d. This whole universe is 2-in,2-out regular.
 universe={(D,E,U,d) for D,U,d in product(range(q),repeat=3) for E in {0,d}}
 indeg=Counter(move(s,c,q) for s in universe for c in 'OZ')
 need(set(indeg)==universe and all(x==2 for x in indeg.values()),'regularity')
 classes=[];owner={}
 for s in sorted(universe):
  if s in owner:continue
  i=len(classes);todo=[s];owner[s]=i
  for v in todo:
   for c in 'OZ':
    z=move(v,c,q)
    if z not in owner:owner[z]=i;todo.append(z)
    need(owner[z]==i,'closed component')
  classes.append(todo)
 # Balanced digraph + forward closure => each such component strongly connected.
 # Independently validate reverse reachability, not only rely on this argument.
 for i,C in enumerate(classes):
  rev=defaultdict(list)
  for s in C:
   for c in 'OZ':rev[move(s,c,q)].append(s)
  seen={C[0]};todo=[C[0]]
  for s in todo:
   for v in rev[s]:
    if v not in seen:seen.add(v);todo.append(v)
  need(len(seen)==len(C),'reverse connectivity')
 mass=defaultdict(Q)
 for c in 'OZ':mass[owner[move((1,0,2,1),c,q)]]+=Q(1,2)
 # Initial class must be aperiodic, witnessed by a self-loop.
 need(len(mass)==1,'initial single class')
 C=classes[next(iter(mass))]
 need(any(move(s,c,q)==s for s in C for c in 'OZ'),'initial aperiodicity')
 for k in range(3):
  out=[Q(0)]*q
  for i,p in mass.items():
   for s in classes[i]:out[(s[0]-s[1]+s[2])%q]+=p/len(classes[i])
  need(sum(out)==1,'probability total')
  print('q=',q,'J=',k,'probability=',','.join(map(str,out)))
  vals=[]
  for v in range(round(__import__('math').log(q,3))):
   inds=[r for r in range(q) if r%(3**v)==0 and r%(3**(v+1))!=0]
   vals.append(sum(out[r] for r in inds)*q/len(inds))
  vals.append(q*out[0])
  print('q=',q,'J=',k,'valuation-local-factors=',','.join(map(str,vals)))
  if k==2:
   precision=0;power=q
   while power>1:precision+=1;power//=3
   predicted=[Q(47,54)]+[Q(157,108)-Q(1,3**v)*(Q(v,6)+Q(85,144)) for v in range(1,precision)]
   predicted.append(Q(157,108)-Q(1,3**precision)*(Q(precision,8)+Q(11,24)))
   need(vals==predicted,'analytic all-k factor')

  new=defaultdict(Q)
  for i,p in mass.items():
   cnt=Counter(owner[move(move(s,'J',q),c,q)] for s in classes[i] for c in 'OZ')
   for j,n in cnt.items():new[j]+=p*Q(n,2*len(classes[i]))
  mass=new
 print('PASS q=',q,'states=',len(universe),'classes=',len(classes),'all classes reverse-connected; initial aperiodic')
for q in (3,9,27):compute(q)
