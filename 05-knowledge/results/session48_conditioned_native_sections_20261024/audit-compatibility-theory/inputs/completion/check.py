#!/usr/bin/env python3
"""Independent brute host/observer enumeration versus block-completion formulas."""
import argparse,itertools,json,math
from collections import Counter
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--output',type=Path);args=p.parse_args();checks=0

def req(x):
 global checks
 checks+=1
 if not x:raise ValueError('exact control '+str(checks))
def blocks(n,paths):
 succ={};pred={}
 for P in paths:
  for u,v in zip(P,P[1:]):
   if u in succ and succ[u]!=v:return None
   if v in pred and pred[v]!=u:return None
   succ[u]=v;pred[v]=u
 seen=set();out=[]
 for v in range(n):
  if v in pred:continue
  P=[];u=v
  while u not in seen:
   seen.add(u);P.append(u)
   if u not in succ:break
   u=succ[u]
  out.append(tuple(P))
 return out if len(seen)==n else None

def formula(n,paths,e,pins=None):
 B=blocks(n,paths)
 if B is None:return Counter(),[]
 c=len(B);free=e-c+1;req(free>=0)
 orders=[];poly=Counter()
 for order in itertools.permutations(B):
  P=tuple(x for b in order for x in b)
  if pins and any(P[k]!=v for v,k in pins.items()):continue
  orders.append(P);a=sum(b[-1]<d[0] for b,d in zip(order,order[1:]))
  for j in range(free+1):poly[a+j]+=math.comb(free,j)
 return poly,orders

def partial_cases(n,cells):
 pairs=list(itertools.combinations(range(n),2));covered=set(e for C in cells for e in itertools.combinations(C,2));cv=sorted(covered);free=[e for e in pairs if e not in covered]
 for code in range(1<<len(cv)):
  base=[[0]*n for _ in range(n)]
  for j,(u,v) in enumerate(cv):base[u][v]=(code>>j)&1;base[v][u]=1-base[u][v]
  menus=[[P for P in itertools.permutations(C) if all(base[u][v] for u,v in zip(P,P[1:]))] for C in cells]
  hosts=[]
  for mask in range(1<<len(free)):
   S=[r[:] for r in base]
   for j,(u,v) in enumerate(free):S[u][v]=(mask>>j)&1;S[v][u]=1-S[u][v]
   HP=[P for P in itertools.permutations(range(n)) if all(S[u][v] for u,v in zip(P,P[1:]))]
   hosts.append((mask.bit_count(),HP))
  for paths in itertools.product(*menus):
   positive=set(tuple(sorted((u,v))) for P in paths for u,v in zip(P,P[1:]));negative=covered-positive
   # If another cell declares a nonconsecutive pair positive, both requirements apply.
   negative=set(e for C,P in zip(cells,paths) for e in itertools.combinations(C,2) if e not in {tuple(sorted(x)) for x in zip(P,P[1:])})
   actual=[]
   for weight,HP in hosts:
    for P in HP:
     adj={tuple(sorted(x)) for x in zip(P,P[1:])}
     if positive<=adj and not negative&adj:actual.append((weight,P))
   expected,orders=formula(n,paths,len(free));got=Counter(w for w,P in actual);req(got==expected)
   B=blocks(n,paths)
   if B is None:req(not actual)
   else:req(len(actual)==math.factorial(len(B))*(1<<(len(free)-len(B)+1)))
   # Converse reverses all local observer orders and every uncovered source arc.
   rev,_=formula(n,[P[::-1] for P in paths],len(free));req(rev==Counter({len(free)-k:v for k,v in expected.items()}))
   if code==0:
    for v in range(n):
     for k in range(n):
      exp,_=formula(n,paths,len(free),{v:k});req(exp==Counter(w for w,P in actual if P[k]==v))
   yield len(actual),bool(B),len(free)

stats={'local_instances':0,'feasible_local_instances':0,'marked_completions_counted':0}
for n in range(1,5):
 cells=[C for k in (2,3) for C in itertools.combinations(range(n),k)]
 families=[()]+[(C,) for C in cells]+list(itertools.combinations(cells,2))
 if n>=3:families.append(tuple((0,i) for i in range(1,n)))
 for family in families:
  for count,feasible,e in partial_cases(n,family):
   stats['local_instances']+=1;stats['feasible_local_instances']+=feasible;stats['marked_completions_counted']+=count
# Empty-cell Eulerian identity, independent full source and Hamiltonian-order census.
euler=[]
for n in range(1,6):
 pairs=list(itertools.combinations(range(n),2));Pall=list(itertools.permutations(range(n)));actual=Counter()
 for code in range(1<<len(pairs)):
  S=[[0]*n for _ in range(n)]
  for j,(u,v) in enumerate(pairs):S[u][v]=(code>>j)&1;S[v][u]=1-S[u][v]
  for P in Pall:
   if all(S[u][v] for u,v in zip(P,P[1:])):actual[code.bit_count()]+=1
 A=Counter(sum(u<v for u,v in zip(P,P[1:])) for P in Pall);f=len(pairs)-n+1;pred=Counter()
 for a,c in A.items():
  for j in range(f+1):pred[a+j]+=c*math.comb(f,j)
 req(actual==pred);req(sum(actual.values())==math.factorial(n)*(1<<f))
 euler.append({'n':n,'Eulerian_ascents':[A[i] for i in range(n)],'marked_tournament_polynomial':[actual[i] for i in range(len(pairs)+1)]})
# Branch and cycle obstruction controls; pairwise coherent sources alone do not cure them.
req(blocks(4,[(0,1),(0,2),(0,3)]) is None)
req(blocks(3,[(0,1),(1,2),(2,0)]) is None)
# Two disjoint fixed two-owner observers:16 marked completions, not16 distinct hosts necessarily.
poly,orders=formula(4,[(0,1),(2,3)],4);req(sum(poly.values())==16 and len(orders)==2)
out={'status':'PASS','checks':checks,**stats,'Eulerian_controls':euler,'scope':'All coherent one/two-cell families of sizes2/3 on at most4 labelled owners, star controls, selected clock pins, converse, and empty-cell full labelled tournament census n<=5. Counts marked source/observer pairs, not isomorphism classes.'}
s=json.dumps(out,indent=2,sort_keys=True)+'\n'
if args.output:args.output.write_text(s)
else:print(s,end='')
