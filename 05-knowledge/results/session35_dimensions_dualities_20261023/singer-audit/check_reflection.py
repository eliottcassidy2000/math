from itertools import permutations,product,combinations
from collections import Counter
from pathlib import Path
import json

def req(x,m):
 if not x:raise RuntimeError(m)
def run(p):
 Q={x*x%p for x in range(1,p)}
 def arc(a,b):return (b-a)%p in Q
 C=[tuple(x) for x in combinations(range(p),3) if all(sum(arc(a,b) for b in x if a!=b)==1 for a in x)]
 m=(p-1)//2;degrees=Counter();fixeds=Counter()
 # Normalize reflection center to zero; all other centers are translations.
 for order in permutations(range(1,m+1)):
  for signs in product((1,-1),repeat=m):
   left=[a*s%p for a,s in zip(order,signs)]
   P=left+[0]+[-a%p for a in left[::-1]]
   if not all(arc(a,b) for a,b in zip(P,P[1:])):continue
   pe={frozenset(e) for e in zip(P,P[1:])}
   A=[x for x in C if not any(frozenset(e) in pe for e in combinations(x,2))]
   fixed=sum(set(x)=={-a%p for a in x} for x in A)
   req(fixed==((p-3)//2 if p%8==3 else 0),'fixed count')
   req(len(A)%2==0,'even admissible degree')
   degrees[len(A)]+=1;fixeds[fixed]+=1
 return {'p':p,'reflection_center':0,'reflection_symmetric_HP_count':sum(degrees.values()),'admissible_degree_histogram':dict(degrees),'fixed_admissible_counts':dict(fixeds)}
P=[6,4,5,3,7,0,1,10,8,9,2];p=11;Q={i*i%p for i in range(1,p)}
def arc(a,b):return (b-a)%p in Q
req(all(arc(a,b) for a,b in zip(P,P[1:])),'witness HP')
blue=all(arc(P[a],P[b])==arc(P[10-b],P[10-a]) for a,b in combinations(range(11),2))
pe={frozenset(e) for e in zip(P,P[1:])}
A=[C for C in combinations(range(11),3) if all(sum(arc(a,b) for b in C if a!=b)==1 for a in C) and not any(frozenset(e) in pe for e in combinations(C,2))]
req(not blue and len(A)==32,'black even witness')
out={'status':'PASS','reflection_controls':[run(p) for p in (3,7,11)],'black_even_witness':{'p':11,'path':P,'blue':blue,'admissible_degree':len(A)}}
Path(__file__).with_name('reflection_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
