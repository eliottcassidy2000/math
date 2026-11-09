import json,math
from collections import defaultdict
from pathlib import Path
root=Path(__file__).parent
banks=defaultdict(set)
for line in (root/'MINIMAL_SEEDS.txt').read_text().splitlines():
 h,q,a,g=map(int,line.split());banks[h,q].add((a,g))
def need(b,s):
 if not b:raise RuntimeError(s)
expected={(33,660):176,(33,1320):704,(33,2640):2112,(341,5456):492,(341,10912):1968}
need({k:len(v) for k,v in banks.items()}==expected,'independent bank differs')
rows=[]
for (h,q),seeds in sorted(banks.items()):
 order=next(d for d in range(1,12) if pow(2,d,h)==1);need(order==10,'order')
 counts=[0,0];returns=[0,0];gains=[set(),set()]
 for a,g in sorted(seeds):
  need(math.gcd(h,math.gcd(a,g))==1,'content')
  for extra in (0,1):
   p=(a%2)^extra;start=a,g,p;stack=[(a,g,p,0,0)]
   while stack:
    aa,gg,pp,d,wt=stack.pop();counts[extra]+=1
    if d==10:
     if (aa,gg,pp)==start:returns[extra]+=1;gains[extra].add(wt)
     continue
    # Enumerate all four bit pairs and solve the raw interface equations;
    # do not import or assume the producer's parity-reduced edge routine.
    for x in (0,1):
     for y in (0,1):
      na=2*aa+h*y-q*x;ng2=gg+q*y-h*x
      if ng2%2:continue
      ng=ng2//2
      if not(-h<na<q and -h<ng<q):continue
      nw=aa*y+h*pp*y-x*ng-h*x*(ng%2)
      stack.append((na,ng,y,d+1,wt+nw))
 need(returns==[0,0],'closed short return')
 rows.append(dict(h=h,q=q,order2=order,physical_AG_seeds=len(seeds),compatible_p_partial_walks=counts[0],opposite_p_partial_walks=counts[1],closed_walks=returns))
out=dict(status='PASS',AG_seeds=sum(len(x) for x in banks.values()),augmented_seeds=2*sum(len(x) for x in banks.values()),compatible_partial_walks=sum(r['compatible_p_partial_walks'] for r in rows),both_p_partial_walks=sum(r['compatible_p_partial_walks']+r['opposite_p_partial_walks'] for r in rows),groups=rows)
(root/'MINIMAL_WALKS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
