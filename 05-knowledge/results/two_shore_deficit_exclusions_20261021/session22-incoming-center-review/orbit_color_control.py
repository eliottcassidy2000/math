import json
from itertools import product
checks=0
cases=0
def ck(ok,why):
 global checks
 checks+=1
 if not ok:raise ValueError(why)
def partitions(n):
 def rec(a):
  if len(a)==n:
   yield tuple(a);return
  for k in range(max(a)+2):yield from rec(a+[k])
 yield from rec([0])
def run(word,color):
 best=cur=0
 for x in word:
  cur=cur+1 if x==color else 0
  best=max(best,cur)
 return best
for q in range(5,9):
 ck(q-(2 if q%2==0 else 1)+4>=8,'four additional black orbits')
 for word in partitions(q):
  colors=range(max(word)+1)
  counts={a:word.count(a) for a in colors}
  for a,c in counts.items():
   if c>4:continue
   # Every actual anti-map image is one of these equal-size classes.
   images=[b for b,d in counts.items() if d==c]
   L=4 if c==1 else 5
   for m in (5,7):
    full=word*m
    bound={b:max(run(full[i:]+full[:i],b) for i in range(q)) for b in images}
    ck(all(v<=c for v in bound.values()),'all cuts cyclic color-run bound')
    prefix=(a,)*L
    ck(run(prefix,a)>bound[a],'constructed original run separates cuts')
    for b in images:
     # Restriction of any induced color permutation to this repeated run.
     mate=tuple(b for x in reversed(prefix))
     ck(run(mate,b)>bound[b],'reversed and permuted run separates cuts')
    cases+=1
rows=[(25,5,5),(35,5,7),(35,7,5)]
ck(all(n==q*m and 5<=q<=8 and m>=5 for n,q,m in rows),'remaining large rows')
print(json.dumps({'checks':checks,'color_class_cases':cases,'q':[5,6,7,8],'m':[5,7],'necessary_vertex_transitive_orders':[25,35],'scope':'abstract partition/color-permutation control, not a tournament census'},sort_keys=True))
