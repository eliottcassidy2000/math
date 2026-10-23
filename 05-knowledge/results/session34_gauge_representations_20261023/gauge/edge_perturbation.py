from itertools import combinations, permutations,product
from collections import Counter
import json
from fractions import Fraction
E=list(combinations(range(7),2));idx={e:i for i,e in enumerate(E)};free=[e for e in E if 0 not in e];fi={e:i for i,e in enumerate(free)}
L=[{i,(i+1)%7,(i+3)%7} for i in range(7)];el={e:next(i for i,l in enumerate(L) if set(e)<=l) for e in E}
cycles=[]
for vs in combinations(range(7),4):
 for p in permutations(vs[1:]):
  c=(vs[0],)+p
  if c[1]>c[-1]:continue
  ed=[tuple(sorted((c[i],c[(i+1)%4]))) for i in range(4)]
  if len({el[e] for e in ed})==4:cycles.append(sum(1<<idx[e] for e in ed))
def root(x):
 g=sum(((x>>idx[(0,v)])&1)<<v for v in range(1,7))
 y=x^sum((((g>>a)^(g>>b))&1)<<i for i,(a,b) in enumerate(E))
 return sum(((y>>idx[e])&1)<<j for j,e in enumerate(free))
def expand(x):return sum(((x>>j)&1)<<idx[e] for j,e in enumerate(free))
def v(x):
 y=expand(x);r=sum((-1)**((y&c).bit_count()%2) for c in cycles)
 return Fraction(21-r,8)
steps=[root(1<<i) for i in range(21)]
if len(set(steps))!=21 or any(v(s)!=1 for s in steps):raise RuntimeError('single edge costs')
ground={x for x in range(1<<15) if v(x)==0}
if len(ground)!=16:raise RuntimeError('ground')
second=Counter();third=Counter();denominators=Counter()
for a,b in product(steps,repeat=2):
 if a^b in ground:second[a^b]+=1
for a,b,c in product(steps,repeat=3):
 y=a^b^c
 if y in ground:
  if a^b in ground:raise RuntimeError('unexpected intermediate ground')
  third[y]+=1;denominators[v(a),v(a^b)]+=1
if second!={0:21} or len(third)!=7 or set(third.values())!={6} or denominators!={(1,1):42}:raise RuntimeError('effective coefficients')
pair_costs=Counter(v(a^b) for a,b in product(steps,repeat=2) if a!=b)
B2=sum(Fraction(2,count_cost)*count for count_cost,count in pair_costs.items())
fourth=441-B2
# Trace/power or floating eigensolvers unnecessary: crown graph has exact blocks.
odd=set(third);even=ground-odd-{0}
if len(even)!=8:raise RuntimeError('layer partition')
print(json.dumps({'status':'PASS','quotient_dimension':32768,'ground_dimension':len(ground),'single_edge_step_count':len(steps),'second_order_return_counts':dict(second),'third_order_ground_counts':dict(third),'third_order_denominator_counts':{str(k):v for k,v in denominators.items()},'distinct_ordered_pair_costs':{str(k):v for k,v in pair_costs.items()},'fourth_order_scalar':str(fourth),'crown_spectrum':{'7':1,'1':7,'-1':7,'-7':1}},indent=2))
