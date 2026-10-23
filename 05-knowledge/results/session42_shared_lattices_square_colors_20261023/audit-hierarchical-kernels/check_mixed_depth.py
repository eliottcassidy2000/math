"""Independent exhaustive W2^3 mixed-prefix counts and owner set cover."""
from itertools import product
from collections import Counter
import json
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
# Chart tuple (outer_rotation, three inner rotations) has exactly81 values.
charts=list(product(range(3),repeat=4));edges=[(0,1),(1,2),(0,2)]
def agreement(x,y):return 2 if x==y else (1 if x[0]==y[0] else 0)
hist=Counter((agreement(a,b),agreement(b,c),agreement(a,c)) for a,b,c in product(charts,repeat=3))
def components(depths,threshold):
 labels=list(range(3))
 for (i,j),l in zip(edges,depths):
  if l>=threshold:
   old=labels[j];new=labels[i];labels=[new if x==old else x for x in labels]
 ren={};return [ren.setdefault(x,len(ren)) for x in labels]
rows=[]
for depths in product(range(3),repeat=3):
 C0=components(depths,1);C1=components(depths,2);c0=len(set(C0));c1=len(set(C1))
 exact=sum(v for h,v in hist.items() if all(x>=y for x,y in zip(h,depths)))
 expected=3**(c0+3*c1)
 need(exact==expected,'mixed exact chart count')
 # Actual owner (cell, first_digit, last_digit) forces the root variable of
 # its G0 component and the matching prefix variable of its G1 component.
 owners={((1<<C0[cell])|(1<<(c0+digit*c1+C1[cell]))) for cell in range(3) for digit in range(3)}
 target=(1<<(c0+3*c1))-1;dp={0:0}
 for owner in owners:
  for mask,cost in list(dp.items()):
   nxt=mask|owner;dp[nxt]=min(dp.get(nxt,99),cost+1)
 need(dp[target]==3*c1,'mixed exact owner cover')
 rows.append({'edge_depths':depths,'group_order':exact,'owner_base':dp[target],'component_counts':[c0,c1]})
print(json.dumps({'status':'PASS','checks':checks,'enumerated_chart_triples':81**3,'all_mixed_depth_patterns':27,'scope':'declared W2^3 chart carrier; base independently checked via exact parameter-fixing set cover','rows':rows},indent=2))
