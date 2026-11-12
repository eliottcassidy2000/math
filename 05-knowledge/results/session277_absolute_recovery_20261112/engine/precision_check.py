"""Independent exhaustive chain controls for the precision-debt lemma."""
from itertools import product
from pathlib import Path
import json
count=0
for length in range(1,5):
 for shifts in product(range(-2,3),repeat=length):
  for b in range(1,7):
   B=[0]*(length+1);B[-1]=b
   for i in range(length-1,-1,-1):B[i]=max(0,B[i+1]-shifts[i])
   # The starting integer is divisible by every required prefix power.
   need=max(0,-min([0]+[sum(shifts[:j]) for j in range(1,length+1)]))
   for c in range(17):
    full=c<<need; projected=full%(1<<B[0])
    for i,s in enumerate(shifts):
     if s<0:
      assert full%(1<<(-s))==0 and projected%(1<<(-s))==0
      full>>=-s;projected>>=-s
     else:full<<=s;projected<<=s
     projected%=1<<B[i+1]
    assert projected==full%(1<<b)
    count+=1
# Exact divisibility does not make a right shift well-defined mod 2**b.
for b in range(1,13):
 assert 0%(1<<b)==(1<<b)%(1<<b)
 assert (0//2)%(1<<b)!=((1<<b)//2)%(1<<b)
out=dict(chain_controls=count,unsafe_fixed_precision_hostile_controls=12)
(Path(__file__).parent/'precision_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
