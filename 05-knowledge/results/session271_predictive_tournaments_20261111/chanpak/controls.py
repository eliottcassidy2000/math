#!/usr/bin/env python3
"""Small controls for the new support-rank and boundary sandwich only."""
from predict import growth_rank,spine,HERE
import itertools,random,json
rng=random.Random(27132);checks=0
for core in [(2,0),(2,4,1),(6,4,0)]:
 k=len(core)
 supports=[[c for c in range(1<<k) if mask>>c&1] for mask in range(1,1<<(1<<k))]
 for S in supports:
  r,w=growth_rank(core,S);p=w['permutation'];A=[i for i,x in enumerate(w['active_slots']) if x]
  u=[next(c for c in S if not j or c>>p[j-1]&1) for j in A]
  v=[next(c for c in S if not(c>>p[j]&1)) for j in A if j<k]
  for _ in range(3):
   a=[rng.choice(S) for _ in range(rng.randrange(4))];z=[rng.choice(S) for _ in range(5)];b=[rng.choice(S) for _ in range(rng.randrange(4))]
   word=a+u+z+v+b;h=spine.value(word,core)
   assert r**len(z)<=h<=__import__('math').factorial(k)*2**(k+1)*r**len(word),(core,S,word,r,h)
   checks+=1
out={'scope':'new support-rank sandwich; not inherited transfer re-audit','cores':3,'nonempty_support_cases':525,'sandwich_checks':checks}
(HERE/'controls.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
