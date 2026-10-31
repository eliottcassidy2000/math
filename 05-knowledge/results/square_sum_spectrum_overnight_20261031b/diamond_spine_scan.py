#!/usr/bin/env python3
"""All-word singleton/strong-four-block image below an explicit ceiling.

The q cutoff retains the ordinal-product tail, unlike an unsafe direct transfer
of singleton-only pruning. No word-depth bound is used.
"""
import argparse
from collections import deque
import hashlib
import json
from diamond_profiles import profiles,apply,single
from audit_feedback import require


def main():
 p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=10000);args=p.parse_args()
 N=args.limit;require(N>=1,'positive ceiling');blocks=profiles()
 roots=[(1,0,0,1),(1,1,1,1)]
 seen=set(roots);queue=deque(roots);support=set();overflow_outputs=set();fringe=0
 def emit(F):
  while F<=N:support.add(F);F*=5
 while queue:
  source=queue.popleft();emit(source[0])
  children=[single(source,b) for b in (0,1)]+[apply(source,b) for b in blocks]
  for state in children:
   F,q,g,h=state
   require(F>=source[0] and q>=source[1] and h>=1 and F>=h and 0<=g<=h,'actual-state monotonicity')
   if F>N:continue
   emit(F)
   if q>N:
    fringe+=1;overflow_outputs.add(F);continue
   if state not in seen:seen.add(state);queue.append(state)
 require(7 not in support and 21 not in support,'known permanent hostile')
 missing=[h for h in range(1,N+1,2) if h not in support and h not in (7,21)]
 bitmap=bytes(int(h in support) for h in range(N+1))
 print(json.dumps({'status':'complete','limit':N,'source_roots':roots,'states':len(seen),
                   'distinct_counts':len(support),'allowed_missing':missing,'q_overflow_edges':fringe,
                   'distinct_q_overflow_output_counts':len(overflow_outputs),
                   'support_sha256':hashlib.sha256(bitmap).hexdigest(),
                   'scope':'Exact all-word image for these two roots and eighteen generators, including all powers-of-five ordinal-product tails; no unrestricted absence claim.'},indent=2,sort_keys=True))


if __name__=='__main__':main()
