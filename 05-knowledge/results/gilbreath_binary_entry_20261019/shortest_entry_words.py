#!/usr/bin/env python3
from itertools import product
from collections import deque
import json

def step(x,u):return tuple(abs(a-b) for a,b in zip(x,x[1:]+(u,)))
rows=[]
for M in range(2,6):
 for F in range(1,5):
  states=frozenset(product(range(M+1),repeat=F))
  q=deque([(states,'')]);seen={states}
  while q:
   states,w=q.popleft()
   if all(max(x)<=1 for x in states):break
   for u in(0,1):
    nxt=frozenset(step(x,u) for x in states)
    if nxt not in seen:seen.add(nxt);q.append((nxt,w+str(u)))
  else:raise RuntimeError('No word')
  # Independent direct replay through every original state.
  for x in product(range(M+1),repeat=F):
   for u in w:x=step(x,int(u))
   if max(x)>1:raise RuntimeError((M,F,w,x))
  rows.append(dict(M=M,F=F,shortest_length=len(w),lexicographic_first_word=w,visited_subset_states=len(seen)))
print(json.dumps({'status':'FINITE-EXACT','universe':'M=2..5,F=1..4; complete breadth-first subset automaton without pruning','rows':rows},indent=2))
