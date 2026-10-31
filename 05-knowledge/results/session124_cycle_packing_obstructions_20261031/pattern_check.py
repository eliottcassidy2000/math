"""Exact finite forbidden-subsequence presentation and monotone recognizer."""
from itertools import product,combinations
from collections import Counter
from pathlib import Path
import json
from check import graph,triangles
P=Path(__file__).resolve().parent
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def packing_bad(w):
 ts=triangles(graph(len(w),sum(v<<(3*i) for i,v in enumerate(w))))
 return any(not(a&b) for a,b in combinations(ts,2))
def subseq(u,w):
 k=0
 for x in w:
  if k<len(u) and x==u[k]:k+=1
 return k==len(u)
def main():
 bad=[]
 for m in (3,4):
  for w in product(range(8),repeat=m):
   if packing_bad(w) and not any(subseq(v,w) for v in bad):bad.append(w)
 need(Counter(map(len,bad))=={3:240,4:35},'minimal basis size')
 stored=json.load(open(P/'patterns.json'))
 need(bad==[tuple(x) for k in ('3','4') for x in stored['minimal_bad_patterns'][k]],'frozen basis')
 prefixes={w[:j] for w in bad for j in range(len(w))}
 need(len(prefixes)==84,'proper-prefix count including empty')
 maximal=[]
 for k in range(1,9):
  for A in combinations(range(8),k):
   if any(set(w)<=set(A) for w in bad):continue
   if not any(set(A)<set(B) for B in maximal):maximal.append(A)
 maximal=[A for A in maximal if not any(set(A)<set(B) for B in maximal)]
 need(maximal==[tuple(a) for a in stored['maximal_free_alphabets']],'maximal alphabets')
 expected=[]
 for A in combinations(range(8),2):
  x,y=A;diff=x^y
  if diff not in (1,2,4):continue
  a=diff.bit_length()-1;b=(a+1)%3;c=(a+2)%3
  # Core is a->b->c->a, hence b->c. Common c=1,b=0 gives forbidden bc-vertex triangle.
  if not(x>>b&1) and x>>c&1:continue
  expected.append(A)
 need(maximal==expected,'analytic binary-edge criterion')
 # Recognizer stores all currently realized proper prefixes as subsequences.
 def run(w):
  state={()};jumps=0
  for x in w:
   if any(v[:-1] in state and v[-1]==x for v in bad):return False,jumps
   fresh=state|{v for v in prefixes if v and v[-1]==x and v[:-1] in state}
   jumps+=fresh!=state;state=fresh
  need(jumps<=83,'strict transition count')
  return True,jumps
 max_jumps=0;accepted=0
 for m in range(6):
  for w in product(range(8),repeat=m):
   ok,jumps=run(w);need(ok==not_bad(w),'automaton versus actual disjoint triangles')
   accepted+=ok;max_jumps=max(max_jumps,jumps)
 # Deliberately drop a minimal forbidden pattern: that word becomes falsely accepted.
 hostile=bad[0];need(not any(subseq(v,hostile) for v in bad[1:]),'basis-deletion negative control')
 return dict(status='PASS',checks=checks,minimal_basis_by_length=dict(Counter(map(len,bad))),
  proper_prefixes=len(prefixes),maximal_free_alphabets=maximal,
  exhaustive_words_through_length5=sum(8**m for m in range(6)),accepted=accepted,
  maximum_observed_strict_updates=max_jumps,
  global_proved_limits=dict(strict_updates=83,free_phases=84,alphabet_size=2),
  boundary='Enumeration validates the finite basis; all-order equivalence follows because two disjoint triangles use at most four spine vertices. Phase count uses monotonicity, not a finite-length extrapolation.')
def not_bad(w):return not packing_bad(w)
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
