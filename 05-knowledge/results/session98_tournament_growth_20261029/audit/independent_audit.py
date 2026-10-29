#!/usr/bin/env python3
"""Independent labelled-action tests. No imports from deliverable checkers."""
import itertools as it, json
from fractions import Fraction as F
from collections import Counter
from math import factorial

def need(v,s):
 if not v: raise RuntimeError(s)
def edge_data(n):return list(it.combinations(range(n),2))
def act(s,p,edges):
 lookup={e:i for i,e in enumerate(edges)};o=[0]*len(s)
 for (i,j),v in zip(edges,s):
  a,b=p[i],p[j]
  if a>b:a,b=b,a;v=-v
  o[lookup[a,b]]=v
 return tuple(o)
def graph_orbits(n,alphabet):
 es=edge_data(n);ps=list(it.permutations(range(n)));seen=set();hist=Counter();primes=[]
 for s in it.product(alphabet,repeat=len(es)):
  if s in seen:continue
  images=[act(s,p,es) for p in ps];seen.update(images);hist[s.count(0)]+=1
  aut=[p for p,z in zip(ps,images) if z==s]
  if 0 in alphabet:
   v2=(len(aut)&-len(aut)).bit_length()-1
   need(v2<=s.count(0),'2-primary budget')
  else:
   lookup={e:v for e,v in zip(es,s)}
   def arrow(i,j):return lookup[i,j] if i<j else -lookup[j,i]
   prime=n>=3 and not any(all(len({arrow(z,u) for u in M})==1 for z in range(n) if z not in M) for k in range(2,n) for M in it.combinations(range(n),k))
   if prime:
    vertexorbits={min(p[i] for p in aut) for i in range(n)}
    primes.append({'state':s,'group_order':len(aut),'surplus_one':len(vertexorbits)})
 return hist,primes

def literal_burnside(n):
 es=edge_data(n);total=0
 for p in it.permutations(range(n)):
  # Walk directed ordered pairs; orientation consistency around unordered orbits.
  unseen=set(es);free=0;bad=False
  while unseen:
   i,j=min(unseen);a,b=i,j
   while True:
    unseen.discard(tuple(sorted((a,b))));a,b=p[a],p[b]
    if {a,b}=={i,j}:break
   if (a,b)!=(i,j):bad=True;break
   free+=1
  if not bad:total+=2**free
 need(total%factorial(n)==0,'integral Burnside')
 return total//factorial(n)
def marked_counts(n):
 es=edge_data(n);ps=list(it.permutations(range(n)));seen=set();e=v=0;hostile=None
 for s in it.product((-1,1),repeat=len(es)):
  if s in seen:continue
  images=[act(s,p,es) for p in ps];seen.update(images)
  aut=[p for p,z in zip(ps,images) if z==s]
  v+=len({min(p[i] for p in aut) for i in range(n)})
  e+=len({min(tuple(sorted((p[i],p[j]))) for p in aut) for i,j in es})
  for k,(i,j) in enumerate(es):
   flip=s[:k]+(-s[k],)+s[k+1:]
   if flip not in images:continue
   def arrow(a,b):return s[es.index((a,b))] if a<b else -s[es.index((b,a))]
   twins=all(arrow(i,z)==arrow(j,z) for z in range(n) if z not in (i,j))
   if not twins:hostile={'n':n,'state':s,'flipped_pair':(i,j)}
 return e,v,hostile

def main():
 counts=[literal_burnside(n) for n in range(8)]
 need(counts==[1,1,1,2,4,12,56,456],'literal all-permutation control')
 partial={n:dict(graph_orbits(n,(-1,0,1))[0]) for n in range(5)}
 marked={n:marked_counts(n) for n in range(1,6)}
 need(marked[4][2] is not None,'unmarked non-twin loop witness')
 for n in range(2,5):need(2*partial[n][1]==marked[n][0]+marked[n-1][1],'one-hole marked identity')
 primes=graph_orbits(5,(-1,1))[1]
 need(len(primes)==3,'prime five census')
 need(sorted(p['surplus_one'] for p in primes)==[1,5,5],'same scalar prime count hostile')
 # C3: a marked + at vertex 0. Stabilizer trivial. Forgetting marks merges two completions.
 completions=[(1,a,b) for a,b in it.product((-1,1),repeat=2)]
 orbit=lambda x:min(x[i:]+x[:i] for i in range(3))
 need(len(completions)==4 and len({orbit(x) for x in completions})==3,'flag hostile')
 # One missing pair: transposition fixes it although it fixes no tournament.
 need(partial[2]=={0:1,1:1},'one-defect hostile')
 # Baseline-one assumption is necessary for a C3 palette with two tagged singletons.
 baseline=(2**3+2*2)//3
 need(baseline==4 and baseline!=1,'multi-singleton hostile')
 print(json.dumps({'status':'PASS','literal_burnside':counts,'partial_graph_histograms':partial,'marked_pair_vertex_counts':{n:r[:2] for n,r in marked.items()},'unmarked_non_twin_loop':next((r[2] for r in marked.values() if r[2]),None),'five_vertex_primes':primes,'hostiles':{'equal_prime_count_different_surplus_one':[1,5],'marked_vs_unmarked_completions':[4,3],'two_tagged_singletons_C3_orbits':baseline}},indent=2))
if __name__=='__main__':main()
