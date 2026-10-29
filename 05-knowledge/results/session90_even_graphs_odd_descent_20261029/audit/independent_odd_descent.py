#!/usr/bin/env python3
"""Nonabelian odd coprime-action audit on the nine-vertex wreath tournament."""
from itertools import product,permutations
from independent_holonomy import mul,inv,power
import json
from collections import Counter
C=0
def ck(x,msg):
 global C
 C+=1
 if not x:raise RuntimeError(msg)
def move(p,x):return sum(((x>>i)&1)<<p[i] for i in range(len(p)))
def orbits(points,G,action):
 todo=set(points);out=[]
 while todo:
  x=min(todo);orb={action(g,x) for g in G};ck(orb<=todo,'orbit partition');todo-=orb;out.append(orb)
 return out

def main():
 A=[tuple(3*((i//3+a)%3)+(i%3+b[i//3])%3 for i in range(9)) for a in range(3) for b in product(range(3),repeat=3)]
 tau=tuple(3*((-(i//3))%3)+(-(i%3))%3 for i in range(9));central=[a for a in A if mul(a,tau)==mul(tau,a)]
 allorb=orbits(range(512),A,move);fixed=[x for x in range(512) if move(tau,x)==x]
 stable=[o for o in allorb if {move(tau,x) for x in o}==o];fo=orbits(fixed,central,move)
 ck(len(stable)==len(fo),'nonabelian fixed-orbit descent')
 repairs=0
 for x in fixed:
  for g in A:
   y=move(g,x)
   if y not in fixed:continue
   ag=mul(mul(tau,g),tau);h=mul(inv(g),ag);a=power(h,41);c=mul(g,a)
   ck(move(a,x)==x,'square root retains joint subset')
   ck(c in central and move(c,x)==y,'corrected transporter');repairs+=1
 def arc(i,j):
  a,b=divmod(i,3);c,d=divmod(j,3)
  return (c-a)%3==1 if a!=c else (d-b)%3==1
 H=0;fixedpaths=[]
 for p in permutations(range(9)):
  if all(arc(x,y) for x,y in zip(p,p[1:])):
   H+=1
   if tuple(tau[x] for x in p[::-1])==p:fixedpaths.append(p)
 pairs=[(i,tau[i]) for i in range(9) if i<tau[i]];fold=0
 for seq in permutations(range(4)):
  for signs in product(range(2),repeat=4):
   p=tuple(pairs[seq[i]][signs[i]] for i in range(4))
   fold+=all(arc(x,y) for x,y in zip(p,p[1:])) and arc(p[-1],0)
 ck(fold==len(fixedpaths),'independent half-word observer')
 blue=orbits(fixedpaths,central,lambda g,p:tuple(g[v] for v in p))
 ck(fold==len(blue)*len(central),'free fixed path orbit')
 pi={v:i for i,pair in enumerate(pairs) for v in pair};vertices=sorted(pi)
 states=[(mask,v) for mask in range(1,16) for v in vertices if mask>>pi[v]&1]
 def act(g,state):
  mask,v=state
  return sum(1<<pi[g[pairs[i][0]]] for i in range(4) if mask>>i&1),g[v]
 so=orbits(states,central,act);owner={s:i for i,o in enumerate(so) for s in o}
 def succ(s):
  mask,v=s
  return [(mask|(1<<pi[w]),w) for w in vertices if not(mask>>pi[w]&1) and arc(v,w)]
 transitions=[]
 for o in so:
  counts=Counter(owner[y] for y in succ(next(iter(o))));transitions.append(counts)
  for state in o:ck(Counter(owner[y] for y in succ(state))==counts,'quotient equitability')
 totals=[len(o) if next(iter(o))[0].bit_count()==1 else 0 for o in so]
 wrong=[1 if next(iter(o))[0].bit_count()==1 else 0 for o in so]
 for i in sorted(range(len(so)),key=lambda i:next(iter(so[i]))[0].bit_count()):
  for j,w in transitions[i].items():totals[j]+=totals[i]*w;wrong[j]+=wrong[i]*w
 # Independent literal partial words, not a second dynamic program.
 literal=Counter()
 for k in range(1,5):
  for seq in permutations(range(4),k):
   for signs in product(range(2),repeat=k):
    word=[pairs[seq[j]][signs[j]] for j in range(k)]
    if all(arc(x,y) for x,y in zip(word,word[1:])):literal[sum(1<<i for i in seq),word[-1]]+=1
 for i,o in enumerate(so):ck(totals[i]==sum(literal[s] for s in o),'every orbit total versus partial words')
 terminal=[i for i,o in enumerate(so) if next(iter(o))[0]==15 and arc(next(iter(o))[1],0)]
 ck(sum(totals[i] for i in terminal)==fold,'weighted quotient final')
 ck(sum(wrong[i] for i in terminal)==19,'wrong unit orbit initialization hostile')
 quotient=dict(states=len(states),orbits=len(so),arcs=sum(len(succ(s)) for s in states),weighted_entries=sum(len(x) for x in transitions),orbit_sizes=dict(Counter(map(len,so))),wrong_unit_result=sum(wrong[i] for i in terminal))
 print(json.dumps(dict(status='PASS',checks=C,ordinary_subset_orbits=len(allorb),fixed_subsets=len(fixed),stable_subset_orbits=len(stable),fixed_subset_orbits=len(fo),repaired_transporters=repairs,H=H,Aut=len(A),centralizer=len(central),fixed_paths=fold,blue_path_orbits=len(blue),quotient=quotient),indent=2))
if __name__=='__main__':main()
