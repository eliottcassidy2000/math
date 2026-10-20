#!/usr/bin/env python3
"""Exhaust all subdirect F3 spaces and their actual ambient conjugacies."""
from itertools import product,combinations
from collections import Counter
from pathlib import Path
import json
import audit_census as a

def span(vectors,d):
 return frozenset(tuple(sum(c*v[i] for c,v in zip(coef,vectors))%3 for i in range(d)) for coef in product(range(3),repeat=len(vectors)))
def spaces(d):
 vec=list(product(range(3),repeat=d));out=set()
 for r in range(1,d+1):
  for vs in combinations(vec[1:],r):out.add(span(vs,d))
 return {V for V in out if all({v[i] for v in V}=={0,1,2} for i in range(d))}
def orbit(V,frob):
 out=set()
 for eps in product((1,2),repeat=2):
  for swap in (False,True):
   W=[]
   for v in V:
    z=list(v[-2:]);z=[x*y%3 for x,y in zip(z,eps)]
    if swap:z=z[::-1]
    W.append(tuple(([v[0]] if frob else [])+z))
   out.add(frozenset(W))
 return out

def main():
 receipts=[];allcanon=[]
 for d in (2,3):
  left=spaces(d);count=len(left);representatives=[]
  while left:
   V=min(left,key=lambda v:tuple(sorted(v)));O=orbit(V,d==3)
   a.m.need(O<=spaces(d),'conjugacy stays subdirect');left-=O
   representatives.append({'order_quotient':len(V),'vectors':sorted(V),'conjugacy_orbit_size':len(O)})
  a.m.need(len(representatives)==(2 if d==2 else 5),'subdirect orbit count')
  allcanon.extend((d,min(tuple(sorted(W)) for W in orbit(frozenset(map(tuple,r['vectors'])),d==3))) for r in representatives)
  receipts.append({'dimension':d,'raw_subdirect_spaces':count,'ambient_classes':len(representatives),'representatives':representatives})
 actual=[]
 for n,name,gens in a.allcases():
  if n!=13 or '(7,3,3)' not in name:continue
  H=a.m.closure(gens,n);vectors=set();frob=False
  for h in H:
   mult=(h[1]-h[0])%7;u={1:0,2:1,4:2}[mult];frob|=u!=0
   vectors.add((u,h[7]-7,h[10]-10))
  d=3 if frob else 2;V=frozenset(v if frob else v[1:] for v in vectors)
  canonical=min(tuple(sorted(W)) for W in orbit(V,frob));actual.append((d,canonical))
 a.m.need(len(actual)==7 and len(set(actual))==7,'seven proposed rows distinct')
 a.m.need(set(actual)==set(allcanon),'seven proposed rows exhaust all conjugacy classes')
 Path(__file__).with_suffix('.json').write_text(json.dumps({'status':'PASS','checks':a.m.CHECKS,'families':receipts},indent=2)+'\n');print([(r['dimension'],r['raw_subdirect_spaces'],r['ambient_classes']) for r in receipts])
if __name__=='__main__':main()
