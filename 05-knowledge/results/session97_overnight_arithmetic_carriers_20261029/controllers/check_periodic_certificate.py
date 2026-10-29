#!/usr/bin/env python3
from fractions import Fraction as F
from itertools import product
from check_integer_termination import require
from check_exceptional_cylinders import WORDS,chart
import json

def main():
 charts=[chart(w) for w in WORDS]
 alpha=[-F(c['C'],3**c['J']-(1<<c['D'])) for c in charts]
 low=min(alpha);high=max(alpha)
 require(low==-F(16177,139) and high==-1,'universal periodic interval')
 edges={}
 for n in range(-115,0,2):
  matches=[c for c in charts if n%c['modulus']==c['residue']]
  require(len(matches)<=1,'disjoint')
  if matches:
   c=matches[0];edges[n]=(3**c['J']*n+c['C'])//(1<<c['D'])
   require(c['word']==(1,) or n==-3,'only two chart types')
 require(edges[-3]==-1 and edges[-1]==-1,'terminal cycle')
 require(all(y<n for n,y in edges.items() if n not in (-3,-1)),'no other cycle')
 # Independent exact fixedpoints for all failurewords of block length≤5.
 integer_words=[];tested=0
 for length in range(1,6):
  for indices in product(range(5),repeat=length):
   A=1;C=D=0
   for i in indices:
    c=charts[i];aa=3**c['J'];C=aa*C+(1<<D)*c['C'];A*=aa;D+=c['D']
   fixed=-F(C,A-(1<<D));tested+=1
   require(low<=fixed<=high,'fixedpoint convex hull')
   if fixed.denominator==1:
    require(fixed==-1 and set(indices)=={0},'integer periodic classification')
    integer_words.append(indices)
 preimages={}
 for y in [-1,-3]:
  preimages[str(y)]=[]
  for i,c in enumerate(charts):
   x=F((1<<c['D'])*y-c['C'],3**c['J'])
   if x.denominator==1:preimages[str(y)].append([i,int(x)])
 require(preimages=={'-1':[[0,-1],[1,-3]],'-3':[]},'integer eventual periodic preimages')
 print(json.dumps({'fixedpoints':[str(a) for a in alpha],'odd_integer_interval':[-115,-1],
  'periodic_block_words_checked':tested,'only_integer_cycle':[-1],
  'eventually_periodic_integer_sources':[-3,-1],'integer_preimages':preimages},sort_keys=True))
if __name__=='__main__':main()
