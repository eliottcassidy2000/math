#!/usr/bin/env python3
"""Independent Hamilton counting by actual feedback ports and a transitive chain."""
from collections import defaultdict
from itertools import combinations,permutations
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261030'))
from audit_feedback import subset_dp,require
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031b'))
from audit_inverse_witnesses import adjacency

def chain_order(rows,ports):
 rest=[v for v in range(len(rows)) if v not in ports]
 order=sorted(rest,key=lambda v:-sum(rows[v]>>u&1 for u in rest))
 return order if all(rows[u]>>v&1 for i,u in enumerate(order) for v in order[i+1:]) else None

def count(rows,ports):
 chain=chain_order(rows,ports);require(chain is not None,'declared ports leave an actual transitive chain')
 k=len(ports);total=0
 for order in permutations(ports):
  states={(0,)*(k+1):1}
  for v in chain:
   following=defaultdict(int)
   for state,multiplicity in states.items():
    for gap in range(k+1):
     if state[gap]==0 and gap>0 and not(rows[order[gap-1]]>>v&1):continue
     updated=list(state);updated[gap]=1 if gap==k or rows[v]>>order[gap]&1 else 2
     following[tuple(updated)]+=multiplicity
   states=following
  for state,multiplicity in states.items():
   if 2 in state:continue
   if any(state[j]==0 and not(rows[order[j-1]]>>order[j]&1) for j in range(1,k)):continue
   total+=multiplicity
 return total

def word_ports(word):
 ports=[1];n=2
 for kind,cut in word:
  if kind=='D':ports.append(n);n+=4
  else:n+=1
 return ports

def main():
 small=0
 for n in range(1,6):
  edges=list(combinations(range(n),2))
  for code in range(1<<len(edges)):
   rows=[0]*n
   for bit,(u,v) in enumerate(edges):
    if code>>bit&1:rows[u]|=1<<v
    else:rows[v]|=1<<u
   found=False
   for k in range(min(3,n)+1):
    for ports in combinations(range(n),k):
     if chain_order(rows,ports) is not None:
      require(count(rows,ports)==subset_dp(rows),'small literal subset DP');small+=1;found=True;break
    if found:break
   require(found,'small feedback coverage')
 receipts=[]
 for file in ('guarded_targets.json','formula_witnesses.json'):
  for r in json.loads((HERE/file).read_text())['results']:
   if r.get('word') is None:continue
   rows=adjacency(r['word']);ports=word_ports(r['word']);actual=count(rows,ports)
   require(actual==r['target'],'large actual-graph port count')
   receipts.append({'target':r['target'],'vertices':len(rows),'ports':ports,'actual':actual})
 # Actual long words, independently checked against the explicit arithmetic form.
 from one_block_formula import formula
 from find_formula_witness import seed_word,suffix_word,block_word
 long_checks=0;largest_order=0;largest_count=0
 for k in [-1,0,1,2,3,4,5]:
  name='D' if k==-1 else 'DD' if k==0 else f'DT{k}D'
  records=[list(map(int,line.split())) for line in(HERE/f'{name}_profiles.tsv').read_text().splitlines()]
  for j in range(16):
   H,E,p,q,s,t,F,mask=records[(j*149+37)%len(records)]
   b=dict(zip(('H','E','p','q','s','t','F'),(H,E,p,q,s,t,F)))
   x=(1<<(20+j%7))+j*137;y=(1<<(19+(3*j)%9))-1-j*97
   word=seed_word(x)+block_word(name,mask)+suffix_word(y);rows=adjacency(word)
   actual=count(rows,word_ports(word));expected=formula(b,x,y)
   require(actual==expected,'long literal feedback-port count')
   long_checks+=1;largest_order=max(largest_order,len(rows));largest_count=max(largest_count,actual)
 print(json.dumps({'status':'passed','small_complete_parent_checks':small,'positive_targets':receipts,
                   'long_word_checks':long_checks,'largest_order':largest_order,'largest_count':largest_count,
                   'scope':'Counts actual adjacency using port permutations and chain-gap partitions, with no block-transfer or digital formula in the counting algorithm.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
