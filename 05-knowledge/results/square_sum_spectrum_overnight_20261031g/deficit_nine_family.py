#!/usr/bin/env python3
"""Actual adjacency witnesses for the complete positive deficit-nine line."""
import json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031d'))
from bank_query import query
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031c'))
from audit_port_words import adjacency
from audit_power_sources import path_counts,marked,append_zero,need
PARENT=[12,121,2,100,45,69,21]
MASK=67

def join(left,right):
 n=len(left);m=len(right)
 return [r|(((1<<m)-1)<<n) for r in left]+[r<<n for r in right]

def canonical(n):
 if n<9:raise ValueError('six-arc family begins at order9')
 rows=[((1<<n)-1)^((1<<(v+1))-1) for v in range(n)]
 for u,v in [(6,0),(4,1),(5,2),(n-1,0),(n-1,4),(n-1,5)]:
  rows[v]^=1<<u;rows[u]^=1<<v
 return rows

def construct(k):
 if k<5:raise ValueError('positive admissible deficit-nine exponents start at5; k4 gives the forbidden count7')
 if k>=10:
  rows,x=marked(PARENT,MASK)
  for _ in range(k-9):rows=append_zero(rows,x)
  return rows
 if k==8:return join(adjacency(query(13)['word']),adjacency(query(19)['word']))
 return adjacency(query((1<<k)-9)['word'])

def main():
 records=[]
 for k in range(5,18):
  rows=construct(k);H,_=path_counts(rows);need(H==(1<<k)-9,'all actual deficit-nine witnesses')
  if k>=10:
   need(len(rows)==k-1,'one vertex per exponent, exact clock')
   p=[1,4,3,5,6,0,2]+list(range(8,len(rows)))+[7]
   relabelled=[sum(((rows[u]>>v)&1)<<j for j,v in enumerate(p)) for u in p]
   need(relabelled==canonical(k-1),'explicit six reversed arcs, with the marked root re-labelled last')
   need(path_counts(canonical(k-1))[0]==H,'canonical direct graph count')
  records.append({'k':k,'raw_H':H,'vertices':len(rows),'rows':rows})
 print(json.dumps({'status':'passed','independent_counted_exponents':[5,17],'all_size_statement':'For every integer k>=5 there is a tournament with raw H=2^k-9. For k>=10 the explicit construction has k-1 vertices. The only positive missing exponent is k4, whose H is7.','canonical_one_based_backward_arcs':['7->1','5->2','6->3','n->1','n->5','n->6'],'witnesses':records,'scope':'The all-size conclusion follows the native zero-tail recurrence; these adjacency counts check its actual source and finite boundary.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
