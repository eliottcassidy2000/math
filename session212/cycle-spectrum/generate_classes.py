#!/usr/bin/env python3
"""Complete tournament representatives via exact isomorphism canonicalization.
Adapted with credit from Session197 lowcycle/verify.py's exact generator.
No nauty, network access, or precomputed representatives are required.
"""
from itertools import permutations, product
import argparse
from pathlib import Path

def canonical(rows):
 n=len(rows);colors=[r.bit_count() for r in rows]
 while True:
  signatures=[(colors[u],tuple(sorted(colors[v] for v in range(n) if rows[u]>>v&1))) for u in range(n)]
  palette={s:i for i,s in enumerate(sorted(set(signatures)))}
  new=[palette[s] for s in signatures]
  if len(set(new))==len(set(colors)):break
  colors=new
 cells=[[v for v in range(n) if new[v]==c] for c in sorted(set(new))]
 best=None
 for blocks in product(*(permutations(c) for c in cells)):
  order=sum(blocks,());inverse=[0]*n
  for j,v in enumerate(order):inverse[v]=j
  code=tuple(sum(1<<inverse[v] for v in range(n) if rows[u]>>v&1) for u in order)
  if best is None or code<best:best=code
 return best

def generate(max_order):
 ts={(0,)};expected=[1,1,2,4,12,56,456,6880]
 for n in range(2,max_order+1):
  nxt=set()
  for rows in ts:
   for mask in range(1<<(n-1)):
    extended=[r|(((mask>>u)&1)<<(n-1)) for u,r in enumerate(rows)]+[((1<<(n-1))-1)^mask]
    nxt.add(canonical(extended))
  ts=nxt
  if len(ts)!=expected[n-1]:raise RuntimeError(('class count',n,len(ts)))
 return sorted(ts)

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--order',type=int,default=8);parser.add_argument('--output',default='tournaments8.txt');args=parser.parse_args()
 if not 1<=args.order<=8:raise RuntimeError('order must be 1..8')
 rows=generate(args.order);Path(args.output).write_text(''.join(' '.join(map(str,r))+'\n' for r in rows))
 print(f'order {args.order}: {len(rows)} exact isomorphism classes')
