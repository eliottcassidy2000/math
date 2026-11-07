#!/usr/bin/env python3
"""Generate every tournament isomorphism class through order 8, from scratch.
No external package or repository data is used. Bit j in rows[i] means i->j.
"""
from itertools import permutations, product
from pathlib import Path
import json, hashlib
HERE=Path(__file__).resolve().parent
COUNTS=[1,1,1,2,4,12,56,456,6880]
def check(test,message):
 if not test: raise RuntimeError(message)
def canonical(rows):
 n=len(rows); colors=[r.bit_count() for r in rows]
 while True:
  signatures=[(colors[u],tuple(sorted(colors[v] for v in range(n) if rows[u]>>v&1))) for u in range(n)]
  palette={s:i for i,s in enumerate(sorted(set(signatures)))}
  new=[palette[s] for s in signatures]
  if len(set(new))==len(set(colors)): break
  colors=new
 cells=[[u for u in range(n) if new[u]==c] for c in sorted(set(new))]
 best=None
 for blocks in product(*(permutations(c) for c in cells)):
  order=sum(blocks,()); inv=[0]*n
  for j,u in enumerate(order): inv[u]=j
  key=tuple(sum(1<<inv[v] for v in range(n) if rows[u]>>v&1) for u in order)
  if best is None or key<best: best=key
 return best

def main():
 classes={(0,)}; receipts=[]
 for n in range(1,9):
  if n>1:
   m=n-1; nxt=set()
   for rows in classes:
    for mask in range(1<<m):
     extended=tuple(r|(((mask>>i)&1)<<m) for i,r in enumerate(rows))+(((1<<m)-1)^mask,)
     nxt.add(canonical(extended))
   classes=nxt
  check(len(classes)==COUNTS[n],f'wrong class count at {n}')
  data=''.join(' '.join(map(str,rows))+'\n' for rows in sorted(classes))
  (HERE/f'classes{n}.rows').write_text(data)
  record={'n':n,'classes':len(classes),'sha256':hashlib.sha256(data.encode()).hexdigest()}
  receipts.append(record); print(json.dumps(record),flush=True)
 (HERE/'input_receipt.json').write_text(json.dumps(receipts,indent=2)+'\n')
if __name__=='__main__': main()
