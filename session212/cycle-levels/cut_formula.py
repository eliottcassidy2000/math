#!/usr/bin/env python3
"""Independent endpoint-cut formula checks on every class/cut through m=6."""
from pathlib import Path
from itertools import permutations
import json
from verify import check,tournament,rooted_paths
HERE=Path(__file__).resolve().parent
checked=0
for m in range(2,7):
 for line in (HERE/f'classes{m}.rows').read_text().splitlines():
  rows=list(map(int,line.split()));P=[[0]*m for _ in range(m)]
  # This counter uses literal permutations, independent of both subset DPs.
  for p in permutations(range(m)):
   if all(rows[x]>>y&1 for x,y in zip(p,p[1:])):P[p[0]][p[-1]]+=1
  for A in range(1<<m):
   extension=[r|((not(A>>i&1))<<m) for i,r in enumerate(rows)]+[A]
   predicted=sum(P[a][b] for a in range(m) for b in range(m) if A>>a&1 and not(A>>b&1))
   tournament(extension)
   check(rooted_paths(extension,m,True)==predicted,f'endpoint-cut formula m={m}, A={A}')
   checked+=1
print(json.dumps({'literal_endpoint_matrix_and_cut_checks':checked,'parent_orders':'2..6','extension_orders':'3..7'},indent=2))
