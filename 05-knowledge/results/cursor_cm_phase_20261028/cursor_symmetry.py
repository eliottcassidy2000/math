"""Complete cursor census via prime-order endpoint symmetries, avoiding rigid inputs."""
from itertools import permutations
from math import lcm
from pathlib import Path
import json
import cursor_probe as c
HERE=Path(__file__).resolve().parent
def invariant_states(n,g):
 pairs=[(i,j) for i in range(n) for j in range(i+1,n)];index={e:k for k,e in enumerate(pairs)}
 seen=set();fixed=0;free=[]
 for k,(i,j) in enumerate(pairs):
  if k in seen:continue
  edges=[];a,b=i,j
  while True:
   u,v=sorted((a,b));ix=index[u,v];edges.append((ix,a<b,u,v));seen.add(ix)
   a,b=g[a],g[b]
   if (a,b)==(i,j):break
  signs={forward for ix,forward,u,v in edges if v==u+1}
  if len(signs)>1:return
  plus=sum(1<<ix for ix,forward,u,v in edges if forward)
  minus=sum(1<<ix for ix,forward,u,v in edges if not forward)
  if signs:fixed|=plus if next(iter(signs)) else minus
  else:free.append((plus,minus))
 masks=[fixed]
 for a,b in free:masks=[v|x for v in masks for x in (a,b)]
 for mask in masks:
  rows=[0]*n
  for k,(i,j) in enumerate(pairs):
   a,b=(i,j) if mask>>k&1 else (j,i);rows[a]|=1<<b
  yield tuple(rows)
def order(g):
 result=1;seen=set()
 for i in range(len(g)):
  if i in seen:continue
  j=i;s=0
  while j not in seen:seen.add(j);s+=1;j=g[j]
  result=lcm(result,s)
 return result
def main():
 output=[]
 for n in (7,8,9):
  bank=set();symmetries=0
  for middle in permutations(range(1,n-1)):
   g=(0,)+middle+(n-1,)
   if order(g) not in (3,5,7):continue
   symmetries+=1
   bank.update(invariant_states(n,g))
  print(json.dumps({'n':n,'prime_order_permutations':symmetries,'nontrivial_endpoint_states':len(bank)}),flush=True)
  dual=[];collisions=[]
  for t in sorted(bank):
   a=c.analyze(t)
   if a is not None:
    dual.append(a)
    if a['collision']:collisions.append(a)
  row={'n':n,'prime_order_permutations':symmetries,'nontrivial_endpoint_states':len(bank),
       'both_endpoint_groups_nontrivial':len(dual),'candidate_paths':sum(a['candidate_paths'] for a in dual),
       'collisions':collisions,'dual_examples':dual[:8]}
  output.append(row)
  print(json.dumps({k:v for k,v in row.items() if k not in ('collisions','dual_examples')}|{'collisions':len(collisions)}),flush=True)
 result={'status':'FINITE-EXACT complete via Cauchy prime-order cover','rows':output}
 (HERE/'cursor_symmetry.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
