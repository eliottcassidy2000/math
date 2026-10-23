#!/usr/bin/env python3
"""Only the true reflection-fixed cube: broader blue arithmetic probe."""
import itertools as it,json,time
from collections import Counter
from pathlib import Path
import sympy as sp
from check_dyadic_network import lap,need,rank2,components
from check_dyadic_lift import dyadic_smith
from check_bockstein import nullbasis2
HERE=Path(__file__).resolve().parent

def run(n):
 free=[(i,j) for i in range(n) for j in range(i+2,n)];pairs=list(it.combinations(range(n),2));orbs=[];seen=set()
 for e in free:
  if e in seen:continue
  orb=sorted(set([e,(n-1-e[1],n-1-e[0])]))
  seen.update(orb);orbs.append(orb)
 full=(1<<len(orbs))-1;codes=[]
 for x in range(full+1):
  a=[[0]*n for _ in range(n)]
  for i in range(n-1):a[i][i+1]=1
  for k,orb in enumerate(orbs):
   for i,j in orb:
    u,v=(i,j) if x>>k&1 else(j,i);a[u][v]=1
  d=[sum(row) for row in a];blocks=[[v for v in range(n) if d[v]==s] for s in sorted(set(d))];best=None
  for chunks in it.product(*(it.permutations(b) for b in blocks)):
   p=tuple(v for chunk in chunks for v in chunk)
   code=sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(pairs))
   if best is None or code<best:best=code
  codes.append((tuple(sorted(d)),best))
 keys=sorted(set(codes));ix={key:i for i,key in enumerate(keys)};nodes=[ix[c] for c in codes];r=len(keys);edges=Counter()
 for x in range(full+1):
  if x<(x^full):edges[tuple(sorted((nodes[x],nodes[x^full])))]+=1
 es=[(i,j,w) for(i,j),w in edges.items()];L=lap(r,es);cs=components(r,es);vals=[];forms=[]
 for c in cs:
  A=L.extract(c[:-1],c[:-1]);v=[e for e in dyadic_smith(A.tolist()) if e];vals.extend(v)
  bs=nullbasis2(A.tolist());V=sp.Matrix(len(bs),A.rows,lambda i,j:bs[i][j]);Q=(V*A*V.T).applyfunc(lambda z:int(z)//2%2)
  forms.append(dict(size=len(c),bicycle=len(bs),half_form_rank=rank2(Q.tolist()),alternating=not any(Q[i,i] for i in range(Q.rows))))
 out=dict(n=n,blue_presentations=full+1,selfconverse_nodes=r,components=len(cs),bicycle_dimension=r-rank2(L.tolist())-len(cs),two_primary_exponents=sorted(vals),half_forms=forms,edges=es)
 if len(cs)==1:
  basis=nullbasis2(L[:-1,:-1].tolist())
  masks=[sum(v<<i for i,v in enumerate(row)) for row in basis]
  if len(masks)<=12 and not forms[0]['alternating']:
   best=None
   for state in range(1<<len(masks)):
    mask=0
    for i,m in enumerate(masks):
     if state>>i&1:mask^=m
    count=sum(w for i,j,w in es if ((mask>>i)^(mask>>j))&1)
    if count%4==2 and (best is None or count<best[0]):best=count,mask
   need(best is not None);count,mask=best
   cut=[e for e in es if ((mask>>e[0])^(mask>>e[1]))&1]
   need(all(sum(w for i,j,w in cut if v in(i,j))%2==0 for v in range(r)))
   out['nonalternating_bicycle_witness']=dict(subset_nodes=[i for i in range(r) if mask>>i&1],crossing_blue_lines=count,cut_edge_multiplicities=cut,every_crossing_degree_even=True,half_energy_mod2=1,minimum_among_all_binary_harmonic_cuts=True,kernel_states_checked=1<<len(masks))
 (HERE/f'blue_only_n{n}.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k not in ('edges','nonalternating_bicycle_witness')},flush=True)

if __name__=='__main__':
 import sys
 for n in map(int,sys.argv[1:] or ['8']):run(n)
