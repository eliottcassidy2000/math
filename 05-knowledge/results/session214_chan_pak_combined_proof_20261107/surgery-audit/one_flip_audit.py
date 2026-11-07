#!/usr/bin/env python3
"""Literal directed-word verification of two physical one-arc escapes."""
import json
from pathlib import Path
import sys

def need(x,why):
 if not x:raise RuntimeError(why)

def count(rows):
 n=len(rows);full=(1<<n)-1
 def walk(last,used):
  if used==full:return 1
  total=0;todo=rows[last]&~used
  while todo:
   b=todo&-todo;todo-=b;total+=walk(b.bit_length()-1,used|b)
  return total
 return sum(walk(i,1<<i) for i in range(n))

source=Path(sys.argv[1]) if len(sys.argv)>1 else Path('../surgery/ONE_FLIP_RESULT.json')
rows=[]
for c in json.loads(source.read_text())['cases']:
 m,q,tau=c['m'],c['q'],c['tau'];L=m+2;a=L;b=L+1
 p=c['parent_rows'];y=c['child_rows'];n=L+2
 for z in (p,y):
  need(len(z)==n,'order')
  for i in range(n):
   need(z[i]>=0 and not(z[i]>>i&1) and z[i]<(1<<n),'row')
   for j in range(i):need(((z[i]>>j)&1)+((z[j]>>i)&1)==1,'edge')
 for i in range(L):
  for j in range(i+1,L):need(p[i]>>j&1,'chain')
 need((p[a]&((1<<L)-1))==2,'comb sole exit')
 need(p[b]&1,'comb x0')
 need(((p[1]>>b)&1)==tau,'comb tau')
 for i in range(m):need(((p[i+2]>>b)&1)==((q>>i)&1),'comb bits')
 need(((p[a]>>b)&1)==(c['apex_orientation']=='+'),'comb orientation')
 r=c['flip_chain_index']
 need([(i,p[i]^y[i]) for i in range(n) if p[i]!=y[i]]==[(r,1<<a),(a,1<<r)],'one flip')
 hp,hy=count(p),count(y)
 need(hp==c['parent_H'] and hy==c['child_H'],'literal counts')
 need((y[a]&((1<<L)-1)).bit_count()==2,'second exit')
 rows.append({'n':n,'parent_H':hp,'child_H':hy,'chain_vertex':r,'changed_arcs':1})
print(json.dumps({'status':'PASS','count_method':'literal directed Hamiltonian words','cases':rows},indent=2,sort_keys=True))
