#!/usr/bin/env python3
"""Direct arc-based audit, independent of C++ bit-mask move templates."""
from pathlib import Path
import json
from verify_realizability import canonical_certificate, burnside
P=Path(__file__).resolve().parent

def check(n):
 f=[(i,j) for i in range(n) for j in range(i+2,n)];ix={p:k for k,p in enumerate(f)};N=1<<len(f)
 parent=list(range(N));size=[1]*N
 def root(a):
  while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
  return a
 def join(a,b):
  a,b=root(a),root(b)
  if a!=b:
   if size[a]<size[b]:a,b=b,a
   parent[b]=a;size[a]+=size[b]
 valid=0
 for x in range(N):
  out=[[0]*n for _ in range(n)]
  for i in range(n):
   for j in range(i+1,n):
    bit=1 if j==i+1 else x>>ix[i,j]&1;out[i][j]=bit;out[j][i]=1-bit
  for a in range(n):
   for b in range(a+2,n):
    word=list(range(n));v=word.pop(a);word.insert(b,v)
    if all(out[word[i]][word[i+1]] for i in range(n-1)):
     valid+=1;y=sum(out[word[i]][word[j]]<<k for k,(i,j) in enumerate(f));join(x,y)
 assert valid*16==N*(n*n+n-4)
 ordinary=len({root(x) for x in range(N)});assert ordinary==burnside(n)['ordinary']
 for x in range(N):
  rx=sum((x>>k&1)<<ix[n-1-j,n-1-i] for k,(i,j) in enumerate(f));join(x,rx)
 labels=[root(x) for x in range(N)];cpp=list(map(int,(P/f'local_assignment_n{n}.txt').read_text().split()))
 canon=[h['node'] for h in canonical_certificate(n)['halfedges']]
 for other in (cpp,canon):
  forward={};reverse={}
  for a,b in zip(labels,other):
   assert forward.setdefault(a,b)==b and reverse.setdefault(b,a)==a
 merged=len(set(labels));assert merged==burnside(n)['merged']
 return {'n':n,'tilings':N,'ordinary_components':ordinary,'merged_components':merged,'valid_moves':valid,'direct_arcs_vs_cpp_partition':'PASS','direct_arcs_vs_permutation_canonical_partition':'PASS'}

if __name__=='__main__':
 reports=[check(n) for n in range(3,8)]
 cpp=[json.loads(s) for s in (P/'local_move_classifier.out').read_text().splitlines()]
 for r in cpp:
  n=r['n'];b=burnside(n)
  assert r['ordinary_components']==b['ordinary']
  assert r['self_converse_components']==b['self_converse']
  assert r['merged_components']==b['merged']
  assert r['ordinary_even_size_components']==0
  assert r['valid_directed_moves']*16==r['tilings']*(n*n+n-4)
  hist={int(s):c for s,c in r['merged_size_histogram'].items()}
  assert sum(hist.values())==b['merged'] and sum(s*c for s,c in hist.items())==r['tilings']
  assert sum(c for s,c in hist.items() if s%2)==b['self_converse']
  assert all(s%2 or s%4==2 for s in hist)
 print(json.dumps({'direct_controls':reports,'burnside_and_inventory_controls_through_n8':'PASS'},indent=2))
