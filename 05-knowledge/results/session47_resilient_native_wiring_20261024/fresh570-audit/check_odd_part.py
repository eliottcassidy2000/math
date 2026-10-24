#!/usr/bin/env python3
"""Exact witness that the full original critical group has nontrivial odd part.
No determinant calculation or prime factorization is needed: a literal weighted
spanning tree is heavier than the independently certified entire 2-primary order.
"""
from pathlib import Path
from collections import Counter
import argparse,csv,gzip,hashlib,json,math
H=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();checks=0
def need(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)
path=H/'inputs/original/metagraph_order8_20261023/merged_edges.tsv.gz'
need(hashlib.sha256(path.read_bytes()).hexdigest()=='67107747944898a23248c2a50fc10f81bac4c217e9aedddac4c8116809ee4f8a','original edge hash')
edges=[]
for r in csv.DictReader(gzip.open(path,'rt'),delimiter='\t'):
 i,j,b,k=(int(r[c]) for c in ['source','target','blue','black']);w=b+k
 need(0<=i<=j<3528 and b>=0 and k>=0,'original edge schema')
 if i!=j and w:edges.append((w,i,j))
parent=list(range(3528))
def find(i):
 while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
 return i
selected=[];product=1
for w,i,j in sorted(edges,key=lambda r:(-r[0],r[1],r[2])):
 u,v=find(i),find(j)
 if u!=v:parent[u]=v;selected.append([i,j,w]);product*=w
 if len(selected)==3527:break
need(len(selected)==3527 and len({find(i) for i in range(3528)})==1,'literal selected weighted spanning tree')
# A second reconstruction authenticates acyclicity and each original conductance.
edge={(i,j):w for w,i,j in edges};adj=[[] for _ in range(3528)]
for i,j,w in selected:
 need(edge[i,j]==w and w>0,'certificate edge is an original positive conductance');adj[i].append(j);adj[j].append(i)
seen={0};todo=[(0,-1)]
while todo:
 i,prev=todo.pop()
 for j in adj[i]:
  if j==prev:continue
  need(j not in seen,'independent DFS acyclicity');seen.add(j);todo.append((j,i))
need(len(seen)==3528,'independent DFS spanning')
fil=json.loads((H/'guard-completeness.normal.json').read_text())['filtration'];ranks=fil['valuation_ranks'];v2=sum(i*r for i,r in enumerate(ranks,1))
need(ranks==[3358,4,0,0,0,0,0,1] and v2==3374,'independently replayed original 2-primary order')
need(product==math.prod(w for i,j,w in selected) and product.bit_length()==9205 and product>2**v2,'tree weight exceeds full 2-primary order')
r={'status':'PASS: original critical group has a nontrivial odd primary part, so its FULL exponent is not256','checks':checks,'vertices':3528,'selected_tree_edges':selected,'tree_conductance_product_decimal':str(product),'tree_product_bit_length':product.bit_length(),'tree_weight_multiplicities':dict(sorted(Counter(w for i,j,w in selected).items())),'critical_two_primary_order_log2':v2,'odd_primary_order_lower_bound_log2':product.bit_length()-1-v2,'proof':'Weighted matrix-tree theorem gives |K|=tau>=product>2^3374. Complete original dyadic filtration gives v2(|K|)=3374. Since product>=2^9204, |K|/2^3374>=2^5830; this integer is odd and therefore is actually greater than2^5830. Any finite abelian group with such an order has nontrivial odd primary part; its exponent cannot be256.','scope':'A lower-bound audit witness, not an exact tree count or an odd-factor classification.'}
a.out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'tree_product_bit_length':product.bit_length(),'v2_critical_order':v2}))
