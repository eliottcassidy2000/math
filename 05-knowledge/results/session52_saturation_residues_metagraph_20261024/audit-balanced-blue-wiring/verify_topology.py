#!/usr/bin/env python3
"""Replay a small spanning two-connected subgraph against literal adjacency."""
import json,sys
from pathlib import Path
P=Path(__file__).resolve().parent
source=Path(sys.argv[1]) if len(sys.argv)>1 else P/'producer'/'normal.json'
c=json.loads((P/'topology-certificate.json').read_text());d=json.loads(source.read_text())
def need(x,msg):
 if not x:raise RuntimeError(msg)
edges={tuple(sorted((a,b))) for a,b,m in d['adjacency_multiplicities'] if m>0 and a!=b}
verts={r['atlas_index'] for r in d['nodes']};t=c['triangle'];need(len(t)==3 and len(set(t))==3,'triangle distinct');need(all(tuple(sorted((t[i],t[(i+1)%3]))) in edges for i in range(3)),'triangle edges');used=set(t)
for a,v,b in c['two_edge_ears']:
 need(a in used and b in used and a!=b and v not in used,'ear joins two distinct old nodes through one new node');need(tuple(sorted((a,v))) in edges and tuple(sorted((v,b))) in edges,'ear edges');used.add(v)
need(used==verts and len(used)==229,'spanning all229');need(c['certificate_edges']==455,'455 edges')
print(json.dumps({'status':'PASS','vertices':229,'edges':455,'conclusion':'spanning two-connected simple subgraph certifies blue induced graph connected with no articulations or bridges'}))
