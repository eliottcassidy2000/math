#!/usr/bin/env python3
import json,math,collections
from pathlib import Path
from sympy import Matrix,QQ
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
def locate(*p):
 for q in p:
  if q.exists():return q
 raise RuntimeError('source missing')
def rank(rows):
 if not rows:return 0
 return DomainMatrix.from_Matrix(Matrix(rows)).convert_to(QQ).rank()
r=json.load(open(locate(ROOT/'session40-square-sum/receipt.json',ROOT/'square-sum/receipt.json')))
expected={23:(2,0),25:(9,0),26:(11,3),27:(20,5),28:(25,7),29:(16,2),30:(17,3),31:(34,9),32:(45,12)};out=[]
for row in r['rows']:
 n=row['n']
 if n not in expected:continue
 E=[(i,j) for i in range(1,n+1) for j in range(i+1,n+1) if math.isqrt(i+j)**2==i+j]
 vectors=[];groups=collections.defaultdict(list)
 for p in row['paths']:
  edges={tuple(sorted(e)) for e in zip(p,p[1:])};v=[int(e in edges) for e in E];vectors.append(v);groups[tuple(sorted((p[0],p[-1])))].append(v)
 D=[[x-y for x,y in zip(v,vectors[0])] for v in vectors[1:]];W=[[x-y for x,y in zip(v,g[0])] for g in groups.values() for v in g[1:]]
 R=(rank(D),rank(W))
 if R!=expected[n]:raise RuntimeError((n,R))
 extra=None
 if n==32:
  B=[[int(v in e) for e in E] for v in range(1,n+1)];F=[[int(e==f) for e in E] for f in [(5,31),(9,27)]]
  if (rank(B),rank(B+F))!=(32,34):raise RuntimeError('basis rank')
  if any(sum(a*b for a,b in zip(w,f)) for w in W for f in F):raise RuntimeError('not annihilator')
  extra={'degree_rank':32,'extended_kernel_rank':34,'kernel_dimension':len(E)-R[1]}
 out.append({'n':n,'total_rank':R[0],'within_endpoint_rank':R[1],'kernel_extension':extra})
output={'status':'PASS','method':'independent SymPy DomainMatrix exact rational rank; rebuilt edge universe','rows':out}
Path(__file__).with_name('path-scatter-audit-receipt.json').write_text(json.dumps(output,indent=2)+'\n');print(output)
