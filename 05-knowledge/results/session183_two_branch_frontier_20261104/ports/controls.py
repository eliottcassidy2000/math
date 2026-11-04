from pathlib import Path
import json
from verify import subset,need
from automorphism import automorphisms
P=Path(__file__).resolve().parent
out=[]
for rows,order,response in zip(json.loads((P/'hostile_input.json').read_text())['rows'],[81,27],[1061775,1395279]):
 n=len(rows);h=subset(rows)[2][-1];child=[r|(((256>>i)&1)<<n) for i,r in enumerate(rows)]+[((1<<n)-1)^256]
 hc=subset(child)[2][-1];autos,nodes=automorphisms(rows)
 need(h==462915 and hc==response and len(autos)==order,'hostile control')
 # No deletion leaves a transitive tournament: every deletion retains a triangle.
 transitive_deletions=[]
 for deleted in range(n):
  rem=[v for v in range(n) if v!=deleted]
  cycle=any((rows[a]>>b&1) and (rows[b]>>c&1) and (rows[c]>>a&1) for a in rem for b in rem for c in rem)
  if not cycle:transitive_deletions.append(deleted)
 need(not transitive_deletions,'hostile is outside almost-transitive family')
 out.append(dict(H=h,birth_cut=256,birth_H=hc,full_automorphism_order=len(autos),automorphism_search_nodes=nodes,transitive_deletions=transitive_deletions))
(P/'controls.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
