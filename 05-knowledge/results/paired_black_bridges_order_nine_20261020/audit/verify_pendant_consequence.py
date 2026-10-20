"""Check the imported leaf-series algebra and its finite bridge consequences.

This does not re-audit the companion theorem D(z)=z^4.
"""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
COMPANION=ROOT.parent/'session15-deficit-centers'
EXPECTED='0123f01652018d17f62b6bf25f3acfa71fa743e55e798edb4daf1f7dfaf9a55c'
assert hashlib.sha256((COMPANION/'PUBLICATION_MANIFEST.json').read_bytes()).hexdigest()==EXPECTED
claim=json.loads((COMPANION/'leaf_series_coefficients.json').read_text())
N=claim['checked_through'];L=[]
num={4:2,5:-1,7:-1};den={1:-1,2:-1,5:1,6:-1,7:1,9:1}
for n in range(N+1):L.append(num.get(n,0)-sum(c*L[n-j] for j,c in den.items() if n>=j))
assert L==[r['all_leaves'] for r in claim['rows']]
P=[v-(n==4) for n,v in enumerate(L)];assert all(v>=0 for v in P)
finite=[]
for n in range(3,10):
 if n==9:
  data=json.loads((ROOT/'census_n9.json').read_text());count=len(data['paired_black']['bridges']);leaves=sum(r['degree']==1 for r in data['nodes'])
 else:
  old=ROOT.parent/'session14-metagraph-census'/('n8' if n==8 else '')/f'census_n{n}.json'
  data=json.loads(old.read_text());count=len(data['graphs']['paired_black']['bridges']);leaves=sum(r['reflection_pairs']==1 for r in data['nodes'])
 assert leaves==L[n] and count==P[n]
 finite.append({'n':n,'leaves':leaves,'pendant_bridges':P[n],'total_bridges':count,'nonleaf_bridges':0})
result={'imported_result':'D(z)=z^4; accepted from separately independently audited companion, not reproved here','companion_manifest_sha256':EXPECTED,'leaf_series_coefficients_checked_through':N,'pendant_series':'L(z)-z^4','total_bridge_series':'L(z)-z^4+E(z)','nonleaf_remainder':'E(z)=O(z^10), with nonnegative integer coefficients; not proved zero','finite_checks':finite,'all_checks':'passed'}
(ROOT/'audit/pendant_consequence_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print('Imported leaf-series algebra and pendant-bridge correction passed through 200 coefficients; exact bridge counts agree through n9.')
