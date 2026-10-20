#!/usr/bin/env python3
from pathlib import Path
import re,json
from collections import Counter
ROOT=Path(__file__).resolve().parent;REPO=ROOT.parent/'math/05-knowledge/results'
REG=re.compile(r'n=\s*(\d+) (.*?): \|H\|=\s*(\d+) \|N\|=\s*(\d+) c_2=\s*(\d+) \|Fix\|=\s*(\d+) \| N-orbits: T\s*(\d+), N-even G\s*(\d+) \(R OK\) \| tau_T=\s*(\d+) tau_E=\s*(\d+) defect=\s*(\d+) \| W=\s*(\d+) E_T=\s*(\d+) E_even=\s*(\d+) E_wit=\s*(\d+)')
fields=['order','normalizer_order','pair_orbits','fixed_objects_each','relative_tournament_orbits','relative_even_graph_orbits','tau_T','tau_E','difference','sign_escape_global_types','tournament_fusion_defect','even_graph_fusion_defect','witness_excess']
mapnames={'F21_fixed':lambda n:f'F21 (7,1^{n-7})','C15_5plus3':lambda n:f'C15 (3,5,1^{n-8})','C21_7plus3':lambda n:f'C21 (3,7,1^{n-10})','F21_diagonal_7plus3':lambda n:f'F21 diag (7,3,1^{n-10})','F21_times_C3':lambda n:f'F21xC3 (7,3,1^{n-10})','F55_fixed':lambda n:f'F55 (11,1^{n-11})','C15_5plus3plus3':lambda n:f'C15 (3,3,5,1^{n-11}) synchronised','C5_times_C3_squared':lambda n:f'C3xC3xC5 (3,3,5,1^{n-11})','C35_7plus5':lambda n:f'C35 (5,7,1^{n-12})','F21_times_C5':lambda n:f'F21xC5 (7,5,1^{n-12})'}
incoming={}
for path in sorted((ROOT/'incoming').glob('*.out')):
 for line in path.read_text().splitlines():
  match=REG.search(line)
  if not match:continue
  n,name,*vs=match.groups();key=(int(n),name.strip());r=dict(zip(fields,map(int,vs)));r['sign_escape_relative_orbits']=r['sign_escape_global_types']+r['witness_excess'];incoming[key]=r
old={(r['n'],r['class']):r for r in json.loads((ROOT/'prior_reference_rows.json').read_text())['rows']}
comparisons=[]
new={}
for file in ('large13.json','small13.json'):
 for r in json.loads((ROOT/file).read_text())['rows']:new[r['n'],r['class']]=r
for label,refs in [('prior_audited',old),('new_independent',new)]:
 count=0
 for key,r in incoming.items():
  if key not in refs:continue
  expected=refs[key]
  for f,v in r.items():
   if f=='witness_excess':w=expected['sign_escape_relative_orbits']-expected['sign_escape_global_types']
   else:w=expected[f]
   if v!=w:raise RuntimeError(f'{label} {key} {f}: {v} != {w}')
  count+=1
 print(label,count)
 comparisons.append({'reference':label,'matching_rows':count,'all_reported_fields_match':True})
combined={k:r for k,r in old.items() if k[0]<=12};combined.update(new)
if len(incoming)!=44 or len(combined)!=48:raise RuntimeError('scope size')
summary=[]
for (n,name),r in sorted(combined.items()):
 summary.append({'n':n,'class':name,**{f:r[f] for f in fields if f!='witness_excess'},'sign_escape_relative_orbits':r['sign_escape_relative_orbits'],'witness_excess':r['sign_escape_relative_orbits']-r['sign_escape_global_types']})
(ROOT/'complete_census_summary.json').write_text(json.dumps({'status':'PASS','classified_rows':48,'incoming_computed_rows':44,'incoming_table_rows_including_borrowed':46,'new_degree13_rows':20,'class_counts_by_degree':dict(Counter(k[0] for k in combined)),'rows':summary},indent=2,sort_keys=True)+'\n')
print('total',len(summary),'failures',[(r['n'],r['class'],r['difference']) for r in summary if r['difference']])

(ROOT/'receipt_comparison.json').write_text(json.dumps({'status':'PASS','comparisons':comparisons,'incoming_rows':len(incoming),'complete_rows':len(combined)},indent=2)+'\n')
