#!/usr/bin/env python3
"""Normalize audited JSON and build a portable atlas. No repository writes."""
import json, hashlib
from pathlib import Path
from collections import Counter
P=Path(__file__).resolve().parent
W=P.parent/'session29-target-wiring'; E=P.parent/'session29-bridge-extension'
read=lambda p:json.loads(p.read_text())
w=read(W/'wiring-results.json'); c=read(W/'cycle-results.json'); a=read(W/'all-length-results.json'); e=read(E/'verify-normal.json')
checks=[]
def check(name,condition):
 if not condition: raise ValueError(name)
 checks.append(name)
rows=w['witness_rows']; targets=w['incidence_table']; sources=['AA','AB/BA','BB']; leaves={tuple(x['target_key']) for x in c['restricted_multigraph_bridges']}
for t in targets:
 rs=[r for r in rows if r['merged_key']==t['root_pair']]
 aut={r['target_aut_order'] for r in rs}; check('consistent target '+str(t['root_pair']),len(aut)==1 and all(r['H']==t['H'] for r in rs))
 t['aut']=aut.pop(); t['degree']=(1 if t['self_converse'] else 2)*t['H']//t['aut']; t['residual']=t['degree']-t['total']; t['leaf']=tuple(t['root_pair']) in leaves
 t['witnesses']=[{'edge':r['row'],'source':r['source'],'tiles':r['tiles'],'blue':r['blue']} for r in rs]
 for s in sources:
  for color in ['blue','black']:
   check('incidence '+str(t['root_pair'])+s+color,t['incidences'][s][color]==sum((r['source'] if r['source'] in ['AA','BB'] else 'AB/BA')==s and r['blue']==(color=='blue') for r in rs))
 check('integral degree '+str(t['root_pair']),(1 if t['self_converse'] else 2)*t['H']%t['aut']==0 and t['residual']>=0)
check('13 roots / 91 targets',len(w['root_classes'])==13 and len(targets)==91)
check('324 lines / 18 blue',len(rows)==324 and sum(r['blue'] for r in rows)==18)
check('source masses',[sum(t['incidences'][s][k] for t in targets for k in ['blue','black']) for s in sources]==[25,130,169])
check('source blue masses',[sum(t['incidences'][s]['blue'] for t in targets) for s in sources]==[5,0,13])
check('13 self-converse / 8 leaves',sum(t['self_converse'] for t in targets)==13 and sum(t['leaf'] for t in targets)==8)
check('H moments',sum(r['H'] for r in rows)==195472 and sum(r['H'] for r in rows if r['blue'])==11522)
check('H range',(min(t['H'] for t in targets),max(t['H'] for t in targets))==(221,1049))
check('all-r alphabet',(a['a'],a['b'],a['f'])==(13,17,5))
check('324 global nonbridges',c['full_metagraph_certified_nonbridge_edges']+e['verified_original_global_nonbridges']==324)
data={'roots':w['root_classes'],'targets':targets,'sources':sources,'status':{'original_boundary_leaves':8,'certified_global_nonbridges':324},'alphabet':[a['a'],a['b'],a['f']]}
(P/'atlas-data.json').write_text(json.dumps(data,indent=2)+'\n')
html=(P/'template.html').read_text().replace('__DATA__',json.dumps(data,separators=(',',':')))
(P/'atlas.html').write_text(html)
report={'passed':True,'checks':len(checks),'targets':91,'edges':324,'blue':18,'black':306,'global_nonbridges':324,'full_merged_target_degree_sum':sum(t['degree'] for t in targets),'captured_incidence_sum':324,'source_sha256':{str(p.relative_to(P.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [W/'wiring-results.json',W/'cycle-results.json',W/'all-length-results.json',E/'verify-normal.json']}}
(P/'validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
