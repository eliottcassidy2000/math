from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
states=[]
a=json.loads((HERE/'neutral_881_patch.json').read_text(encoding='utf-8'))
for source in a['sources']:
    for row in source['global_probe']['two_arc_H881_endpoints']:
        states.append(('width2_A' if source['parent_order']==9 else 'width2_B',row['arcs'],row['state']))
    for i,z in enumerate(source['all_partial_patch_states']):
        if z[0]==881:states.append(('five_patch_cube',i,z))
for row in json.loads((HERE/'neutral_881_three.json').read_text(encoding='utf-8'))['H881_endpoints']:
    states.append(('width3_B',row['arcs'],row['state']))
for row in json.loads((HERE/'neutral_881_four.json').read_text(encoding='utf-8'))['validated_H881_endpoints']:
    states.append(('width4_B',row['arcs'],row['state']))
found=[dict(source=tag,patch=patch,state=z,first_birth=z[0]+z[3]-z[1]) for tag,patch,z in states if z[3]-z[1]==810]
summary=dict(status='PASS',H881_endpoint_records=len(states),distinct_marked_profiles=len({tuple(z) for _,_,z in states}),target_first_birth=1691,matches=found,nearest=sorted({z[0]+z[3]-z[1] for _,_,z in states},key=lambda h:(abs(h-1691),h))[:10])
(HERE/'next_1691_probe.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(summary))
