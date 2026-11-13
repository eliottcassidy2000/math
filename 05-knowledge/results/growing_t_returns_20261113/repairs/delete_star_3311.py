from pathlib import Path
import importlib.util,json
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('owned_star_probe',HERE/'star_successor_probe.py');p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
T=p.T1691;cases=[]
for deleted in range(len(T)):
    keep=[u for u in range(len(T)) if u!=deleted];U=p.induced(T,keep);row=p.probe(U,3311);row['deleted_old_vertex']=deleted;row['retained_old_names']=keep;cases.append(row)
result=dict(status='PASS',source_rows=T,source_root=p.readyroot,target=3311,cases=cases)
(HERE/('delete_star_3311.json' if __debug__ else 'delete_star_3311.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
for row in cases:
    print(json.dumps(dict(deleted=row['deleted_old_vertex'],source_H=row['source_H'],choices=row['unique_start_star_choices'],hits=[dict(port_old=row['retained_old_names'][x['port']],star_flips_old=[row['retained_old_names'][u] for u in x['star_flips']],parent_state=x['parent_state'],strong=x['output_strong']) for x in row['hits']])))
