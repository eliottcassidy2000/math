from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
x=json.loads((HERE/'chain_881_1691_3311_6551.json').read_text(encoding='utf-8'))
stages=[('source881','source_rows',None,881),('helper1589','helper_rows',None,1589),('parent351','first_parent_rows',None,351),('ready1691','first_ready_rows',None,1691),('parent853','second_parent_rows',None,853),('ready3311','second_ready_rows',None,3311),('parent495','third_parent_rows',None,495),('ready6551_weak','third_ready_rows',None,6551),('ready6551_strong','final_rows',x['retained_names'],6551)]
records=[]
for name,key,names,H in stages:
 T=x[key];names=list(range(len(T))) if names is None else names;idx={old:i for i,old in enumerate(names)};p,q=idx[0],idx[2];exterior=[i for i in range(len(T)) if i not in(p,q)]
 order=sorted(exterior,key=lambda u:-sum(bool(T[u]>>v&1) for v in exterior))
 if not all(T[u]>>v&1 for i,u in enumerate(order) for v in order[i+1:]):raise RuntimeError('nontransitive residual')
 word=['%d%d'%((T[u]>>p)&1,(T[u]>>q)&1) for u in order]
 records.append(dict(stage=name,H=H,backbone_old_names=[names[u] for u in order],port_edge='0->2' if T[p]>>q&1 else '2->0',word=word))
for patchkey in ('first_patch','second_patch','third_patch'):
 if not all(0 in pair or 2 in pair for pair in x[patchkey]):raise RuntimeError('repair leaves fixed two-star menu')
result=dict(status='PASS',ports=[0,2],letter_definition='([backbone_vertex->0],[backbone_vertex->2])',records=records,scope='Exact finite four-letter trace; no uniform selector for later ray targets')
(HERE/('two_port_trace.json' if __debug__ else 'two_port_trace.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
for row in records:print(row['stage'],row['port_edge'],' '.join(row['word']))
