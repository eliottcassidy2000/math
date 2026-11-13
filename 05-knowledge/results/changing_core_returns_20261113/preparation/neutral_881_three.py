from pathlib import Path
from itertools import combinations
import importlib.util,json,time
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('independent_preparation_counts',HERE/'check.py');e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
spec=importlib.util.spec_from_file_location('audited_selector_wheel',HERE.parent/'selector/selector.py');selector=importlib.util.module_from_spec(spec);spec.loader.exec_module(selector)
source=json.loads((HERE.parent/'selector/receipt.json').read_text(encoding='utf-8'))['sources'][1]
T=tuple(source['rows']);root=source['root'];arcs=list(combinations(range(len(T)),2));endpoints=[];count=0;started=time.monotonic()
for patch in combinations(arcs,3):
    U=list(T)
    for a,b in patch:U[a]^=1<<b;U[b]^=1<<a
    U=tuple(U);H=e.tables(U)[0][-1]
    if H==881:
        z=e.state(U,root);wheel=selector.wheel(z)
        row=dict(arcs=patch,rows=U,state=z,strong=e.strong(U),root_indegree=sum(bool(U[v]>>root&1) for v in range(len(U))),ready_defect=z[2]-z[1],wheel=wheel,entry_possible=bool(wheel['good']))
        if wheel['good']:
            k=next(k for k in range(1,wheel['period']+1) if k%wheel['period'] in wheel['good'])
            row['first_good_branch']=k
            # Count a short physical return only; the wheel proves all branch lengths.
            if k<=3:
                V,r=e.returned(U,root,k);actual=e.state(V,r)
                e.need(actual==e.state_return(z,k),'independent literal return')
                row['literal_return']=dict(k=k,rows=V,root=r,state=actual)
        endpoints.append(row)
    e.tables.cache_clear();count+=1
    if count%5000==0:print(json.dumps(dict(done=count,H881=len(endpoints),elapsed=round(time.monotonic()-started,2))),flush=True)
result=dict(status='PASS',scope='All three-arc endpoint patches from the second named strong unit-ready H881 source, 55 choose3',patch_count=count,H881_endpoints=endpoints,entry_repairs=[r for r in endpoints if r['entry_possible']],count_engine='Independent subset-path engine from preparation/check.py',wheel_engine='Audited exact selector.wheel; successful short returns independently counted')
(HERE/('neutral_881_three.json' if __debug__ else 'neutral_881_three.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(status='PASS',patch_count=count,H881_endpoints=len(endpoints),entry_repairs=len(result['entry_repairs']),elapsed=round(time.monotonic()-started,2))))
