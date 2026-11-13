from pathlib import Path
from math import comb
import importlib.util,json,shutil,subprocess,tempfile
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('four_independent_counts',HERE/'check.py');e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
spec=importlib.util.spec_from_file_location('four_audited_wheel',HERE.parent/'selector/selector.py');selector=importlib.util.module_from_spec(spec);spec.loader.exec_module(selector)
source=json.loads((HERE.parent/'selector/receipt.json').read_text(encoding='utf-8'))['sources'][1]
T=tuple(source['rows']);root=source['root'];input_rows=' '.join(map(str,T))+'\n'
compiler=shutil.which('g++');e.need(compiler is not None,'g++ compiler available')
with tempfile.TemporaryDirectory(prefix='chan_pak_881_four_') as build:
    executable=Path(build)/'counter.exe'
    subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'neutral_881_four.cpp'),'-o',str(executable)],check=True,capture_output=True,text=True)
    streams={width:subprocess.run([str(executable),str(width)],input=input_rows,check=True,capture_output=True,text=True).stdout for width in (3,4)}
receipts={width:[json.loads(line) for line in text.splitlines()] for width,text in streams.items()}
old=json.loads((HERE/'neutral_881_three.json').read_text(encoding='utf-8'))
def canonical(rows):return sorted((tuple(tuple(a) for a in x['arcs']),tuple(x['state']),x.get('indegree',x.get('root_indegree')),x['strong']) for x in rows)
old_filtered=[x for x in old['H881_endpoints'] if x['root_indegree']>=2]
e.need(canonical(receipts[3][:-1])==canonical(old_filtered),'independent C++ three-arc universe equals exhaustive Python scan')
for width in (3,4):e.need(receipts[width][-1]['all_arc_patches']==comb(55,width),'complete patch enumeration')
validated=[]
for row in receipts[4][:-1]:
    U=list(T)
    for a,b in row['arcs']:U[a]^=1<<b;U[b]^=1<<a
    U=tuple(U);e.need(e.tournament(U),'reconstructed tournament')
    z=e.state(U,root);e.need(z==tuple(row['state']),'independent full C++ candidate state')
    e.need(e.strong(U)==row['strong'],'independent strongness')
    indegree=sum(bool(U[v]>>root&1) for v in range(len(U)));e.need(indegree==row['indegree'],'root degree')
    wheel=selector.wheel(z);record=dict(row,rows=U,wheel=wheel,entry_possible=bool(wheel['good']))
    if wheel['good']:
        reach=[U[v]|(1<<v) for v in range(len(U))]
        for v in range(len(U)):
            for a in range(len(U)):
                if reach[a]>>v&1:reach[a]|=reach[v]
        core=[v for v in range(len(U)) if (reach[root]>>v&1) and (reach[v]>>root&1)];tail=[v for v in range(len(U)) if v not in core]
        V=tuple(sum(1<<j for j,v in enumerate(core) if U[u]>>v&1) for u in core);rr=core.index(root)
        tail_rows=tuple(sum(1<<j for j,v in enumerate(tail) if U[u]>>v&1) for u in tail)
        e.need(all(U[u]>>v&1 for u in core for v in tail) and e.tables(tail_rows)[0][-1]==1,'dominated transitive tail')
        e.need(e.strong(V) and e.state(V,rr)==z,'strong core preserves complete profile')
        k=next(k for k in range(1,wheel['period']+1) if k%wheel['period'] in wheel['good']);e.need(k==1,'all selected first exits are RNR')
        W,wr=e.returned(V,rr,k);actual=e.state(W,wr);e.need(actual==e.state_return(z,k) and e.strong(W),'literal strong core return')
        record['core_repair']=dict(retained_old_names=core,deleted_old_names=tail,root_old_name=root,root_new_index=rr,rows=V,state=z,literal_return=dict(k=k,rows=W,root=wr,state=actual))
    validated.append(record);e.tables.cache_clear()
best=next(x for x in validated if x['arcs']==[[0,9],[2,6],[2,10],[3,6]])
e.need(best['core_repair']['rows']==(354,228,425,179,231,192,140,1,250),'selected nine-vertex core literal')
e.need(best['state']==[881,169,201,977] and best['wheel']['period']==2 and best['wheel']['good']==[1],'selected core exact wheel')
result=dict(status='PASS',scope='Complete width4 search with necessary indegree>=2 filter; all matching states independently recounted; width3 counter independently agrees with full Python scan',input_rows=T,root=root,three_arc_cpp_crosscheck=receipts[3][-1],four_arc_cpp_counts=receipts[4][-1],validated_H881_endpoints=validated,entry_repair_count=sum(x['entry_possible'] for x in validated),strong_endpoint_entry_repairs=sum(x['entry_possible'] and x['strong'] for x in validated),selected_nine_vertex_core=best)
(HERE/'neutral_881_four.raw.jsonl').write_text(streams[4],encoding='utf-8',newline='\n')
(HERE/('neutral_881_four.json' if __debug__ else 'neutral_881_four.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(status='PASS',three_arc_crosscheck=receipts[3][-1],four_arc_counts=receipts[4][-1],entry_repair_count=result['entry_repair_count'],strong_endpoint_entry_repairs=result['strong_endpoint_entry_repairs'],selected_core_state=best['state'],literal_exit=best['core_repair']['literal_return']['state'])))
