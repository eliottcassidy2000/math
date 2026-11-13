from pathlib import Path
import importlib.util,json,gzip
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parent.parent
p=RESULTS/'changing_core_returns_20261113/preparation/check.py';s=importlib.util.spec_from_file_location('jet_independent_counter',p);e=importlib.util.module_from_spec(s);s.loader.exec_module(e)
p=RESULTS/'berggren_target_return_20261111/reset/funnel_six_state.py';s=importlib.util.spec_from_file_location('jet_complete_compiler',p);f=importlib.util.module_from_spec(s);s.loader.exec_module(f)
head=json.loads(gzip.decompress((RESULTS/'multicore_metagraph_20261113/next_target.json.gz').read_bytes()))
def inflate(T,v,m):
    n=len(T);rows=list(T)+[0]*(m-1);module=[v]+list(range(n,n+m-1))
    for i in range(n):
        if i==v:continue
        for w in module[1:]:
            if T[i]>>v&1:rows[i]|=1<<w
            else:rows[w]|=1<<i
    for i,u in enumerate(module):
        for w in module[i+1:]:rows[u]|=1<<w
    return tuple(rows)
records=[]
for hit in head['hits']:
    if hit['H']!=3311:continue
    T,r=f.a.birth(f.a.unrank(hit['order'],hit['rank']),0)
    e.need(e.state(T,r)[0]==3311,'literal core H')
    for v in range(len(T)):
        V=inflate(T,v,2);H2=e.tables(V)[0][-1];e.need((H2-3311)%2==0,'integral K2 coefficient')
        rec=dict(parent_order=hit['order'],parent_rank=hit['rank'],order=len(T),ready_root=r,module_vertex=v,H2=H2,K2=(H2-3311)//2,K2_match=H2==9607)
        if rec['K2_match']:
            W=inflate(T,v,3);H3=e.tables(W)[0][-1];rec.update(H3=H3,K3=(H3-3311-6*3148)//6,K3_match=H3==27815)
        records.append(rec);e.tables.cache_clear()
result=dict(status='PASS',scope='Every resource vertex of the complete three strong unit-ready3311 cores; inherited complete finite head',target_coefficients=[3311,3148,936],vertices=len(records),K2_matches=[r for r in records if r['K2_match']],records=records)
(HERE/('jet_match.json' if __debug__ else 'jet_match.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(status='PASS',vertices=len(records),matches=result['K2_matches'],K2_values=sorted({r['K2'] for r in records}))))
