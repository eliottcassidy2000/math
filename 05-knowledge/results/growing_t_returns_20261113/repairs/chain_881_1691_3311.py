from pathlib import Path
from functools import lru_cache
import importlib.util,json
HERE=Path(__file__).resolve().parent;p=HERE.parent.parent/'changing_core_returns_20261113/preparation/check.py';s=importlib.util.spec_from_file_location('chain_recount',p);e=importlib.util.module_from_spec(s);s.loader.exec_module(e)
A=(518,388,456,1011,487,455,387,257,1,502)
first=json.loads((HERE/'a2_alt298_two_star.raw.jsonl').read_text(encoding='utf-8-sig').splitlines()[0]);P=tuple(first['rows']);e.need(e.state(P,0)==(183,1,19,1509),'first parent independent state')
@lru_cache(None)
def paths(last,left):
    if not left:return 1
    return sum(paths(v,left^(1<<v)) for v in range(len(P)) if left>>v&1 and P[last]>>v&1)
e.need(paths(0,((1<<len(P))-1)^1)==1,'canonical unique start path')
order=[0];left=((1<<len(P))-1)^1
while left:
    last=order[-1];v=next(v for v in range(len(P)) if left>>v&1 and P[last]>>v&1 and paths(v,left^(1<<v)))
    order.append(v);left^=1<<v
T,r=e.birth(P,0,True);e.need(e.state(T,r)==(1691,182,183,1692) and e.strong(T),'first actual ready target')
canonical=json.loads((HERE/'alternate_1691_inputs.json').read_text(encoding='utf-8'))[0]['rows'];mapping=order+[r]
relabeled=tuple(sum(1<<j for j,v in enumerate(mapping) if T[u]>>v&1) for u in mapping)
e.need(relabeled==tuple(canonical),'physical canonical bridge retains all names')
second=json.loads((HERE/'alternate298_two_star.raw.jsonl').read_text(encoding='utf-8-sig').splitlines()[0]);Q=list(T);patch=[]
for a,b in second['changed_arcs']:
    u,v=mapping[a],mapping[b];Q[u]^=1<<v;Q[v]^=1<<u;patch.append((u,v))
Q=tuple(Q);port=mapping[second['root']];e.need(e.state(Q,port)==(205,1,65,3107),'second physical parent independent state')
W,wr=e.birth(Q,port,True);e.need(e.state(W,wr)==(3311,204,205,3312),'second actual ready target')
deleted=mapping[9];keep=[v for v in range(len(W)) if v!=deleted];V=tuple(sum(1<<j for j,v in enumerate(keep) if W[u]>>v&1) for u in keep);vr=keep.index(wr)
e.need(e.state(V,vr)==(3311,204,205,3312) and e.strong(V),'neutral sink deletion strong return')
result=dict(status='PASS',source_rows=A,first_patch=first['changed_arcs'],first_parent_rows=P,first_parent_state=e.state(P,0),first_output_rows=T,first_output_root=r,first_output_state=e.state(T,r),canonical_to_actual_names=mapping,second_patch=patch,second_parent_rows=Q,second_parent_root=port,second_parent_state=e.state(Q,port),second_output_rows=W,second_output_root=wr,second_output_state=e.state(W,wr),deleted_old_name=deleted,retained_names=keep,final_rows=V,final_root=vr,final_state=e.state(V,vr),scope='Literal two-step changing-core program; no all-k renewal claimed')
(HERE/('chain_881_1691_3311.json' if __debug__ else 'chain_881_1691_3311.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result))
