from pathlib import Path
from functools import lru_cache
import importlib.util,json
HERE=Path(__file__).resolve().parent
p=HERE.parent.parent/'changing_core_returns_20261113/preparation/check.py';spec=importlib.util.spec_from_file_location('three_step_independent_counter',p);e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
def hit(name):return json.loads((HERE/name).read_text(encoding='utf-8-sig').splitlines()[0])
def apply(T,patch):
    V=list(T)
    for a,b in patch:V[a]^=1<<b;V[b]^=1<<a
    return tuple(V)
def induced(T,keep):return tuple(sum(1<<j for j,v in enumerate(keep) if T[u]>>v&1) for u in keep)
def unique_path(T,root):
    @lru_cache(None)
    def count(last,left):
        if not left:return 1
        return sum(count(v,left^(1<<v)) for v in range(len(T)) if left>>v&1 and T[last]>>v&1)
    left=((1<<len(T))-1)^(1<<root);e.need(count(root,left)==1,'unique path for named isomorphism');path=[root]
    while left:
        last=path[-1];v=next(v for v in range(len(T)) if left>>v&1 and T[last]>>v&1 and count(v,left^(1<<v)))
        path.append(v);left^=1<<v
    return path
A=(518,388,456,1011,487,455,387,257,1,502);helper,hr=e.birth(A,3,True)
e.need(e.state(helper,hr)[0]==1589,'actual helper birth')
h1=hit('helper_alt303_two_star.raw.jsonl');P=apply(helper,h1['changed_arcs']);e.need(P==tuple(h1['rows']) and e.state(P,0)==(351,1,3,1341) and e.strong(P),'first parent literal')
T1,r1=e.birth(P,0,True);e.need(e.state(T1,r1)==(1691,350,351,1692) and e.strong(T1),'first strong ready target')
mapping=unique_path(P,0)+[r1];canonical=json.loads((HERE/'alternate_1691_inputs.json').read_text(encoding='utf-8'))[2]['rows'];e.need(induced(T1,mapping)==tuple(canonical),'actual name map into alternate303')
h2=hit('alternate303_two_star.raw.jsonl');patch2=[(mapping[a],mapping[b]) for a,b in h2['changed_arcs']];Q=apply(T1,patch2);port=mapping[h2['root']]
e.need(e.state(Q,port)==(853,1,5,2459) and e.strong(Q),'second parent literal')
T2,r2=e.birth(Q,port,True);e.need(e.state(T2,r2)==(3311,852,853,3312) and e.strong(T2),'second strong ready target');mapping2=mapping+[r2]
h3=hit('alternate6551_two_star.raw.jsonl');patch3=[(mapping2[a],mapping2[b]) for a,b in h3['changed_arcs']];R=apply(T2,patch3);port3=mapping2[h3['root']]
e.need(e.state(R,port3)==(495,1,33,6057),'third parent literal')
T3,r3=e.birth(R,port3,True);e.need(e.state(T3,r3)==(6551,494,495,6552),'third ready target');deleted=mapping2[10]
e.need(T3[deleted]==0,'named final sink');keep=[v for v in range(len(T3)) if v!=deleted];final=induced(T3,keep);finalroot=keep.index(r3)
e.need(e.state(final,finalroot)==(6551,494,495,6552) and e.strong(final),'strong6551 after neutral sink deletion')
result=dict(status='PASS',scope='Actual finite three-step target chain with one helper; no all-k closure',source_rows=A,helper_port=3,helper_rows=helper,helper_state=e.state(helper,hr),first_patch=h1['changed_arcs'],first_parent_rows=P,first_parent_state=e.state(P,0),first_ready_rows=T1,first_ready_root=r1,first_ready_state=e.state(T1,r1),canonical_to_actual_names=mapping,second_patch=patch2,second_parent_rows=Q,second_parent_root=port,second_parent_state=e.state(Q,port),second_ready_rows=T2,second_ready_root=r2,second_ready_state=e.state(T2,r2),third_patch=patch3,third_parent_rows=R,third_parent_root=port3,third_parent_state=e.state(R,port3),third_ready_rows=T3,third_ready_root=r3,third_ready_state=e.state(T3,r3),deleted_old_name=deleted,retained_names=keep,final_rows=final,final_root=finalroot,final_state=e.state(final,finalroot))
(HERE/('chain_881_1691_3311_6551.json' if __debug__ else 'chain_881_1691_3311_6551.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(status='PASS',helper_H=result['helper_state'][0],ready_states=[result['first_ready_state'],result['second_ready_state'],result['final_state']],deleted=deleted,final_order=len(final),final_rows=final)))
