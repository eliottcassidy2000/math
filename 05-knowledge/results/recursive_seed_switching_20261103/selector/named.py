"""Literal same-order/clock repair: retain old seed names, promote two old tails."""
from pathlib import Path
import sys,json,hashlib,itertools
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
OLD=RESULTS/'zero_ready_component_returns_20261103'
sys.path.insert(0,str(OLD/'inverse'));import reader
from reader import need

def main():
    tree=json.loads((OLD/'continuation/named_tree.json').read_text())
    parent=next(x for x in tree['targets'] if x['H']==4611815)
    source_seed=parent['seed_reader'];old_seed=tuple(source_seed['rows']);old_seed_names=source_seed['names']
    need(old_seed==(38,60,24,49,33,4) and old_seed_names==[6,4,7,8,14,31],'Actual source seed and names')
    rows=reader.birth(reader.birth(old_seed,37),55)
    need(rows==(38,124,24,241,97,4,165,55),'Two local seed births')
    root=1;word='1111000100';prefix=list(parent['prefix']);tail=parent['component_names'][len(old_seed):]
    promoted=tail[:2];seed_names=old_seed_names+promoted;tail_names=tail[2:]
    need(prefix==[50,51,1,2] and len(tail_names)==len(word),'Same physical prefix/clock')
    component_names=seed_names+tail_names;component=reader.graph(rows,root,word);fields=reader.all_fields(rows)[root]
    output=reader.run(fields,word);target=3*parent['H']+2
    need(output[0]+(1<<len(prefix))*output[3]==target==13835447,'Exact target')
    names=list(parent['names']);edges=set()
    for i,a in enumerate(component_names):
        for j,b in enumerate(component_names):
            if component[i]>>j&1:edges.add((a,b))
    for i,a in enumerate(prefix):
        for b in prefix[i+1:]+component_names:edges.add((a,b))
    mark=seed_names[root]
    for v in names:
        if v!=39:edges.add((mark,39) if v==mark else (39,v))
    full=tuple(sum(1<<j for j,b in enumerate(names) if (a,b) in edges) for a in names)
    need(reader.strong(component) and reader.strong(full),'Strong component and full target')
    loc={v:i for i,v in enumerate(names)}
    changes=[[a,b] for i,a in enumerate(names) for b in names[i+1:] if bool(parent['rows'][i]>>loc[b]&1)!=bool(full[i]>>loc[b]&1)]
    need(all(a not in prefix and b not in prefix for a,b in changes),'All prefix incidence preserved')
    need(all([a,b] not in changes and [b,a] not in changes for i,a in enumerate(old_seed_names) for b in old_seed_names[i+1:]),'Every old6seed incidence preserved literally')
    required={39,mark,31,6};extras=[v for v in seed_names if v not in required];feedback=None
    for k in range(len(extras)+1):
        for add in itertools.combinations(extras,k):
            choice=required|set(add);rest=[loc[v] for v in names if v not in choice]
            chain=sorted(rest,key=lambda v:-sum(full[v]>>u&1 for u in rest))
            if all(full[u]>>v&1 for i,u in enumerate(chain) for v in chain[i+1:]):feedback=[v for v in names if v in choice];break
        if feedback is not None:break
    need(feedback is not None,'A structural feedback deletion containing retained ports')
    target_record=dict(H=target,parent_H=parent['H'],offset=2,names=names,rows=full,root=39,mark=mark,previous_mark=31,spare=6,prefix=prefix,
        component_names=component_names,component_rows=component,seed_reader=dict(rows=rows,root=root,names=seed_names,fields=fields,word=word),component_eight_fields=output,
        structural_feedback_names=feedback,promoted_tail_names=promoted,new_names=[],deleted_old_names=[],old_name_identity=[[v,v] for v in names],old_arc_differences=changes,
        seed_birth_recipe=[dict(name=promoted[0],newborn_out_mask=37),dict(name=promoted[1],newborn_out_mask=55)],
        seed_recipe_scope='These vertices are fresh relative to the6-vertex seed only; both are existing named full-source tail vertices. The full graph changes by the recorded old-arc rewrites.')
    pins=[OLD/'continuation/named_tree.json',OLD/'inverse/reader.py',Path(__file__)]
    report=dict(status='FINITE_EXACT_NAMED_READER_TEMPLATE',scope='One fixed-bank miss repaired by local seed enlargement and selecting the existing physical vertex4; no unlimited closure.',source=parent,source_H=parent['H'],source_names=names,source_rows=parent['rows'],target=target_record,source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    (HERE/'named_repair.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(H=target,order=len(names),mark=mark,previous_mark=31,prefix=prefix,promoted=promoted,changed_old_arcs=len(changes),seed_names=seed_names)),flush=True)
if __name__=='__main__':main()
