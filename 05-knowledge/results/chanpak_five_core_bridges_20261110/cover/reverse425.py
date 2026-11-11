"""Reverse one-birth probes in exactly the inherited normalized target inverse."""
import importlib.util
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'session251_creative_completion_20261111/ternary/compiler.py'
spec=importlib.util.spec_from_file_location('inherited',SOURCE)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
CORES=((0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11))

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def induced(rows,ids):
    return tuple(sum(((rows[v]>>w)&1)<<j for j,w in enumerate(ids)) for v in ids)

def strong(rows):
    full=(1<<len(rows))-1
    for v in range(len(rows)):
        seen=1<<v
        while True:
            after=seen
            for w,r in enumerate(rows):
                if seen>>w&1:after|=r
            if after==seen:break
            seen=after
        if seen!=full:return False
    return True

def main():
    graphs={};raw=0
    for ci,core in enumerate(CORES):
        for cut in range(1<<len(core)):
            for value,A,B in c.generate(c.signature(core,cut),425):
                if value!=425:continue
                raw+=1;rows=tuple(c.construct(core,cut,A,B))
                graphs.setdefault(rows,(ci,cut,A,B))
    hits=[]
    for rows,params in sorted(graphs.items(),key=lambda kv:(len(kv[0]),kv[0])):
        hs=c.hp_subsets(rows);full=(1<<len(rows))-1
        need(hs[-1]==425,'actual target')
        for v in range(len(rows)):
            if hs[full^(1<<v)]!=141:continue
            ids=[w for w in range(len(rows)) if w!=v];source=induced(rows,ids)
            source_full=(1<<len(source))-1
            global_sources=[j for j,r in enumerate(source) if r==source_full^(1<<j)]
            global_sinks=[j for j,r in enumerate(source) if r==0]
            neutral=global_sources+global_sinks
            need(len(neutral)==1,'exact neutral padding')
            strong_ids=[j for j in range(len(source)) if j not in neutral]
            ancestor=induced(source,strong_ids)
            O=sum(((rows[v]>>w)&1)<<j for j,w in enumerate(ids))
            child=tuple(r|(0 if O>>j&1 else 1<<len(source)) for j,r in enumerate(source))+(O,)
            need(child==induced(rows,ids+[v]),'explicit reinsertion labelling')
            need(c.hp_subsets(ancestor)[-1]==141 and strong(ancestor),'strong ancestor')
            need(not strong(source) and strong(rows),'weak/strong boundary')
            hits.append({'core_index':params[0],'core_rows':CORES[params[0]],
                         'cut':params[1],'A':params[2],'B':params[3],
                         'target_rows':rows,'deleted_target_vertex':v,
                         'source_vertex_labels_in_target':ids,'source_rows':source,
                         'new_outmask_in_source_labels':O,
                         'source_h':141,'source_is_strong':False,
                         'global_source_vertices':global_sources,'global_sink_vertices':global_sinks,
                         'strong_ancestor_vertices_in_source':strong_ids,'strong_ancestor_rows':ancestor,
                         'trace':[141,141,425]})
    need(len(graphs)==42 and len(hits)==6,'normalized witness scope')
    return {'status':'PASS','raw_inverse_witnesses':raw,'distinct_literal_targets':len(graphs),
            'one_birth_hits':hits,
            'scope':'Only distinct literal target425 graphs produced by the inherited normalized finite inverse over90 raw cuts. Normalization is complete for values, not for every graph embedding. All six H141 parents here are weak neutral paddings; no new strong one-birth claim.'}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
