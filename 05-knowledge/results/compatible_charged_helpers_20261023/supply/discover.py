"""Bounded discovery from the existing atlas; only common-neighbor fibres."""
from collections import defaultdict,Counter
from itertools import combinations,product
from pathlib import Path
import csv,gzip,json,sys
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parents[1]
sys.path.insert(0,str(RESULTS/'injective_native_blocks_20261023/search'))
import explore as e
TARGETS=(570,1965,2234)
REV=[e.CHORDS.index((7-j,7-i)) for i,j in e.CHORDS]
def orbit(mask):
    rev=sum(((mask>>j)&1)<<REV[j] for j in range(21))
    return sorted({mask,rev,e.FULL^mask,e.FULL^rev})
def label(s):return [e.X[s]%128]+[e.Y[j][s]%2 for j in range(1,5)]
def add(a,b):return tuple((x+y)%m for x,y,m in zip(a,b,(128,2,2,2,2)))
def main():
    adj={t:{} for t in TARGETS}
    with gzip.open(RESULTS/'metagraph_order8_20261023/merged_edges.tsv.gz','rt') as f:
        for z in csv.DictReader(f,delimiter='\t'):
            a,b=int(z['source']),int(z['target']);weights=[int(z['blue']),int(z['black'])]
            if a in adj:adj[a][b]=weights
            if b in adj:adj[b][a]=weights
    candidates=sorted((set(adj[570])&set(adj[1965]))|(set(adj[570])&set(adj[2234])))
    old=json.loads((RESULTS/'native_packet_decomposition_20261023/native/domains.json').read_text())['sources']
    oldids={r['source'] for r in old}
    used={min(orbit(r['mask'])) for z in old for r in z['exact_all_shared_clocks_target_domain']}
    records=[];moves=[]
    for s in candidates:
        rows=e.paths(s)
        targetrows=[r for r in rows if r['target'] in TARGETS]
        for r in targetrows:r['orbit']=orbit(r['mask']);r['resource']=min(r['orbit'])
        rr=dict(source=s,H=len(rows),old_source=s in oldids,label=label(s),target_paths=targetrows,moves=[])
        for u in (1965,2234):
            for a,b in product([r for r in targetrows if r['target']==570],[r for r in targetrows if r['target']==u]):
                fixed=[[j,x] for j,(x,y) in enumerate(zip(a['path'],b['path'])) if x==y]
                complete=[r for r in rows if all(r['path'][j]==x for j,x in fixed)]
                retained=[r for r in complete if r['target'] in (570,u)]
                move=dict(source=s,targets=[570,u],before=a,after=b,common_clocks=fixed,
                    cost=8-len(fixed),complete_pinned_count=len(complete),target_pinned_count=len(retained),
                    target_pinned_masks=[r['mask'] for r in retained],pins_target_count=dict(Counter(r['target'] for r in complete)),
                    old_resource_collision=bool({a['resource'],b['resource']}&used),same_resource=a['resource']==b['resource'])
                rr['moves'].append(move);moves.append(move)
        records.append(rr)
    newmoves=[m for m in moves if m['source'] not in oldids and not m['old_resource_collision'] and not m['same_resource']]
    circuits=[]
    for u in (1965,2234):
        bysource=defaultdict(list)
        for m in newmoves:
            if m['targets'][1]==u:bysource[m['source']].append(m)
        best={s:min(ms,key=lambda m:(m['cost'],m['target_pinned_count'],m['before']['path'],m['after']['path'])) for s,ms in bysource.items()}
        bins=defaultdict(list)
        for a,b in combinations(sorted(best),2):bins[add(label(a),label(b))].append((a,b))
        for pairs in bins.values():
            for p,q in combinations(pairs,2):
                if set(p)&set(q):continue
                slots=[best[s] for s in p+q];resources=[m[z]['resource'] for m in slots for z in ('before','after')]
                if len(set(resources))!=8:continue
                circuits.append(dict(targets=[570,u],positive=list(p),negative=list(q),cost=sum(m['cost'] for m in slots),slots=slots,
                                     integer_X_debt=sum(e.X[s] for s in p)-sum(e.X[s] for s in q)))
    circuits.sort(key=lambda c:(c['cost'],abs(c['integer_X_debt']),c['targets'],c['positive'],c['negative']))
    result=dict(targets=TARGETS,target_labels={t:label(t) for t in TARGETS},X={t:e.X[t] for t in TARGETS},
        common_neighbor_sources=candidates,old_intersection=sorted(set(candidates)&oldids),old_resources=len(used),
        records=records,new_moves=newmoves,circuits=circuits,
        summary=dict(candidate_sources=len(candidates),new_sources=len(set(candidates)-oldids),all_H=sum(r['H'] for r in records),
            new_moves=len(newmoves),cost_histogram=dict(Counter(m['cost'] for m in newmoves)),four_source_circuits=len(circuits)))
    (HERE/'discovery.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('target_labels','X','old_intersection','summary')},indent=2))
    for c in circuits[:8]:print('circuit',c['targets'],c['positive'],c['negative'],'cost',c['cost'],'integerdebt',c['integer_X_debt'])
    print('best new moves')
    for m in sorted(newmoves,key=lambda m:(m['cost'],m['source']))[:10]: print(m['source'],m['targets'],m['cost'],m['common_clocks'],m['target_pinned_count'])
if __name__=='__main__':main()
