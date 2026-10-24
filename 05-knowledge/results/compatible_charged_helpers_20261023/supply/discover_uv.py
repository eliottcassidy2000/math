"""New U/V common-clock menus beyond the inherited28; discovery only."""
from discover import *
def main():
    adj={t:set() for t in (1965,2234)}
    with gzip.open(RESULTS/'metagraph_order8_20261023/merged_edges.tsv.gz','rt') as f:
        for z in csv.DictReader(f,delimiter='\t'):
            a,b=int(z['source']),int(z['target'])
            if a in adj:adj[a].add(b)
            if b in adj:adj[b].add(a)
    old=json.loads((RESULTS/'native_packet_decomposition_20261023/native/domains.json').read_text())['sources']
    oldids={r['source'] for r in old}; used={min(orbit(r['mask'])) for z in old for r in z['exact_all_shared_clocks_target_domain']}
    newsources=sorted((adj[1965]&adj[2234])-oldids)
    menus=[]
    for s in newsources:
        rows=e.paths(s);targetrows=[r for r in rows if r['target'] in (1965,2234)]
        for r in targetrows:r['orbit']=orbit(r['mask']);r['resource']=min(r['orbit'])
        moves=[]
        for a,b in product([r for r in targetrows if r['target']==1965],[r for r in targetrows if r['target']==2234]):
            fixed=[[j,x] for j,(x,y) in enumerate(zip(a['path'],b['path'])) if x==y]
            full=[r for r in rows if all(r['path'][j]==x for j,x in fixed)]
            retained=[r for r in full if r['target'] in (1965,2234)]
            moves.append(dict(before=a,after=b,common_clocks=fixed,cost=8-len(fixed),
                complete_pinned_count=len(full),target_pinned_count=len(retained),
                target_pinned_masks=[r['mask'] for r in retained],pins_target_count=dict(Counter(r['target'] for r in full)),
                old_resource_collision=bool({a['resource'],b['resource']}&used),same_resource=a['resource']==b['resource']))
        good=[m for m in moves if not m['old_resource_collision'] and not m['same_resource']]
        if not good:continue
        best=min(good,key=lambda m:(m['cost'],m['target_pinned_count'],m['complete_pinned_count'],m['before']['path'],m['after']['path']))
        menus.append(dict(source=s,canon=int(e.NODES[s]['canon_a']),H=len(rows),X=e.X[s],label=label(s),
                         full_UV_count=len(targetrows),moves=moves,best=best))
    oldlabels={s:label(s) for s in oldids}
    for m in menus:m['equal_charge_old_sources']=[s for s,v in oldlabels.items() if v==m['label']]
    result=dict(targets=[1965,2234],old_sources=sorted(oldids),new_candidates=newsources,menus=menus,
        summary=dict(new_sources=len(menus),all_H=sum(m['H'] for m in menus),
        best_cost_histogram=dict(Counter(m['best']['cost'] for m in menus)),
        complete_binary_target_menus=sum(m['best']['target_pinned_count']==2 for m in menus),
        raw_binary_menus=sum(m['best']['complete_pinned_count']==2 for m in menus),
        old_charge_matches=sum(bool(m['equal_charge_old_sources']) for m in menus)))
    (HERE/'uv_discovery.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result['summary'],indent=2))
    for m in sorted(menus,key=lambda m:(m['best']['cost'],m['source']))[:20]:print(m['source'],m['best']['cost'],m['best']['common_clocks'],m['label'],m['equal_charge_old_sources'])
if __name__=='__main__':main()
