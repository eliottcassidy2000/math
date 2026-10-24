"""Recover complete actual domains of the inherited28-source two-target supply."""
import csv,gzip,json,sys
from collections import Counter
from itertools import combinations,permutations
from pathlib import Path
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
sys.path.insert(0,str(RESULTS/'injective_native_blocks_20261023/search'))
import explore as e
from check import adjacency,mask_on
TARGETS=(1965,2234)
def main():
    incoming=json.loads((RESULTS/'session45_injective_multisource_moves_20261023/source-circuits/all_menus.json').read_text())
    incoming=[m for m in incoming if tuple(m['targets'])==TARGETS]
    out=[]
    for menu in incoming:
        s=menu['source'];best=menu['best'];p=best['before']['path'];q=best['after']['path'];pins=best['pins']
        first,last,k,a,b=pins;fixed={i:u for i,(u,v) in enumerate(zip(p,q)) if u==v}
        rows=e.paths(s)
        arc=[r for r in rows if r['path'][0]==first and r['path'][-1]==last and r['path'][k:k+2]==[a,b] and r['target'] in TARGETS]
        raw_full=[r for r in rows if all(r['path'][i]==u for i,u in fixed.items())]
        full=[r for r in raw_full if r['target'] in TARGETS]
        any_target=[r for r in rows if r['target'] in TARGETS]
        out.append(dict(source=s,H=len(rows),X=e.X[s],pins=pins,fixed_clocks=sorted(fixed.items()),old_before=p,old_after=q,
                        exact_arc_target_domain=arc,exact_all_shared_clocks_target_domain=full,
                        raw_all_shared_clocks_domain=raw_full,
                        unpinned_target_domain=any_target,
                        code=[e.X[s]%128]+[e.Y[j][s]%2 for j in range(1,5)]))
    receipt=dict(targets=TARGETS,source_count=len(out),sources=out,
                 arc_domain_sizes=dict(Counter(len(z['exact_arc_target_domain']) for z in out)),
                 all_shared_domain_sizes=dict(Counter(len(z['exact_all_shared_clocks_target_domain']) for z in out)),
                 raw_all_shared_domain_sizes=dict(Counter(len(z['raw_all_shared_clocks_domain']) for z in out)),
                 unpinned_domain_sizes=dict(Counter(len(z['unpinned_target_domain']) for z in out)))
    (HERE/'domains.json').write_text(json.dumps(receipt,indent=2))
    print({k:v for k,v in receipt.items() if k!='sources'})
    for z in out:
        if len(z['exact_arc_target_domain'])>2 or len(z['exact_all_shared_clocks_target_domain'])>2:
            print(z['source'],len(z['exact_arc_target_domain']),len(z['exact_all_shared_clocks_target_domain']),Counter(r['target'] for r in z['exact_arc_target_domain']))
if __name__=='__main__':main()
