"""Complete stable-twin contraction graph for the declared strong bank <=8."""
from pathlib import Path
from itertools import permutations,product
from collections import defaultdict,Counter
import json
HERE=Path(__file__).resolve().parent
def need(x,msg):
    if not x:raise RuntimeError(msg)
def key(rows,mask):
    cells=defaultdict(list)
    for i,r in enumerate(rows):cells[(mask>>i&1,int(r).bit_count())].append(i)
    groups=[cells[c] for c in sorted(cells)]
    best=None
    for chunks in product(*(permutations(g) for g in groups)):
        order=sum((tuple(p) for p in chunks),())
        newrows=tuple(sum(1<<j for j,v in enumerate(order) if rows[u]>>v&1) for u in order)
        mark=sum(1<<j for j,v in enumerate(order) if mask>>v&1)
        value=(len(rows),mark,*newrows)
        if best is None or value<best:best=value
    return best
def run():
    source=json.loads((HERE/'cyclic_port_receipt.json').read_text())
    records=[r for r in source['positive_screen_certificates'] if r['individually_stable']]
    ids={key(r['rows'],r['selected']):i for i,r in enumerate(records)}
    need(len(ids)==len(records),'marked canonical keys')
    roots={};depth={};out=[];edges=0
    for i,r in enumerate(records):
        rows=r['rows'];mask=r['selected'];children=set()
        for a,b in r['stable_twin_contractions']:
            kept=[v for v in range(len(rows)) if v!=a]
            child=[sum(1<<j for j,v in enumerate(kept) if rows[u]>>v&1) for u in kept]
            selected=sum(1<<j for j,v in enumerate(kept) if mask>>v&1)
            k=key(child,selected);need(k in ids,'contracted state absent from complete bank')
            target=ids[k];need(records[target]['n']==r['n']-1 and target<i,'strict graded descent')
            children.add(target)
        if children:
            roots[i]=set().union(*(roots[j] for j in children))
            depth[i]=(1+min(depth[j][0] for j in children),1+max(depth[j][1] for j in children))
        else:roots[i]={i};depth[i]=(0,0)
        edges+=len(children)
        out.append({'record':i,'n':r['n'],'rows':rows,'selected':mask,'children':sorted(children),'terminal_birth_records':sorted(roots[i]),'contraction_depth_min_max':depth[i]})
    levels=[]
    for n in range(1,9):
        inds=[i for i,r in enumerate(records) if r['n']==n]
        levels.append({'n':n,'states':len(inds),'births':sum(not out[i]['children'] for i in inds),'nonunique_terminal_births':sum(len(roots[i])>1 for i in inds),'terminal_birth_order_histogram':dict(Counter(records[j]['n'] for i in inds for j in roots[i]))})
    result={'scope':'Every strong marked unit seed with at least one individually stable unselected port through8; exact marked-isomorphism quotient of every eligible stable-twin contraction. Finite confluence check only; no universal confluence inferred.','states':len(records),'edges_between_marked_types':edges,'levels':levels,'records':out}
    tag='.optimized' if not __debug__ else ''
    (HERE/('contraction_dag'+tag+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'states':len(records),'edges':edges,'levels':levels}),flush=True)
if __name__=='__main__':run()
