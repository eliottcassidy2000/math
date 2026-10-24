#!/usr/bin/env python3
"""Regular-nine classes from 85 pointed states and parent-only reroot trace.
Does not import preparation, count target regular classes, or canonicalize
any child. Every canonical call is runtime guarded to order eight.
"""
import argparse, hashlib, itertools, json
from pathlib import Path
from collections import Counter
CHECKS=Counter();OPS=Counter()
def require(c,m):
    CHECKS[m]+=1
    if not c: raise RuntimeError(m)
def stable_colors(rows):
    n=len(rows); c=tuple(r.bit_count() for r in rows)
    while True:
        masks=[sum(1<<v for v in range(n) if c[v]==x) for x in sorted(set(c))]
        sig=[(c[v],tuple((rows[v]&m).bit_count() for m in masks)) for v in range(n)]
        names={s:i for i,s in enumerate(sorted(set(sig)))};new=tuple(names[s] for s in sig)
        if len(set(new))==len(set(c)):return new
        c=new
def canonical(rows):
    require(len(rows)==8,'canonicalization strictly on order-eight parents');OPS['parent_canonicalizations']+=1
    c=stable_colors(rows);cells=[tuple(v for v in range(8) if c[v]==x) for x in sorted(set(c))]
    best=None
    for pieces in itertools.product(*(itertools.permutations(cell) for cell in cells)):
        p=sum(pieces,());value=tuple(sum(((rows[v]>>w)&1)<<j for j,w in enumerate(p)) for v in p)
        OPS['parent_canonical_label_candidates']+=1
        if best is None or value<best:best=value
    return best
def extension(rows):
    require(set(r.bit_count() for r in rows)=={3,4},'almost-regular parent scores')
    child=list(rows)+[0]
    for v in range(8):
        if rows[v].bit_count()==3:child[v]|=1<<8
        else:child[8]|=1<<v
    require(all(r.bit_count()==4 for r in child),'unique extension is regular nine')
    OPS['regular_child_extensions']+=1
    return tuple(child)
def deletion(rows,v):
    order=[w for w in range(9) if w!=v];OPS['regular_child_deletions']+=1
    return tuple(sum(((rows[a]>>b)&1)<<j for j,b in enumerate(order)) for a in order)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='regular-normal.json');args=parser.parse_args()
    here=Path(__file__).resolve().parent;raw=(here/'parents-eight.json').read_bytes()
    require(hashlib.sha256(raw).hexdigest()=='1cd97666c0ee6294bf379b6d77f95bea5e24f8f5b27749ef804c46a45892ac7b','pinned complete parent input')
    payload=json.loads(raw);parents=[tuple(r['rows']) for r in payload['parents'] if set(v.bit_count() for v in r['rows'])=={3,4}]
    index={rows:i for i,rows in enumerate(parents)};D=[];trace=0
    for i,rows in enumerate(parents):
        child=extension(rows);destinations=Counter()
        for v in range(9):
            key=canonical(deletion(child,v));require(key in index,'all reroot parents remain in complete almost-regular atlas');destinations[index[key]]+=1
        require(sum(destinations.values())==9,'regular reroot row sum is nine');require(destinations[i]>0,'regular reroot has self root')
        D.append(destinations);trace+=destinations[i]
    require(trace%9==0,'regular trace divisible by nine')
    remaining=set(range(len(parents)));components=[]
    while remaining:
        seed=min(remaining);component=set(D[seed]);require(seed in component,'own state in regular child component')
        require(component<=remaining,'regular child components disjoint')
        require(all(D[j]==D[seed] for j in component),'regular reroot rows identical within child component')
        remaining-=component
        diag_sum=sum(D[j][j] for j in component);require(diag_sum==9,'each regular child trace is nine')
        orbit_sizes=sorted(D[seed].values());error=sum(9-s for s in orbit_sizes)
        require(error==9*(len(component)-1),'regular error formula per component')
        components.append({'states':sorted(component),'root_orbit_sizes':orbit_sizes,'error':error})
    require(len(components)==trace//9,'independent components agree with regular trace')
    error=9*(len(parents)-trace//9);require(sum(c['error'] for c in components)==error,'sum regular component errors')
    result={'status':'PASS','pointed_regular_states':len(parents),'regular_trace':trace,'regular_classes':trace//9,'regular_error':error,'root_orbit_size_multisets':sorted(Counter(tuple(c['root_orbit_sizes']) for c in components).items()),'components':components,'states':[{'parent':list(rows),'reroot_row':sorted(D[i].items())} for i,rows in enumerate(parents)],'operations':dict(OPS),'checks':dict(CHECKS),'child_canonicalizations':0}
    (here/args.output).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps({k:result[k] for k in ('status','pointed_regular_states','regular_trace','regular_classes','regular_error','operations')},sort_keys=True))
if __name__=='__main__':main()
