"""Complete coprime common-source atlas below target 125, using inherited P1."""
from collections import defaultdict,Counter
from functools import lru_cache
from itertools import permutations,product,combinations
from math import gcd
from pathlib import Path
import csv,gzip,json
from enumerate import scc,hcount,child,require

HERE=Path(__file__).resolve().parent

def pareto(pairs):
    return sorted({p for p in pairs if not any(q!=p and q[0]<=p[0] and q[1]<=p[1] for q in pairs)})

@lru_cache(None)
def normalize(rows):
    comps=scc(rows)
    a=0
    while a<len(comps) and comps[a].bit_count()==1:
        a+=1
    require(a<len(comps),'nontrivial source')
    last=len(comps)-1
    while comps[last].bit_count()==1:
        last-=1
    b=len(comps)-last-1
    mask=sum(comps[a:last+1])
    vs=[i for i in range(len(rows)) if mask>>i&1]
    degrees={i:(rows[i]&mask).bit_count() for i in vs}
    groups=[[i for i in vs if degrees[i]==d] for d in sorted(set(degrees.values()))]
    best=None;best_order=None
    for pp in product(*(permutations(group) for group in groups)):
        order=sum((tuple(p) for p in pp),())
        code=0
        for i in range(len(order)):
            for j in range(i+1,len(order)):
                code=2*code+((rows[order[i]]>>order[j])&1)
        if best is None or code<best:
            best=code;best_order=order
    prefix=[c.bit_length()-1 for c in comps[:a]]
    suffix=[c.bit_length()-1 for c in comps[last+1:]]
    order=prefix+list(best_order)+suffix
    core=[sum(1<<j for j,v in enumerate(best_order) if rows[u]>>v&1) for u in best_order]
    return (len(vs),best),a,b,order,core

def join_source(core,a,b):
    n=a+len(core)+b
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    for i,r in enumerate(core):
        rows[a+i]=(r<<a)|sum(1<<j for j in range(a+len(core),n))
    return rows

def run():
    atlas=defaultdict(dict)
    tested=0
    with (HERE/'primitive_births.tsv').open(encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            h,g=int(r['source_H']),int(r['target_H'])
            if h<=1 or gcd(h,g)!=1:
                continue
            rows=tuple(map(int,r['source_rows'].split(',')))
            key,a,b,order,core=normalize(rows)
            old_out=int(r['new_vertex_out_mask'])
            out=sum(1<<i for i,v in enumerate(order) if old_out>>v&1)
            entry=atlas[h].setdefault(key,dict(core_rows=core,targets={}))
            entry['targets'].setdefault(g,{}).setdefault((a,b),out)
    exported={}
    summary={}
    for h,cores in sorted(atlas.items()):
        facets=set()
        entries=[]
        for key,entry in sorted(cores.items()):
            targets={}
            for g,pads in sorted(entry['targets'].items()):
                minimal=pareto(set(pads))
                targets[g]=[dict(a=a,b=b,out_mask=pads[a,b]) for a,b in minimal]
                # Every retained normalized certificate is checked literally.
                for a,b in minimal:
                    rows=join_source(entry['core_rows'],a,b)
                    require(hcount(rows)==h and hcount(child(rows,pads[a,b]))==g,'canonical padding witness')
                    tested+=1
            entries.append(dict(core_key=list(key),core_rows=entry['core_rows'],targets=targets))
            facets.add(frozenset(targets))
        maximal=sorted([sorted(f) for f in facets if not any(f<g for g in facets)])
        vertices=sorted(set().union(*facets))
        incompatible=[list(p) for p in combinations(vertices,2) if not any(set(p)<=f for f in facets)]
        exported[str(h)]=dict(cores=entries,maximal_common_target_sets=maximal,incompatible_target_pairs=incompatible)
        summary[str(h)]=dict(cores=len(cores),targets=len(vertices),maximal_facets=len(maximal),
                             incompatible_pairs=len(incompatible))
    ternary={}
    for h in range(3,40,2):
        g1,g2=3*h+2,3*h+4
        numerator=Counter()
        compatible=[]
        for entry in exported.get(str(h),{}).get('cores',[]):
            targets=entry['targets']
            if g1 not in targets or g2 not in targets:
                continue
            minima=pareto({(max(p['a'],q['a']),max(p['b'],q['b'])) for p in targets[g1] for q in targets[g2]})
            m=len(entry['core_rows'])
            for i,(a,b) in enumerate(minima):
                numerator[m+a+b]+=1
                if i+1<len(minima):
                    numerator[m+minima[i+1][0]+b]-=1
            # One common physical source for both branches, checked directly.
            a,b=minima[0]
            rows=join_source(entry['core_rows'],a,b)
            for g in [g1,g2]:
                p=next(p for p in targets[g] if p['a']<=a and p['b']<=b)
                shift=a-p['a']
                out=(p['out_mask']<<shift)|sum(1<<j for j in range(a+m+p['b'],a+m+b))
                require(hcount(child(rows,out))==g,'literal common-source padding')
                tested+=1
            compatible.append(dict(core_key=entry['core_key'],minimal_padding_pairs=minima))
        ternary[str(h)]=dict(targets=[g1,g2],common_cores=compatible,
            source_class_GF_numerator={str(k):v for k,v in sorted(numerator.items()) if v},
            denominator='(1-x)^2',eventual_slope=len(compatible))
    (HERE/'coprime_source_compatibility_below125.json').write_text(json.dumps(exported,indent=2)+'\n',encoding='utf-8')
    (HERE/'ternary_common_source_families_through39.json').write_text(json.dumps(ternary,indent=2)+'\n',encoding='utf-8')
    observed_checks=0
    observed33=Counter()
    prev=HERE.parent/'metagraph_aspect_sequences_20261111/records.jsonl.gz'
    with gzip.open(prev,'rt',encoding='utf-8') as f:
        for line in f:
            r=json.loads(line)
            h=r['H']
            if not 1<h<125:
                continue
            key,a,b,order,core=normalize(tuple(r['rows']))
            targets=atlas.get(h,{}).get(key,{}).get('targets',{})
            predicted={g for g,pairs in targets.items() if any(p<=a and q<=b for p,q in pairs)}
            actual={g for g in r['extension'] if g<125 and gcd(h,g)==1}
            require(predicted==actual,'independent complete source availability')
            observed_checks+=1
            if h==33 and {101,103}<=actual:
                observed33[r['n']]+=1
    require(dict(observed33)=={7:2,8:5},'H33 simultaneous-source counts')
    result=dict(status='PASS',nontrivial_source_counts=len(atlas),source_core_types=sum(x['cores'] for x in summary.values()),
        literal_normalized_or_common_witness_checks=tested,by_source_H=summary,
        independent_source_availability_checks=observed_checks,H33_small_common_counts=dict(observed33),
        h33_common=ternary['33'])
    (HERE/'compatibility_checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['by_source_H','h33_common']}))
    print('H33 common source numerator:',ternary['33']['source_class_GF_numerator'])

if __name__=='__main__':
    run()
