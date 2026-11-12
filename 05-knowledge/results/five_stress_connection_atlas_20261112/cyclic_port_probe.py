"""Can an individually stable set of unselected ports contain a directed cycle?

Complete inherited ordinary class universe through8, all unit selections.
Screen every port by literal two-clone defect; only candidate cyclic triples
need full all-length coefficient tests. A tournament is acyclic iff it has
no cyclic triple. This is a finite-order statement about arbitrary cores.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
import gzip,json,hashlib,sys,time
import numpy as np
from extend_round import tables
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'endpoint_coverage_compiler_20261111'))
import inflation_bank as ib
def need(x,msg):
    if not x:raise RuntimeError(msg)
def triangles(rows,verts):
    return [list(tri) for tri in combinations(verts,3) if all(sum(bool(rows[v]>>w&1) for w in tri if v!=w)==1 for v in tri)]
def run():
    source=HERE.parent/'metagraph_aspect_sequences_20261111/records.jsonl.gz'
    with gzip.open(source,'rt') as f:bank=[json.loads(line) for line in f]
    out=[];controls=0;cert=[];t=time.time()
    for n in range(1,9):
        corecount=0;marks=0;ports_checked=0;candidates=[];survivors=[];hist=Counter();stablehist=Counter();false_screen=0;births=0;descendants=0
        for item in bank:
            if item['n']!=n or not item['strong'] or not item['seeds']:continue
            corecount+=1;rows=item['rows'];H,U,W,unit=tables(np.array(rows,dtype=np.int64))
            need(int(H)==item['H'],'input H')
            need(set(unit)=={q[0] for q in item['seeds']},'complete unit selections')
            selections=sorted({item['attachment_orbit'][q[0]] for q in item['seeds']});marks+=len(selections)
            eligible={v for v in range(n) if any(not m>>v&1 for m in selections)}
            children={}
            for v in eligible:
                labels,child=ib.inflate(rows,v,2)
                h,u,w,_=tables(np.array(child,dtype=np.int64));children[v]=(labels,int(h),w);ports_checked+=1
            exact={}
            for mask in selections:
                O=[v for v in range(n) if mask>>v&1];passing=[]
                for v in eligible:
                    if v in O:continue
                    labels,h,w=children[v];sel=[i for i,q in enumerate(labels) if q in O]
                    if h-sum(int(w[a,b]) for a,b in combinations(sel,2))==1:passing.append(v)
                hist[len(passing)]+=1
                tri=triangles(rows,passing)
                if tri:candidates.append({'rep':item['rep'],'selected':mask,'ports':passing,'triangles':tri})
                if not passing:stablehist[0]+=1;continue
                stable=[];coefs={}
                for v in passing:
                    if v not in exact:exact[v]=ib.coefficients(rows,v)
                    h,u,w=exact[v];delta=[h[j]-sum(w[j][a][b] for a,b in combinations(O,2)) for j in range(1,len(h))]
                    coefs[v]=delta
                    if delta==[1]+[0]*(n-1):stable.append(v)
                    else:false_screen+=1
                stablehist[len(stable)]+=1
                genuine=triangles(rows,stable)
                contractions=[]
                for a,b in combinations([v for v in range(n) if v not in O],2):
                    if not (a in stable or b in stable):continue
                    if not all(bool(rows[a]>>v&1)==bool(rows[b]>>v&1) for v in range(n) if v not in [a,b]):continue
                    kept=[v for v in range(n) if v!=a]
                    child=[sum(1<<j for j,v in enumerate(kept) if rows[u]>>v&1) for u in kept]
                    sel=[kept.index(v) for v in O];port=kept.index(b)
                    h,u,w=ib.coefficients(child,port)
                    dd=[h[j]-sum(w[j][i][k] for i,k in combinations(sel,2)) for j in range(1,len(h))]
                    need(dd==[1]+[0]*(n-2),'stable twin contraction')
                    contractions.append([a,b])
                if stable:
                    if contractions:descendants+=1
                    else:births+=1
                cert.append({'n':n,'rep':item['rep'],'rows':rows,'selected':mask,'defect_coefficients':coefs,'individually_stable':stable,'cyclic_triples':genuine,'stable_twin_contractions':contractions})
                if genuine:survivors.append(cert[-1])
        level={'n':n,'strong_classes_with_unit_selections':corecount,'marked_selection_orbits':marks,'literal_port_clones':ports_checked,'two_clone_preserving_port_histogram':dict(hist),'all_length_stable_port_histogram':dict(stablehist),'stable_clone_births':births,'stable_clone_descendants':descendants,'false_positive_port_screens':false_screen,'cyclic_screen_candidates':len(candidates),'all_length_cyclic_survivors':len(survivors)}
        out.append(level);print(json.dumps({**level,'seconds':round(time.time()-t,2)}),flush=True)
    # Positive controls: every earlier stable strong-core port passes literally.
    old=json.loads((HERE.parent/'endpoint_coverage_compiler_20261111/inflation_bank_7.json').read_text())['records']
    for r in old:
        rows=r['rows'];O=[v for v in range(len(rows)) if r['selected']>>v&1]
        labels,child=ib.inflate(rows,r['inflated_vertex'],2);h,u,w,_=tables(np.array(child,dtype=np.int64))
        sel=[i for i,q in enumerate(labels) if q in O]
        need(int(h)-sum(int(w[a,b]) for a,b in combinations(sel,2))==1,'inherited positive control');controls+=1
    result={'scope':'Every inherited ordinary tournament class through8, restricted to strong classes with a unit selection; all selections modulo actual automorphism action. One literal doubling per eligible port, followed by full exact coefficients for every port passing that screen. Individual stability does not assert joint stability.','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'levels':out,'positive_screen_certificates':cert,'old_stable_port_controls':controls}
    tag='.optimized' if not __debug__ else ''
    (HERE/('cyclic_port_receipt'+tag+'.json')).write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':run()
