"""Lift one-vertex extensions of synchronized births to new diagonal families."""
from common import *
from joint_nine import kernel,expo
from synchrony_probe import expand
from itertools import product
from collections import defaultdict,Counter
import time

def diagonal(E):
    g=defaultdict(int)
    for (a,b),c in E.items():g[a*b]+=c
    return {b:c for b,c in g.items() if c}

def run():
    sources={r['birth_id']:r for r in read('birth_sources_synchrony_12.json.gz')['records']};families=read('synchrony_kernels.json.gz')['records'];out=[];stats=Counter();start=time.time()
    for idx,r in enumerate(read('birth_extensions_13_synchrony.json.gz')['records']):
        src=sources[r['parent_birth_id']];base=families[src['source_family']];clones=[[v for v,q in enumerate(src['labels']) if q==p] for p in base['ports'][:2]]
        for P in product(*clones):
            K=kernel(r['rows'],r['selected'],list(P));E=expo(K);stats['pair_kernels']+=1
            if diagonal(E)!={1:1}:continue
            stats['diagonal_unit']+=1
            if E=={(1,1):1}:stats['independent_growth_only']+=1;continue
            gate=set(E)=={(1,1),(2,1),(1,2)} and E[(1,1)]==1 and E[(2,1)]==-E[(1,2)]
            if not gate:stats['different_diagonal_gate']+=1;continue
            free=[]
            for p in r['stable_ports']:
                if p in P:continue
                pp=list(P)+[p]
                if sorted(sum(r['rows'][v]>>w&1 for w in pp) for v in pp)!=[0,1,2]:continue
                triple=expo(kernel(r['rows'],r['selected'],pp));stats['triple_kernels']+=1
                td=defaultdict(int)
                for (a,b,z),c in triple.items():td[a*b,z]+=c
                if {b:c for b,c in td.items() if c}=={(1,1):1}:free.append({'port':p,'triple_exponential':[[*b,c] for b,c in sorted(triple.items())]})
            if not free:continue
            _,_,C=quotient(r['rows'],r['selected']);frozen=[]
            for block in C:
                if len(block)<2:continue
                if set(block)&set(P):continue
                labs,rr=ib.inflate(r['rows'],block[0],2);O=sum(1<<v for v,q in enumerate(labs) if r['selected']>>q&1);pp=[labs.index(p) for p in P]
                g=diagonal(expo(kernel(rr,O,pp)));g[1]=g.get(1,0)-1
                frozen.append({'block':block,'increment_defect_minus_one':[[b,c] for b,c in sorted(g.items()) if c]})
            out.append({'parent_extension_index':idx,'parent_synchrony_family':src['source_family'],'rows':r['rows'],'selected':r['selected'],'coupled_ports':list(P),'free_ports':free,'defect_exponential':[[*b,c] for b,c in sorted(E.items())],'frozen_tests':frozen})
        if idx%10==0:
            save('mixed_family_probe.json.gz',{'scope':'Discovery/full signed-kernel certificates for pairwise equal growth in88 one-vertex children of synchronized births. Birth nonvanishing and family intersections need separate proof.','complete':False,'statistics':dict(stats),'records':out});print(json.dumps({'at':idx,'families':len(out),'stats':dict(stats),'seconds':round(time.time()-start,1)}),flush=True)
    save('mixed_family_probe.json.gz',{'scope':'Complete declared352 pair-growth tests in88 one-vertex children of synchronized births, plus full free-port kernels and old-block increment tests. Birth nonvanishing and all-length intersections require the separate proof.','complete':True,'statistics':dict(stats),'records':out});print(json.dumps({'families':len(out),'stats':dict(stats),'seconds':round(time.time()-start,1)}),flush=True)
if __name__=='__main__':run()
