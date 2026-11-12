"""Full kernels for the coupled frozen blocks and independent free port."""
from common import *
from joint_nine import cones,kernel,expo,evaluate
from synchrony_probe import expand
from collections import defaultdict,Counter
from itertools import product
import time

def diagonal(E):
    out=defaultdict(int)
    for (a,b,z),c in E.items():out[a*b,z]+=c
    return {q:c for q,c in out.items() if c}

def run():
    out=[];start=time.time();keys=set();stats=Counter()
    for i,r in enumerate(read('synchrony_probe.json.gz')['hits']):
        rows=r['rows'];O=r['selected'];P=r['frozen_ports']+[r['free_port']]
        need(sorted(sum(rows[v]>>w&1 for w in P) for v in P)==[0,1,2],'acyclic triple')
        K=kernel(rows,O,P);D,U=cones.kernel(rows,P,[v for v in range(len(rows)) if O>>v&1]);need(K==D,'independent complete defect kernels')
        E=expo(K);UE=expo(U);DE=diagonal(E);need(DE=={(1,1):1},'all-length diagonal unit/free-port identity')
        rr,m,_=expand(rows,O,{P[0]:2,P[1]:2});keys.add(key(rr,m))
        q,om,comps=quotient(rows,O);frozen=[]
        for block in comps:
            if len(block)<2:continue
            need(not set(block)&set(P),'existing frozen blocks distinct from growing vertices')
            lab,child=ib.inflate(rows,block[0],2);cm=sum(1<<j for j,v in enumerate(lab) if O>>v&1);pp=[lab.index(p) for p in P[:2]]
            X=expo(kernel(child,cm,pp));g=defaultdict(int)
            for (a,b),c in X.items():g[a*b]+=c
            g[1]-=1;g={b:c for b,c in g.items() if c}
            frozen.append({'base_block':block,'increment_defect_minus_one':[[b,c] for b,c in sorted(g.items())]})
        grid=[]
        for a,b,z in product(range(1,4),range(1,4),range(1,3)):
            rr,m,_=expand(rows,O,dict(zip(P,[a,b,z])));h,u,w,_=tables(np.array(rr,dtype=np.int64));ss=[j for j in range(len(rr)) if m>>j&1]
            d=defect(h,w,ss);s=sum(int(u[j]) for j in ss);need(d==evaluate(E,[a,b,z]) and s==evaluate(UE,[a,b,z]),'literal complete observable control');grid.append([a,b,z,int(h),d,s])
        out.append({**r,'ports':P,'defect_Stirling':[[*k,c] for k,c in sorted(K.items())],'defect_exponential':[[*k,c] for k,c in sorted(E.items())],'U_exponential':[[*k,c] for k,c in sorted(UE.items())],'diagonal_U':[[*k,c] for k,c in sorted(diagonal(UE).items())],'existing_frozen_block_tests':frozen,'literal_grid':grid})
        stats[tuple(sorted(E.values()))]+=1
        print(json.dumps({'at':i+1,'total':24,'seconds':round(time.time()-start,1),'defect_terms':len(E),'frozen_tests':frozen}),flush=True)
    result={'scope':'All24 synchrony candidates; independent complete acyclic triple kernels, exact diagonal identities, literal grids and tests of existing frozen blocks. Birth status at all lengths additionally requires the sign/nonvanishing proof recorded with these tests.','families':len(out),'distinct_births_at_order12':len(keys),'records':out}
    save('synchrony_kernels.json.gz',result);print(json.dumps({'families':len(out),'distinct_births_at_order12':len(keys)}),flush=True)
if __name__=='__main__':run()
