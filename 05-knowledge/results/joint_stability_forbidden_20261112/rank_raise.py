"""New-port births allowing selection changes and paired insertions.

The one-port universe uses every orientation, including cyclic port sets.
The two-port universe preserves the transitive order of all named clocks.
Both search every selection on the fixed exterior, not just the old mark.
"""
from pathlib import Path
from itertools import combinations,product
from collections import Counter
import sys,time,json
import numpy as np
from census_nine import need,save,tables,ib,strong
from stable_census import selected_defects
from joint_nine import kernel,key,quotient
from cyclic_port_probe import triangles
BASE=[670,668,8,528,36,847,671,60,735,148]
PORTS=[1,2,8,9]

def graphs(extra):
    n=len(BASE)
    if extra==1:
        for sig in range(1<<n):
            yield {'signature':sig},[r|((1<<n) if not sig>>v&1 else 0) for v,r in enumerate(BASE)]+[sig]
        return
    fixed=[i for i in range(n) if i not in PORTS]
    old=sorted(PORTS,key=lambda v:-sum(bool(BASE[v]>>w&1) for w in PORTS))
    need(all(BASE[old[i]]>>old[j]&1 for i in range(4) for j in range(i+1,4)),'old clock order')
    for slots in combinations(range(6),2):
        order=[];oi=0;ni=0
        for k in range(6):
            if k in slots:order.append(n+ni);ni+=1
            else:order.append(old[oi]);oi+=1
        for sigs in product(range(64),repeat=2):
            rr=list(BASE)+[0,0]
            for j,sig in enumerate(sigs):
                v=n+j
                for i,w in enumerate(fixed):
                    if sig>>i&1:rr[v]|=1<<w
                    else:rr[w]|=1<<v
            for i,v in enumerate(order):
                for w in order[i+1:]:
                    if v>=n or w>=n:rr[v]|=1<<w
            yield {'slots':slots,'signatures':sigs},rr

def run(extra):
    need(extra in [1,2],'declared insertion universe');P=PORTS+list(range(10,10+extra));stats=Counter();out=[];seen=set();start=time.time()
    for idx,(construction,rows) in enumerate(graphs(extra),1):
        if not strong(rows):continue
        stats['strong_graph_presentations']+=1
        h,u,w,masks=tables(np.array(rows,dtype=np.int64));pm=sum(1<<v for v in P)
        active=[int(m) for m in masks if not m&pm]
        if not active:continue
        stats['unit_selection_presentations']+=len(active)
        # Reject repeated clocks, even when they look like different vertices.
        keep=[]
        for m in active:
            _,_,comps=quotient(rows,m)
            if any(len(set(c)&set(P))>1 for c in comps):continue
            keep.append(m)
        active=keep
        if not active:continue
        stats['distinct_clock_unit_presentations']+=len(active);cert={m:{} for m in active}
        for v in P:
            R=1+min(rows[v].bit_count(),len(rows)-1-rows[v].bit_count());vals={m:[1] for m in active}
            for size in range(2,R+1):
                if not active:break
                labels,rr=ib.inflate(rows,v,size);sm=[sum(1<<i for i,q in enumerate(labels) if m>>q&1) for m in active]
                ds=selected_defects(np.array(rr,dtype=np.int64),np.array(sm,dtype=np.int64));stats['axis_tables']+=1;remaining=[]
                for m,d in zip(active,ds):
                    vals[m].append(int(d))
                    if int(d)==1:remaining.append(m)
                active=remaining
            for m in active:cert[m][v]={'R':R,'values':vals[m]}
            if not active:break
        for m in active:
            k=key(rows,m)
            if k in seen:continue
            seen.add(k);cyclic=triangles(rows,P);K=None if cyclic else kernel(rows,m,P)
            rec={'construction':construction,'rows':rows,'selected':m,'ports':P,'H':int(h),'U':sum(int(u[v]) for v in range(len(rows)) if m>>v&1),'axis_certificates':cert[m],'cyclic_triples':cyclic,'joint_stable':None if K is None else K=={(1,)*len(P):1},'Stirling_terms':None if K is None else [[*a,c] for a,c in sorted(K.items())]}
            out.append(rec);print('NEW_CLOCK_FAMILY '+json.dumps(rec),flush=True)
        if idx%5000==0:print(json.dumps({'at':idx,'seconds':round(time.time()-start,1),'stats':dict(stats),'new_types':len(out)}),flush=True)
    result={'scope':f'All {idx} declared {extra}-vertex extensions of the inherited Session288 four-clock core; all selections on the six fixed exterior vertices. For two new vertices the whole named port set is transitive. Stable axes certified through the exact endpoint run bound. Different vertices in the same twin block are not counted as new clocks. Not a complete tournament universe.','base_rows':BASE,'base_ports':PORTS,'graph_presentations':idx,'statistics':dict(stats),'new_marked_types':len(out),'records':out}
    save(f'rank_raise_{extra}.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 1)
