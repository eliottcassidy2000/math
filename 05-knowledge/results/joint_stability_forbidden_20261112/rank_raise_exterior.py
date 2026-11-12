"""Add one new fixed exterior vertex and one fifth clock to the four-clock core.

Every orientation of the new fixed vertex, every transitive insertion slot
for the new clock, all incidences with the seven fixed vertices, and every
selection on that exterior. This is a complete declared 655,360 graph search.
"""
from itertools import combinations
from collections import Counter
import json,time
import numpy as np
from rank_raise import BASE,PORTS
from census_nine import need,save,tables,strong,ib
from stable_census import selected_defects,certify
from joint_nine import key,quotient,kernel
from cyclic_port_probe import triangles

def run():
    fixed=[i for i in range(10) if i not in PORTS]+[10];order=sorted(PORTS,key=lambda v:-sum(bool(BASE[v]>>w&1) for w in PORTS));P=PORTS+[11];pm=sum(1<<v for v in P)
    stats=Counter();out=[];seen=set();start=time.time();idx=0
    for ext in range(1024):
        prefix=[r|((1<<10) if not ext>>v&1 else 0) for v,r in enumerate(BASE)]+[ext]
        for slot in range(5):
            for sig in range(128):
                idx+=1
                if idx%20000==0:print(json.dumps({'at':idx,'total':655360,'new_types':len(out),'stats':dict(stats),'seconds':round(time.time()-start,1)}),flush=True)
                rows=list(prefix)+[0]
                for i,v in enumerate(fixed):
                    if sig>>i&1:rows[11]|=1<<v
                    else:rows[v]|=1<<11
                for i,v in enumerate(order):
                    if i<slot:rows[v]|=1<<11
                    else:rows[11]|=1<<v
                if not strong(rows):continue
                stats['strong_presentations']+=1;h,u,w,masks=tables(np.array(rows,dtype=np.int64));active=[]
                for m in masks:
                    m=int(m)
                    if m&pm:continue
                    stats['unit_selection_presentations']+=1;_,_,comps=quotient(rows,m)
                    if any(len(set(c)&set(P))>1 for c in comps):continue
                    active.append(m)
                if not active:continue
                stats['distinct_clock_unit_presentations']+=len(active);cert={m:{} for m in active}
                for v in P:
                    R=1+min(rows[v].bit_count(),11-rows[v].bit_count());vals={m:[1] for m in active}
                    for size in range(2,R+1):
                        if not active:break
                        labels,rr=ib.inflate(rows,v,size);sm=[sum(1<<i for i,q in enumerate(labels) if m>>q&1) for m in active]
                        ds=selected_defects(np.array(rr,dtype=np.int64),np.array(sm,dtype=np.int64));stats['axis_tables']+=1;keep=[]
                        for m,d in zip(active,ds):
                            vals[m].append(int(d))
                            if int(d)==1:keep.append(m)
                        active=keep
                    for m in active:cert[m][v]={'R':R,'values':vals[m]}
                    if not active:break
                for m in active:
                    k=key(rows,m)
                    if k in seen:continue
                    seen.add(k);K=kernel(rows,m,P);allp,allcert,_,_=certify(rows,[m]);_,_,comps=quotient(rows,m)
                    rec={'fixed_vertex_signature':ext,'clock_slot':slot,'clock_signature':sig,'rows':rows,'selected':m,'named_ports':P,'all_stable_ports':allp[m],'eligible_blocks':[c for c in comps if len(c)>1 and set(c)&set(allp[m])],'H':int(h),'U':sum(int(u[v]) for v in range(12) if m>>v&1),'axis_certificates':allcert[m],'cyclic_stable_triples':triangles(rows,allp[m]),'joint_named_ports':K=={(1,)*len(P):1},'Stirling_terms':[[*a,c] for a,c in sorted(K.items())]}
                    out.append(rec);print('FIVE_CLOCK_HIT '+json.dumps(rec),flush=True)
    need(idx==655360,'complete declared extension universe')
    result={'scope':'All 655,360 extensions of the inherited four-clock core by one fixed exterior vertex and one new unselected clock. Five named clocks remain in transitive order. All selections on the seven fixed exterior vertices. Full all-length axis certification and joint tensors on hits. This is a template search, not a complete order-twelve census.','graph_presentations':idx,'statistics':dict(stats),'new_marked_types':len(out),'records':out}
    save('rank_raise_exterior.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
