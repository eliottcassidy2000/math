"""Complete order-eleven stable-clone descendants, without enumerating births.

Every eligible contraction at eleven lands in the complete order-ten stable
bank. Reversing that one contraction therefore enumerates this universe.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
import json,time
import numpy as np
from census_nine import read,save,need,tables,defect,ib
from stable_census import certify
from joint_nine import key,quotient,compress,kernel,expo
from cyclic_port_probe import triangles
from joint_ten import faces_at
HERE=Path(__file__).resolve().parent

def remove(rows,mask,v):
    keep=[i for i in range(len(rows)) if i!=v]
    return [sum(1<<j for j,w in enumerate(keep) if rows[u]>>w&1) for u in keep],sum(1<<j for j,w in enumerate(keep) if mask>>w&1)

def run():
    old=json.loads((HERE.parent/'five_stress_connection_atlas_20261112/cyclic_port_receipt.json').read_text())['positive_screen_certificates']
    banks=[{'rows':r['rows'],'selected':r['selected'],'stable_ports':r['individually_stable']} for r in old if r['individually_stable']]
    banks+=read('stable_nine.json.gz')['records']+read('stable_ten_selected.json.gz')['records']
    lower={key(r['rows'],r['selected']):r for r in banks};roots={};edges={}
    for k,r in sorted(lower.items()):
        rows=r['rows'];m=r['selected'];_,_,comps=quotient(rows,m);kids=set()
        for c in comps:
            if len(c)>1 and set(c)&set(r['stable_ports']):
                rr,mm=remove(rows,m,c[0]);kk=key(rr,mm);need(kk in roots,'complete smaller contraction bank');kids.add(kk)
        edges[k]=kids;roots[k]=set().union(*(roots[j] for j in kids)) if kids else {k}
        need(len(roots[k])==1,'unique roots through ten')
    candidates={};presented=0
    for k,r in lower.items():
        if len(r['rows'])!=10:continue
        for v in r['stable_ports']:
            labels,rr=ib.inflate(r['rows'],v,2);mm=sum(1<<i for i,q in enumerate(labels) if r['selected']>>q&1);kk=key(rr,mm);presented+=1
            candidates.setdefault(kk,(rr,mm))
    union=read('birth_union_ten.json.gz');need(len(candidates)==union['counts_n1_to40'][10],'predicted initial-family coefficient')
    out=[];hist=Counter();jointfail=[];start=time.time()
    for i,(k,(rows,m)) in enumerate(sorted(candidates.items())):
        h,u,w,_=tables(np.array(rows,dtype=np.int64));need(defect(h,w,[v for v in range(11) if m>>v&1])==1,'literal generated unit seed')
        pp,cert,false,calls=certify(rows,[m]);P=pp[m];_,_,comps=quotient(rows,m);eligible=[c for c in comps if len(c)>1 and set(c)&set(P)];kids=set()
        for c in eligible:
            rr,mm=remove(rows,m,c[0]);kk=key(rr,mm);need(kk in lower,'every reverse edge lands in complete bank');kids.add(kk)
        terminal=set().union(*(roots[j] for j in kids));hist[len(terminal)]+=1
        rr,mm,ps,blocks=compress(rows,m,P);cyclic=triangles(rr,ps);K=None if cyclic else kernel(rr,mm,ps)
        joint=K=={(1,)*len(ps):1} if K is not None else None
        if joint is not True:jointfail.append(i)
        out.append({'rows':rows,'selected':m,'H':int(h),'U':sum(int(u[v]) for v in range(11) if m>>v&1),'stable_ports':P,'eligible_blocks':eligible,'axis_certificates':cert[m],'children_keys':sorted(kids),'terminal_birth_keys':sorted(terminal),'joint_stable':joint,'cyclic_stable_triples':cyclic,'Stirling_terms':None if K is None else [[*a,c] for a,c in sorted(K.items())]})
        if i%100==0:print(json.dumps({'at':i,'total':len(candidates),'root_histogram':dict(hist),'joint_failures':len(jointfail),'seconds':round(time.time()-start,1)}),flush=True)
    result={'scope':'Complete ordinary strong marked order-eleven unit seeds with at least one eligible stable-twin contraction. All generated from the complete order-ten stable bank; order-eleven births excluded. Every stable port and contraction recertified.','through_ten_states':len(lower),'through_ten_contraction_edges':sum(len(v) for v in edges.values()),'one_clone_presentations':presented,'marked_descendants':len(out),'edges_to_order_ten':sum(len(r['children_keys']) for r in out),'terminal_birth_count_histogram':dict(hist),'nonjoint_records':jointfail,'records':out}
    save('descendants_eleven.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
