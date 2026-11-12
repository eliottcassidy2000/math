"""Find coupled frozen-length moves preserving an independent stable port.

Discovery screens only. Every positive requires a full multivariate kernel
and an all-length birth audit before it is promoted to a grammar rule.
"""
from common import *
from itertools import combinations
from collections import Counter
import time

def expand(rows,O,lengths):
    labels=[v for v in range(len(rows)) for _ in range(lengths.get(v,1))]
    rr=[sum(1<<j for j,w in enumerate(labels) if (v==w and i<j) or (v!=w and rows[v]>>w&1)) for i,v in enumerate(labels)]
    return rr,sum(1<<i for i,v in enumerate(labels) if O>>v&1),labels

def delta(rows,O,lengths):
    rr,m,_=expand(rows,O,lengths)
    return int(selected_defects(np.array(rr,dtype=np.int64),np.array([m],dtype=np.int64))[0])

def run(limit=10):
    bank=read('birth_union_ten.json.gz',True)['records'];hits=[];stats=Counter();start=time.time()
    for r in bank:
        rows=r['rows'];n=len(rows);O=r['selected'];P=r['ports']
        if n>limit:continue
        for a,b in combinations([v for v in range(n) if not O>>v&1 and v not in P],2):
            stats['pair_screens']+=1
            if delta(rows,O,{a:2,b:2})!=1:continue
            stats['unit_double_screens']+=1
            pp=[p for p in P if delta(rows,O,{a:2,b:2,p:2})==1]
            if not pp:continue
            stats['free_port_screens']+=len(pp)
            if delta(rows,O,{a:3,b:3})!=1:continue
            for p in pp:
                if delta(rows,O,{a:3,b:3,p:2})!=1:continue
                rr,m,_=expand(rows,O,{a:2,b:2});sp,cert,_,_=certify(rr,[m]);_,_,comps=quotient(rr,m)
                eligible=[c for c in comps if len(c)>1 and set(c)&set(sp[m])]
                hits.append({'birth_id':r['birth_id'],'rows':rows,'selected':O,'frozen_ports':[a,b],'free_port':p,'screen_lengths':[1,2,3],'at_2_rows':rr,'at_2_selected':m,'at_2_stable_ports':sp[m],'at_2_eligible_blocks':eligible,'at_2_is_birth':not eligible})
        if r['birth_id']%100==0:print(json.dumps({'at':r['birth_id'],'hits':len(hits),'stats':dict(stats),'seconds':round(time.time()-start,1)}),flush=True)
    out={'scope':f'All stable-clone births through{limit} in the1150 bank, every pair of originally unstable unselected vertices, equal length2/3 screens, each original stable port. Only sampled synchrony; length2 graphs fully certify stable ports and birth status.','statistics':dict(stats),'hits':hits}
    save('synchrony_probe.json.gz',out);print(json.dumps({'hits':len(hits),'birth_hits':sum(x['at_2_is_birth'] for x in hits),'stats':dict(stats)}),flush=True)
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 10)
