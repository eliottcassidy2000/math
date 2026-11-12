"""Complete mixed stable-face complexes for the order-ten stable bank."""
from collections import Counter
from itertools import combinations
import time,json
import numpy as np
from census_nine import read,save,need,tables,defect,ib
from joint_nine import compress,kernel,expo,evaluate,key,cones

def faces_at(E,lengths):
    t=len(lengths);faces=[]
    for mask in range(1<<t):
        p=Counter()
        for bases,c in E.items():
            active=[]
            for i,b in enumerate(bases):
                if mask>>i&1:active.append(b)
                else:c*=b**lengths[i]
            p[tuple(active)]+=c
        p={a:c for a,c in p.items() if c}
        if p=={(1,)*mask.bit_count():1}:faces.append(mask)
    facets=[m for m in faces if not any(m!=q and m&q==m for q in faces)]
    minimal_nonfaces=[m for m in range(1<<t) if m not in faces and all(m^(1<<i) in faces for i in range(t) if m>>i&1)]
    return faces,facets,minimal_nonfaces

def run():
    bank=read('stable_ten_selected.json.gz')['records'];out=[];hist=Counter();birthhist=Counter();fail=[];start=time.time()
    for idx,r in enumerate(bank):
        rr,mm,P,blocks=compress(r['rows'],r['selected'],r['stable_ports']);lengths=[len(c) for c in blocks]
        K=kernel(rr,mm,P);E=expo(K);ok=K=={(1,)*len(P):1};faces,facets,nonfaces=faces_at(E,lengths)
        need(0 in faces and all(1<<i in faces for i in range(len(P))),'every certified original axis is a face')
        need(ok==((1<<len(P))-1 in faces),'full face identity on positive tails')
        if not ok:fail.append(idx)
        hist[len(P),ok]+=1;birth=not r['eligible_blocks']
        if birth:birthhist[len(P),tuple(sorted(m.bit_count() for m in facets))]+=1
        # An uneven mixed literal control, using the actual block base lengths.
        sizes=[x+1+(i%2) for i,x in enumerate(lengths)]
        expanded,O,_=cones.expand(rr,[i for i in range(len(rr)) if mm>>i&1],P,sizes)
        actual=None
        if len(expanded)<=18:
            h,u,w,_=tables(np.array(expanded,dtype=np.int64));actual=defect(h,w,O);need(actual==evaluate(E,sizes),'independent mixed literal control')
        out.append({'code':r['code'],'selected':r['selected'],'birth':birth,'stable_blocks':blocks,'reduced_rows':rr,'reduced_selected':mm,'clock_ports':P,'base_lengths':lengths,'joint_stable':ok,'stable_faces':faces,'maximal_faces':facets,'minimal_nonfaces':nonfaces,'Stirling_terms':[[*a,c] for a,c in sorted(K.items())],'exponential_terms':[[*a,c] for a,c in sorted(E.items())],'literal_lengths':sizes,'literal_defect':actual})
        if idx%100==0:print(json.dumps({'at':idx,'total':len(bank),'mixed_failures':len(fail),'seconds':round(time.time()-start,1)}),flush=True)
    parents=read('stable_nine.json.gz')['records'];expected=set();presentations=0
    for r in parents:
        for v in r['stable_ports']:
            labels,rows=ib.inflate(r['rows'],v,2);mask=sum(1<<i for i,q in enumerate(labels) if r['selected']>>q&1);expected.add(key(rows,mask));presentations+=1
    actualkeys={key(r['rows'],r['selected']) for r in bank if r['eligible_blocks']}
    need(expected==actualkeys and len(actualkeys)==403,'complete descendant key equality')
    result={'scope':'Every complete strong marked order-ten unit seed with a stable port, all-length signed mixed tensors after justified stable-block compression; stable-face complexes anchored at actual block lengths. All stable sets acyclic.','states':len(bank),'joint_failures':fail,'clock_rank_joint_histogram':[[a,b,c] for (a,b),c in sorted(hist.items())],'birth_rank_facet_histogram':[[a,list(b),c] for (a,b),c in sorted(birthhist.items())],'descendant_key_equality':True,'one_clone_presentations_from_order9':presentations,'records':out}
    save('joint_ten.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
