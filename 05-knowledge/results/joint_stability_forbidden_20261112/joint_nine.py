"""Full mixed tensors for the complete order-nine stable bank.

Stable twin blocks are compressed solely as polynomial identity tests.
The uncompressed family occupies an infinite product tail, so exponential
independence makes the compressed and original identity questions equivalent.
No contraction stability assumption enters this reduction.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
import sys,json,gzip,time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'five_stress_connection_atlas_20261112'))
from birth_cones import quotient
from contraction_dag import key
from joint_port_extension import cones
from extend_round import tables
sys.path.insert(0,str(HERE.parent/'session288_twin_lengths_20261112/stable'))
from tensor import kernel,expo,evaluate
from census_nine import need,save,read,defect

def compress(rows,mask,ports):
    rr,mm,components=quotient(rows,mask)
    kept=[];blocks=[]
    for c in components:
        hit=set(c)&set(ports)
        if hit:
            need(hit==set(c),'stability is uniform in a twin block')
            kept.append(c[0]);blocks.append(c)
        else:kept.extend(c)
    kept.sort()
    reduced=[sum(1<<j for j,v in enumerate(kept) if rows[u]>>v&1) for u in kept]
    selected=sum(1<<j for j,v in enumerate(kept) if mask>>v&1)
    reps=[kept.index(c[0]) for c in blocks]
    return reduced,selected,reps,blocks

def run():
    bank=read('stable_nine.json.gz')['records'];out=[];start=time.time();hist=Counter();fail=[]
    for idx,r in enumerate(bank):
        rr,mm,P,blocks=compress(r['rows'],r['selected'],r['stable_ports'])
        K=kernel(rr,mm,P);E=expo(K);ok=K=={(1,)*len(P):1}
        if not ok:fail.append(idx)
        hist[(len(P),ok)]+=1
        # Literal nonuniform multi-inflation controls exercise the mixed tensor.
        lengths=[2+(i%2) for i in range(len(P))]
        expanded,O,_=cones.expand(rr,[i for i in range(len(rr)) if mm>>i&1],P,lengths)
        if len(expanded)<=15:
            h,u,w,_=tables(np.array(expanded,dtype=np.int64));literal=defect(h,w,O)
            need(literal==evaluate(E,lengths),'independent literal mixed control')
        else:literal=None
        out.append({'code':r['code'],'selected':r['selected'],'birth':not r['eligible_twins'],'stable_blocks':blocks,'reduced_rows':rr,'reduced_selected':mm,'clock_ports':P,'joint_stable':ok,'stirling_terms':[[*a,c] for a,c in sorted(K.items())],'exponential_terms':[[*a,c] for a,c in sorted(E.items())],'literal_lengths':lengths,'literal_defect':literal})
        if idx%20==0:print(json.dumps({'at':idx,'total':len(bank),'seconds':round(time.time()-start,2),'failures':len(fail)}),flush=True)
    prior=json.loads(gzip.decompress((HERE.parent/'five_stress_connection_atlas_20261112/descendants_nine.json.gz').read_bytes()))
    newkeys={key(r['rows'],r['selected']) for r in bank if r['eligible_twins']}
    oldkeys={key(r['rows'],r['selected']) for r in prior['records']}
    need(newkeys==oldkeys,'complete descendant graph-key equality')
    result={'scope':'Full all-length mixed-defect tensors for every strong marked order-nine unit seed with a stable unselected port. Stable twin-block compression is justified by identity on infinite product tails; no joint-stability or confluence hypothesis used.','states':len(bank),'joint_failures':fail,'clock_rank_joint_histogram':[[a,b,c] for (a,b),c in sorted(hist.items())],'birth_rank_histogram':dict(Counter(len(r['clock_ports']) for r in out if r['birth'])),'complete_prior_descendant_key_equality':True,'records':out}
    save('joint_nine.json.gz',result)
    print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
