"""Focused extraction of strong five-port birth candidates, not a census."""
from common import *
from birth_extensions import quadratic
from extension_heads import add
from joint_nine import kernel,expo,cones
import time

def run():
    start=time.time();found={}
    for idx,r in enumerate(read('birth_extensions_12.json.gz')['records'][:26]):
        if len(r['stable_ports'])<3:continue
        n=len(r['rows']);O=r['selected'];counts=quadratic(r['rows'],O)
        for marked in [0,1]:
            sm=O|(marked<<n)
            for cut in np.flatnonzero(counts[marked]==1):
                rr=add(r['rows'],int(cut))
                if not strong(rr):continue
                P,cert,_,_=certify(rr,[sm]);P=P[sm]
                if len(P)<5:continue
                _,_,C=quotient(rr,sm)
                if any(len(c)>1 and set(c)&set(P) for c in C):continue
                k=key(rr,sm)
                if k in found:continue
                found[k]={'rows':rr,'selected':sm,'stable_ports':P,'parent_id':r['birth_id'],'cut':int(cut),'new_mark':marked,'acyclic_ports':sorted(sum(rr[v]>>w&1 for w in P) for v in P)==list(range(len(P)))}
                save('rank_candidates.json.gz',{'scope':'Focused candidates from first26 retained order12 parents with at least3 stable ports. No minimality or completeness claim.','records':list(found.values())})
                print(json.dumps(found[k]),flush=True)
    print('candidates',len(found),'seconds',time.time()-start,flush=True)
if __name__=='__main__':run()
