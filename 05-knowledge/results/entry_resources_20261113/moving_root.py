"""Attempt renewable one-exit core entry and audit every re-marking."""
from common import *
from numba import njit
from math import factorial

@njit(cache=True)
def unit_masks(h,w,n):
    W=np.zeros(1<<n,dtype=np.int64);out=[]
    for mask in range(1,1<<n):
        bit=mask&-mask;i=0;t=bit
        while t>1:i+=1;t>>=1
        smaller=mask^bit;value=W[smaller]
        for j in range(n):
            if smaller>>j&1:value+=w[i,j]
        W[mask]=value
        if h-value==1:out.append(mask)
    return out

def run():
    rr=[32,5,33,7,15,986,415,287,543,3295,2559,511];O=10;root=0;records=[]
    for k in range(8):
        need(2*(len(rr)+1)*factorial(len(rr))<2**63,'all path/cover accumulators safely fit int64')
        h,u,w,_=tables(np.array(rr,dtype=np.int64));n=len(rr)
        masks=list(map(int,unit_masks(h,w,n)))
        need((O in masks)==(defect(h,w,[1,3])==1),'inherited matrix is already symmetrized')
        # Exact nonnegative pair weights: adding markings cannot repair delta<1.
        records.append({'step':k,'rows':rr,'fixed_selected':O,'fixed_state':selected_state(rr,O),'all_unit_selected_masks':masks,'unit_markings_retaining_vertex1':[s for s in masks if s&2],'unit_markings_retaining_vertex3':[s for s in masks if s&8]})
        print(k,n,'fixed',records[-1]['fixed_state'],'all markings',masks[:30],'total',len(masks),flush=True)
        rr=[r|((1<<n) if v!=root else 0) for v,r in enumerate(rr)]+[1<<root];root=n
    need(records[0]['fixed_state'][1]==records[1]['fixed_state'][1]==1 and records[2]['fixed_state'][1]==-31,'first core move succeeds; second loses selected gate')
    save('moving_root',{'scope':'Eight actual one-exit root-changing extensions from the inherited12-vertex parent. All selected subsets tested at each order12..19; exact finite experiment, not an all-order obstruction.','records':records})
if __name__=='__main__':run()
