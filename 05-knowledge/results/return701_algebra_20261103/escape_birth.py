"""All one-native/one-birth programs from the supplied701 source."""
from collections import Counter
import sys
from common import HERE,ROOT,need,junction,save,source701
from return_native import paths
sys.path.insert(0,str(ROOT/'05-knowledge/results/defect_offset_dynamics_20261103'))
from ear_cuts import ear_values,tables

def main():
    U,r,b,names=source701();n=len(U);hits=[];ready=[];checks=0;count=0;hist=Counter();good=[]
    for number,P in enumerate(paths(U)):
        V=junction.native(U,P);h=junction.H(V)
        if h>701:continue
        vals=ear_values(V);count+=len(vals)
        for mask,hh in enumerate(vals):
            if hh!=701:continue
            W=junction.ear(V,mask);need(junction.H(W)==701,'independent literal701 birth count');checks+=1
            prof=junction.profile(W,r);strong=junction.strong(W)
            row={'path':P,'intermediate_H':h,'birth_mask':mask,'rows':W,'retained_root_profile':prof,'strong':strong}
            hits.append(row);hist[prof[3]-701]+=1
            if strong and prof[3]-701==1:ready.append(row)
            if strong:
                for rr in range(len(W)):
                    pp=prof if rr==r else junction.profile(W,rr)
                    if pp[3]-701==1:good.append({**row,'ready_root':rr,'ready_root_profile':pp})
        tables.cache_clear()
        if number%100==0:print('native path',number,'701 births',len(hits),'same-root ready',len(ready),'any-root ready',len(good),flush=True)
    save('escape_birth.json',{'status':'FINITE_EXACT_ALL_NATIVE_ONE_BIRTH','source_rows':U,'source_names':names,'root':r,'spare':b,
      'source_paths':701,'evaluated_birth_masks':count,'pruning':'H(native)>701 excluded by proved insertion monotonicity',
      'literal_counts_checked':checks,'retained_root_defect_histogram':dict(sorted(hist.items())),'hits':hits,'ready_returns':ready,'remarked_ready_returns':good})
    print('PASS all native/birth programs: targets',len(hits),'same-root ready',len(ready),'any ready',len(good),flush=True)
if __name__=='__main__':main()
