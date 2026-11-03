"""Actual native-complement endpoints and deletion returns at source701."""
from itertools import combinations,product
from collections import Counter
from common import HERE,old,junction,need,save,source701
from return_native import paths

def main():
    U,r,b,names=source701();n=len(U);ps=list(paths(U));need(len(ps)==701,'701 literal source paths')
    ready=[];hist=Counter();same=[];deleted=[];root_ready=[];seen={}
    protected={r,b,names.index(5),names.index(4),names.index(7)}
    for number,P in enumerate(ps):
        V=junction.native(U,P);h=junction.H(V);hist[h]+=1
        if h==701:same.append({'path':P,'rows':V})
        if h<=701:
            prof=junction.profile(V,r)
            if prof[3]-h==1:root_ready.append({'path':P,'rows':V,'H':h,'profile':prof})
        # A return by deletion may have H(V)>701; no H filter is used here.
        for k in range(n):
            if k in protected:continue
            keep=tuple(i for i in range(n) if i!=k);W=junction.sub(V,keep)
            if junction.H(W)!=701:continue
            rr=keep.index(r);prof=junction.profile(W,rr)
            row={'path':P,'deleted_vertex':k,'deleted_name':names[k],'rows':W,'profile':prof,'root':rr}
            deleted.append(row)
            if prof[3]-701==1:ready.append(row)
    save('escape_probe.json',{'status':'FINITE_EXACT_NATIVE_AND_ONE_DELETION_BANK','source_rows':U,'source_names':names,
      'root':r,'spare':b,'protected_vertices':sorted(protected),'paths':len(ps),'count_histogram':dict(sorted(hist.items())),
      'same_count_complements':same,'native_ready_roots_below701':root_ready,'same_count_after_deletion':deleted,'ready_returns':ready})
    print('native701 endpoints',len(same),'root-ready<=701',len(root_ready),'701 deletions',len(deleted),'ready',len(ready),flush=True)
if __name__=='__main__':main()
