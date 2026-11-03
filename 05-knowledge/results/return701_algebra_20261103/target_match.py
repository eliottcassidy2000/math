"""Distance to the actual2-ready701 source after a supplied native complement."""
from itertools import permutations
import numpy as np
from common import HERE,need,save,source701,junction
from ready701 import inv
from return_native import paths

def code(T):
    return sum(((T[i]>>j)&1)<<k for k,(i,j) in enumerate((i,j) for i in range(len(T)) for j in range(i+1,len(T))))

def main():
    target=None
    for T in inv.two_start(8):
        h,ss,q,F=inv.rooted_data(T,0)
        if h+F-2==701:target,rr=inv.birth(T,0);break
        inv.tables.cache_clear()
    need(target is not None and rr==8,'actual nine-vertex2-ready701 target')
    U,r,b,names=source701();ps=list(paths(U));protected={r,b,names.index(5),names.index(4),names.index(7)}
    native=[junction.native(U,P) for P in ps];best=100;witnesses=[];tested=0
    for deleted in range(10):
        if deleted in protected:continue
        keep=tuple(i for i in range(10) if i!=deleted);slot=keep.index(r)
        actual=[junction.sub(V,keep) for V in native];codes=np.array([code(V) for V in actual],dtype=np.uint64)
        oldslots=tuple(i for i in range(9) if i!=slot);batch=[];maps=[]
        def flush():
            nonlocal best,witnesses,tested
            if not batch:return
            ds=np.bitwise_count(np.array(batch,dtype=np.uint64)[:,None]^codes[None,:]);minimum=int(ds.min());tested+=ds.size
            if minimum<best:best=minimum;witnesses=[];print('new minimum',best,'delete',names[deleted],flush=True)
            if minimum==best:
                for i,j in np.argwhere(ds==best):
                    if len(witnesses)>=80:break
                    mapping=maps[int(i)];W=tuple(sum(((target[a]>>bb)&1)<<mapping[bb] for bb in range(9)) for a in range(9))
                    # W above is in target row order; move rows to the mapped indices.
                    Z=[0]*9
                    for a,row in enumerate(W):Z[mapping[a]]=row
                    V=actual[int(j)];flips=[(a,bb) for a in range(9) for bb in range(a+1,9) if ((V[a]>>bb)^(Z[a]>>bb))&1]
                    witnesses.append({'native_path':ps[int(j)],'deleted_index':deleted,'deleted_name':names[deleted],
                      'remaining_names':[names[t] for t in keep],'intermediate_rows':V,'target_rows':Z,'root':slot,
                      'spare':keep.index(b),'target_to_remaining':mapping,'flips':flips,'target_profile':junction.profile(tuple(Z),slot)})
            batch.clear();maps.clear()
        for perm in permutations(oldslots):
            mapping=perm+(slot,);Z=[0]*9
            for a in range(9):Z[mapping[a]]=sum(((target[a]>>bb)&1)<<mapping[bb] for bb in range(9))
            batch.append(code(Z));maps.append(mapping)
            if len(batch)==512:flush()
        flush()
    save('target_match.json',{'status':'FINITE_EXACT_MINIMUM_IN_FIXED_TARGET_BANK','source_rows':U,'source_names':names,
       'target_rows':target,'target_root':rr,'protected_indices':sorted(protected),'tested_pairs':tested,
       'minimum_arc_discrepancy':best,'witnesses_first80':witnesses})
    print('minimum old arc repairs',best,'tested',tested,flush=True)
if __name__=='__main__':main()
