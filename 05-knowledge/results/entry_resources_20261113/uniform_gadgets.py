"""Core-changing selected-cut gadgets, modulo isomorphism of the new core."""
from common import *
from itertools import permutations,product
from collections import Counter
import time

def canonical(rr):
    n=len(rr);best=None
    for p in permutations(range(n)):
        code=sum(((rr[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(( (i,j) for i in range(n) for j in range(i+1,n))))
        best=code if best is None else min(best,code)
    return best

def cores(q):
    seen=set()
    for code in range(1<<(q*(q-1)//2)):
        rr=labelled(q,code);c=canonical(rr)
        if c not in seen:seen.add(c);yield rr

def attach(rr,O,gg,types,s):
    n=len(rr);out=list(rr)+[r<<n for r in gg]
    for k,t in enumerate(types):
        v=n+k
        for i in range(n):
            if t>>(1 if O>>i&1 else 0)&1:out[v]|=1<<i
            else:out[i]|=1<<v
    return out,O|(s<<n)

def run(maxq=4):
    start=time.time();alive=[];universe=Counter()
    for q in range(1,maxq+1):
        for gg in cores(q):
            for types in product(range(4),repeat=q):
                for s in range(1<<q):alive.append((tuple(gg),types,s));universe[q]+=1
    fails=Counter();witnesses={};tested=0
    for n in range(2,6):
        for code in range(1<<(n*(n-1)//2)):
            rr=labelled(n,code);h,u,w,_=tables(np.array(rr,dtype=np.int64))
            for O in range(1,(1<<n)-1):
                if defect(h,w,[i for i in range(n) if O>>i&1])!=1:continue
                tested+=1;next_alive=[]
                for gg,ts,s in alive:
                    aa,oo=attach(rr,O,gg,ts,s);state=selected_state(aa,oo)
                    if state[1]==1:next_alive.append((gg,ts,s))
                    else:
                        fails[len(gg),n]+=1
                        # Every rejection receives its actual shortest-order source
                        # and after defect; groups share witness source data.
                        key=(tuple(rr),O)
                        witnesses.setdefault(key,[]).append([list(gg),list(ts),s,state[1]])
                alive=next_alive
        print('source order',n,'survivors',Counter(len(g) for g,t,s in alive),'seconds',round(time.time()-start,1),flush=True)
    records=[]
    for gg,ts,s in alive:
        # Padding means every new vertex is unselected and is either before
        # all old vertices or after all of them, consistently with new arcs.
        padding=s==0 and all(t in [0,3] for t in ts) and all(not(ts[i]==0 and ts[j]==3) or not(gg[i]>>j&1) for i in range(len(gg)) for j in range(len(gg)))
        transitive=sorted(r.bit_count() for r in gg)==list(range(len(gg)))
        records.append({'core':gg,'types':ts,'selected_new':s,'ordinary_transitive_padding':padding and transitive})
    out={'scope':f'Every new-core isomorphism class of orders1..{maxq}, every selected/complement uniform orientation for its vertices, every new selected subset. All nonempty proper selected unit sets in old labelled tournaments through5. Rejection certificates exhaustive for this finite gadget universe; surviving nonpadding patterns would still need a proof.','universe_by_order':dict(universe),'unit_source_presentations':tested,'survivors':records,'rejections_by_core_and_source_order':[[*k,v] for k,v in sorted(fails.items())],'rejection_groups':[{'rows':rr,'selected':O,'rejected':v} for (rr,O),v in witnesses.items()]}
    save('uniform_gadgets',out)
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 4)
