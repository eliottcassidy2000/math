"""Exhaustive two-vertex selected-cut gadget screen; not a universal proof."""
from common import *
from itertools import product
import time

def attach(rr,O,types,new_selected):
    n=len(rr);out=list(rr)+[0]*len(types)
    for k,t in enumerate(types):
        v=n+k
        for i in range(n):
            if t>>(1 if O>>i&1 else 0)&1:out[v]|=1<<i
            else:out[i]|=1<<v
        for j in range(k):out[n+j]|=1<<v
    return out,O|(new_selected<<n)

def run():
    candidates=[(ts,s) for size in [1,2] for ts in product(range(4),repeat=size) for s in range(1<<size)]
    alive=set(candidates);failures={};samples=[];tested=0;start=time.time()
    for n in range(2,6):
        for code in range(1<<(n*(n-1)//2)):
            rr=labelled(n,code);h,u,w,_=tables(np.array(rr,dtype=np.int64))
            for O in range(1,(1<<n)-1):
                picks=[i for i in range(n) if O>>i&1]
                if defect(h,w,picks)!=1:continue
                before=(int(h),1,sum(int(u[i]) for i in picks));tested+=1
                if len(samples)<30:samples.append((rr,O,before))
                for ts,s in sorted(alive):
                    aa,oo=attach(rr,O,ts,s);after=selected_state(aa,oo)
                    if after[1]!=1:
                        alive.remove((ts,s));failures[ts,s]={'rows':rr,'selected':O,'before':before,'after':after,'extension_rows':aa,'extension_selected':oo}
        print('order',n,'unit_sources',tested,'survivors',len(alive),'seconds',round(time.time()-start,1),flush=True)
    records=[]
    for ts,s in sorted(alive):
        examples=[]
        for rr,O,before in samples:
            aa,oo=attach(rr,O,ts,s);after=selected_state(aa,oo)
            examples.append({'before':before,'after':after,'strong':strong(aa)})
        records.append({'types':ts,'new_selected':s,'examples':examples})
    out={'scope':'All 72 uniform selected/complement-cut gadgets adding one or two vertices (new vertices in a transitive order), all nonempty proper selected unit sets in labelled tournaments of orders 2 through 5. Passing is finite evidence only.','unit_presentations':tested,'survivors':records,'first_failures':[{'types':ts,'new_selected':s,**f} for (ts,s),f in failures.items()]}
    save('cut_gadgets',out)
if __name__=='__main__':run()
