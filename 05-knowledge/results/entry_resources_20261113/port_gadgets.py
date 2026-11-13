"""One-vertex rules seeing an unselected stable port's actual adjacency."""
from common import *
from itertools import product
from collections import Counter
import time
from stable_census import certify

def extend(rr,O,p,rule):
    fR,fO,edge,sel=rule;n=len(rr);out=list(rr)+[0]
    for i in range(n):
        if i==p:beats=edge
        else:
            f=fO if O>>i&1 else fR;a=(rr[p]>>i)&1
            beats=[0,1,a,1-a][f]
        if beats:out[n]|=1<<i
        else:out[i]|=1<<n
    return out,O|(sel<<n)

def run():
    bank=load(HERE.parent/'endpoint_coverage_compiler_20261111/inflation_bank_7.json')['records']
    rules=list(product(range(4),range(4),range(2),range(2)));alive=rules[:];fails={};start=time.time()
    for r in bank:
        rr,O,p=r['rows'],r['selected'],r['inflated_vertex'];new=[]
        for rule in alive:
            aa,oo=extend(rr,O,p,rule);h,d,u=selected_state(aa,oo)
            if d==1:new.append(rule)
            else:fails[rule]={'rows':rr,'selected':O,'port':p,'after_defect':d}
        alive=new
    records=[]
    for rule in alive:
        stable_counts=Counter();strong_counts=Counter();examples=[]
        for r in bank:
            rr,O,p=r['rows'],r['selected'],r['inflated_vertex'];aa,oo=extend(rr,O,p,rule)
            st,cert,_,_=certify(aa,[oo]);ports=st[oo]
            stable_counts['old_root']+=int(p in ports);stable_counts['new_root']+=int(len(rr) in ports)
            if strong(rr):
                strong_counts['input']+=1;strong_counts['output']+=int(strong(aa))
                if len(examples)<3:examples.append({'before':selected_state(rr,O),'after':selected_state(aa,oo),'rows':rr,'selected':O,'port':p})
        records.append({'rule':rule,'output_stable_counts':dict(stable_counts),'strong_counts':dict(strong_counts),'examples':examples})
    out={'scope':'All64 pointwise adjacency rules (constant0,constant1,copy,reverse independently on O and remaining vertices), either arc to old port, either new selected flag; all634 inherited stable-port presentations through7. Stable outputs have finite all-length certificates; universality beyond this bank is not inferred.','inputs':len(bank),'survivors':records,'failures':[{'rule':rule,**f} for rule,f in fails.items()]}
    save('port_gadgets',out);print(records,'seconds',round(time.time()-start,1),flush=True)
if __name__=='__main__':run()
