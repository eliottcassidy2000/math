"""All-defect sparse next-target census; lift the smallest positive defect stratum."""
from itertools import combinations
from collections import Counter
import json,time
import sparse_repair as s
from next_repair import Q

def main():
    X=s.c.birth(Q,9);edges=list(combinations(range(11),2));all_hits=[];stats=Counter();tic=time.time()
    for radius in range(3):
        for flips in combinations(edges,radius):
            parent=s.flip(X,flips);h,starts,F,counts=s.c.probe.profiles(parent);stats[radius]+=1
            for port,H in enumerate(counts):
                if H!=881:continue
                child=s.c.birth(parent,port)
                if not s.c.strong(child):continue
                all_hits.append(dict(radius=radius,flips=flips,port=port,defect=starts[port],rows=child))
        print('radius',radius,'parents',stats[radius],'strong target hits',len(all_hits),'seconds',round(time.time()-tic,2),flush=True)
    minimum=min((h['defect'] for h in all_hits if h['defect']>0),default=None)
    selected=[h for h in all_hits if h['defect'] and (h['defect']<=6 or h['defect']==minimum)]
    print('minimum positive defect',minimum,'selected hits',len(selected),'spectrum',dict(Counter(h['defect'] for h in all_hits)),flush=True)
    lifted=[];checks=0
    for hit in selected:
        child=tuple(hit['rows']);parent=tuple(row&2047 for row in child[:11]);port=hit['port']
        for module in range(11):
            if module==port:continue
            checks+=1
            hc=s.exp_coeff(s.module_coeff(child,module));dc=s.exp_coeff(s.module_coeff(parent,module,start=port))
            if hc=={1:71,2:405}:lifted.append(dict(base=hit,module=module,H=hc,defect=dc))
    out=dict(status='FINITE-EXACT sparse all-defect census and exact all-height module checks in minimum/small-positive defect strata',scope='Q standard advance, all parent flip subsets radius<=2, every sole-entry root port, strong child H881; lift every vertex except root/predecessor only for positive defects<=6 or the minimum positive defect',stats=dict(stats),positive_minimum=minimum,defect_spectrum=dict(Counter(h['defect'] for h in all_hits)),all_hits=all_hits,module_checks=checks,lifted=lifted)
    (s.HERE/'growing_successor_probe.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('PASS module checks',checks,'uniform lifts',len(lifted),'seconds',round(time.time()-tic,2),flush=True)
if __name__=='__main__':main()
