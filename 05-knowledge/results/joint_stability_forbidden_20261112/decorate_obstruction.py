"""Complete one-vertex uniform-block extensions of the nonconfluent core."""
from pathlib import Path
from collections import Counter,defaultdict
from itertools import combinations
import json
import numpy as np
from census_nine import tables,defect,need,save,strong
from joint_nine import cones,kernel,expo,key
from contractible_obstruction import describe

def axis(E,index,frozen=2):
    c=Counter()
    for bases,value in E.items():c[bases[index]]+=value*bases[1-index]**frozen
    return {b:v for b,v in c.items() if v}=={1:1}

def run():
    Q=[2,100,49,7,75,409,173,31,223];P=[0,8];stats=Counter();records=[];memo={};details={}
    def terminal(rows,mask):
        k=key(rows,mask)
        if k in memo:return memo[k]
        d=describe(rows,mask);need(d['defect']==1,'recursive unit contraction');details[k]=d;roots=set()
        for block in d['eligible_blocks']:
            keep=[i for i in range(len(rows)) if i!=block[0]]
            rr=[sum(1<<j for j,v in enumerate(keep) if rows[u]>>v&1) for u in keep];mm=sum(1<<j for j,v in enumerate(keep) if mask>>v&1)
            roots.update(terminal(rr,mm))
        if not roots:roots={k}
        memo[k]=roots;return roots
    for sig in range(512):
        q=[r|((1<<9) if not sig>>i&1 else 0) for i,r in enumerate(Q)]+[sig]
        if not strong(q):continue
        for mark_new in [False,True]:
            O=50|((1<<9) if mark_new else 0);sel=[i for i in range(10) if O>>i&1];stats['strong_presentations']+=1
            corners=[];okay=True
            for sizes in [[1,2],[2,1],[2,2]]:
                rr,oo,_=cones.expand(q,sel,P,sizes);h,u,w,_=tables(np.array(rr,dtype=np.int64));d=defect(h,w,oo);corners.append(d)
                if d!=1:okay=False;break
            if not okay:continue
            stats['three_unit_corners']+=1;K=kernel(q,O,P);E=expo(K)
            if not axis(E,0) or not axis(E,1):continue
            stats['two_stable_blocks']+=1
            if K=={(1,1):1}:continue
            stats['mixed_failures']+=1
            target,oo,_=cones.expand(q,sel,P,[2,2]);m=sum(1<<i for i in oo);roots=terminal(target,m);merged=set()
            for k in roots:
                r=details[k];n=len(r['rows']);rev=[((1<<n)-1)^r['rows'][i]^(1<<i) for i in range(n)]
                merged.add(min(k,key(rev,r['selected'])))
            entry={'signature':sig,'mark_new':mark_new,'quotient_rows':q,'quotient_selected':O,'ports':P,'Stirling_terms':[[*a,c] for a,c in sorted(K.items())],'target':details[key(target,m)],'ordinary_terminal_births':len(roots),'converse_terminal_births':len(merged),'terminal_births':[details[k] for k in sorted(roots)]}
            records.append(entry);print(json.dumps({'signature':sig,'mark_new':mark_new,'ordinary_roots':len(roots),'converse_roots':len(merged)}),flush=True)
    out={'scope':'All 512 adjacency signatures and both marks for one new vertex on the order-nine nonunit quotient; then double its two named ports. Retain strong presentations. Full mixed kernels and complete terminal contraction sets. Not a complete order-twelve universe.','statistics':dict(stats),'records':records}
    save('decorated_obstructions.json.gz',out);print(json.dumps({'statistics':dict(stats),'counterexamples':len(records)}),flush=True)
if __name__=='__main__':run()
