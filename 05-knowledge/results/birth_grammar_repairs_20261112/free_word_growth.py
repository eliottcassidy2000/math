"""Search insertion of a new independent stable port into a transitive port word."""
from common import *
from birth_extensions import quadratic
from extension_heads import add
import time
from collections import Counter

def signature(rows,P):
    fixed=[v for v in range(len(rows)) if v not in P]
    order=sorted(P,key=lambda v:-sum(rows[v]>>w&1 for w in P))
    need([sum(rows[v]>>w&1 for w in P) for v in order]==list(reversed(range(len(P)))),'acyclic word')
    return {'fixed_vertices':fixed,'port_order':order,'port_out_letters':[sum((rows[v]>>w&1)<<i for i,w in enumerate(fixed)) for v in order]}

def run():
    root=read('rank_candidates.json.gz')['records'][0];front=[root];records=[];start=time.time();stats=Counter();rejected=[]
    for level in range(3):
        nxt={}
        for r in front[:12]:
            rows=r['rows'];O=r['selected'];P=r['stable_ports'];n=len(rows);counts=quadratic(rows,O)[0]
            for cut in np.flatnonzero(counts==1):
                stats['unit_orientations']+=1
                cut=int(cut);rr=add(rows,cut);pp=P+[n]
                if sorted(sum(rr[v]>>w&1 for w in pp) for v in pp)!=list(range(len(pp))):continue
                if not strong(rr):continue
                stats['strong_acyclic_candidates']+=1
                stable,cert,_,_=certify(rr,[O]);stable=stable[O]
                if not set(pp)<=set(stable):continue
                stats['all_old_and_new_ports_stable']+=1
                _,_,C=quotient(rr,O)
                if any(len(c)>1 and set(c)&set(stable) for c in C):
                    rejected.append({'cut':cut,'eligible_blocks':[c for c in C if len(c)>1 and set(c)&set(stable)]});continue
                k=key(rr,O)
                if k in nxt:continue
                nxt[k]={'rows':rr,'selected':O,'stable_ports':stable,'parent_rows':rows,'parent_ports':P,'cut':cut,'signature':signature(rr,stable)}
            print(json.dumps({'order':n+1,'retained':len(nxt),'seconds':round(time.time()-start,2)}),flush=True)
        records+=list(nxt.values());front=list(nxt.values())
        save('free_word_growth.json.gz',{'scope':'Up to3 unselected one-vertex steps from the explicit five-port birth, retaining every old stable port and requiring an acyclic enlarged set. At most12 parents per next layer. Individual stability; joint stability requires its separate tensor certificate.','root_signature':signature(root['rows'],root['stable_ports']),'statistics':dict(stats),'rejected_contractions':rejected,'records':records,'final_layer':len(front)})
        if not front:break
if __name__=='__main__':run()
