"""Complete one-vertex extension step on a declared stable-birth layer.

The new vertex occurs once in each Hamilton word/two-path cover, hence its
orientation response is quadratic. Evaluate0,singletons,pairs, then all cuts.
Selections retain the parent marking, with either status for the new vertex.
"""
from common import *
from extension_heads import add
from itertools import combinations
from collections import Counter
from numba import njit
import time

@njit(cache=True)
def eval_cube(c,linear,pairs):
    n=len(linear);out=np.empty(1<<n,dtype=np.int64);out[0]=c
    for mask in range(1,1<<n):
        bit=mask&-mask;v=0
        while 1<<v!=bit:v+=1
        rest=mask^bit;d=out[rest]+linear[v]
        for w in range(v+1,n):
            if rest>>w&1:d+=pairs[v,w]
        out[mask]=d
    return out

def quadratic(rows,O):
    n=len(rows);selections=np.array([O,O|(1<<n)],dtype=np.int64)
    c=selected_defects(np.array(add(rows,0),dtype=np.int64),selections)
    linear=np.zeros((2,n),dtype=np.int64);pairs=np.zeros((2,n,n),dtype=np.int64)
    for v in range(n):linear[:,v]=selected_defects(np.array(add(rows,1<<v),dtype=np.int64),selections)-c
    for v,w in combinations(range(n),2):pairs[:,v,w]=selected_defects(np.array(add(rows,(1<<v)|(1<<w)),dtype=np.int64),selections)-c-linear[:,v]-linear[:,w]
    out=[eval_cube(c[q],linear[q],pairs[q]) for q in range(2)]
    return out

def run(n=10,min_rank=0,variant=''):
    if variant:
        bank=read(f'birth_sources_{variant}_{n}.json.gz')['records'];parent_complete=False
    elif n<=10:
        bank=[r for r in read('birth_union_ten.json.gz',True)['records'] if len(r['rows'])==n];parent_complete=True
    else:
        bank=read(f'birth_extensions_{n}.json.gz')['records'];parent_complete=False
    if min_rank:
        bank=[r for r in bank if len(r.get('ports',r.get('stable_ports',[])))>=min_rank];parent_complete=False
    found={};unit=set();stats=Counter();start=time.time();rank=Counter();parents=[];weak={}
    for idx,r in enumerate(bank):
        rows=r['rows'];O=r['selected'];counts=quadratic(rows,O);stats['orientation_mark_presentations']+=2**(n+1);stats['literal_interpolation_calls']+=1+n+n*(n-1)//2
        for cut in [0,(1<<n)-1,((idx+17)*73)%2**n,((idx+31)*277)%2**n]:
            got=selected_defects(np.array(add(rows,cut),dtype=np.int64),np.array([O,O|(1<<n)],dtype=np.int64));need([int(x) for x in got]==[int(c[cut]) for c in counts],'literal quadratic cube control');stats['literal_cube_controls']+=2
        parent_new=0
        for marked in [0,1]:
            sm=O|(marked<<n)
            for cut in np.flatnonzero(counts[marked]==1):
                cut=int(cut);rr=add(rows,cut);stats['unit_presentations']+=1
                need(int(selected_defects(np.array(rr,dtype=np.int64),np.array([sm],dtype=np.int64))[0])==1,'every selected unit consequence literal')
                # Exact row/mark duplicate screen first, canonical quotient later.
                raw=(tuple(rr),sm)
                if raw in unit:continue
                unit.add(raw)
                P,cert,_,_=certify(rr,[sm]);P=P[sm]
                if not P:continue
                _,_,comps=quotient(rr,sm);eligible=[c for c in comps if len(c)>1 and set(c)&set(P)]
                stats['stable_presentations']+=1
                if eligible:continue
                k=key(rr,sm)
                if not strong(rr):
                    need(cut in [0,(1<<n)-1],'strong parent has only source/sink weak extensions')
                    weak[k]=len(P);continue
                parent_new+=k not in found
                if k not in found:
                    h,u,w,_=tables(np.array(rr,dtype=np.int64));need(defect(h,w,[v for v in range(n+1) if sm>>v&1])==1,'independent full W control')
                    found[k]={'rows':rr,'selected':sm,'H':int(h),'U':sum(int(u[v]) for v in range(n+1) if sm>>v&1),'stable_ports':P,'axis_certificates':cert[sm],'canonical_key':k,'parent_birth_id':r['birth_id'],'new_vertex_out_mask':cut,'new_vertex_selected':marked,'presentations':[]}
                    rank[len(P)]+=1
                found[k]['presentations'].append([r['birth_id'],cut,marked])
        parents.append({'birth_id':r['birth_id'],'new_classes_first_seen':parent_new})
        if idx%25==0:print(json.dumps({'parents':idx,'total':len(bank),'new_births':len(found),'ranks':dict(rank),'stats':dict(stats),'seconds':round(time.time()-start,1)}),flush=True)
    for i,r in enumerate(found.values()):r['birth_id']=f'{n+1}:{i}';r['ports']=r['stable_ports']
    need(len(weak)==2*len(bank),'complete weak source/sink padding controls')
    out={'scope':f'Every one-vertex orientation extending the{len(bank)} retained strong order{n} stable-clone births, inherited old selection with both markings of the new vertex. Exact quadratic unit response, full all-length axis certificates, strong child, no eligible unselected contraction, ordinary marked canonical keys. Weak source/sink paddings counted separately. Complete within this extension grammar, not a complete order{n+1} birth census.','parent_order':n,'child_order':n+1,'parent_births':len(bank),'parent_layer_complete':parent_complete,'weak_padding_births':len(weak),'weak_padding_port_histogram':dict(Counter(weak.values())),'statistics':dict(stats),'distinct_births':len(found),'stable_port_histogram':dict(rank),'parents':parents,'records':list(found.values())}
    if min_rank:out['scope']+=f' Parent restriction: at least{min_rank} stable ports in the retained layer.';out['parent_minimum_rank']=min_rank
    suffix=f'_rank{min_rank}' if min_rank else ''
    if variant:
        suffix+='_'+variant;out['scope']+=f' Parent source variant: {variant}.'
    save(f'birth_extensions_{n+1}{suffix}.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k not in ['parents','records']}),flush=True)
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 10,int(sys.argv[2]) if len(sys.argv)>2 else 0,sys.argv[3] if len(sys.argv)>3 else '')
