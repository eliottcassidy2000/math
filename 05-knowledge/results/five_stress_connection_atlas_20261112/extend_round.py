"""Exact compiled continuation of the round census; explicit int64 bounds.

For n<=18 every partial path count is <=n!, and every two-path cover sum
is <=(n+1)n!. Even twice this bound is <2**63. No floating arithmetic.
Controls compare the full H,U,W arrays with the independent Python routine.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
import json,sys,time
import numpy as np
from numba import njit
from verify_tournaments import representatives,matrix,aut_shifts,move,orbit,changes,tc,need
HERE=Path(__file__).resolve().parent

@njit
def tables(rows):
    n=len(rows);size=1<<n;full=size-1
    start=np.zeros((size,n),dtype=np.int64)
    end=np.zeros((size,n),dtype=np.int64)
    for mask in range(1,size):
        for v in range(n):
            if not (mask>>v)&1:continue
            rest=mask^(1<<v)
            if rest==0:start[mask,v]=1;end[mask,v]=1
            else:
                for u in range(n):
                    if (rows[v]>>u)&1:start[mask,v]+=start[rest,u]
                    if (rows[u]>>v)&1:end[mask,v]+=end[rest,u]
    q=np.zeros((n,n),dtype=np.int64)
    for mask in range(1,full):
        rest=full^mask
        for i in range(n):
            a=end[mask,i]
            if a==0:continue
            for j in range(n):
                if start[rest,j]:q[i,j]+=a*start[rest,j]
    W=q+q.T;U=np.copy(start[full]);H=np.sum(end[full])
    for i in range(n):
        for j in range(n):U[j]+=q[i,j]
    weights=np.zeros(size,dtype=np.int64)
    masks=[]
    for mask in range(1,size):
        v=0
        while not (mask>>v)&1:v+=1
        rest=mask^(1<<v);value=weights[rest]
        for u in range(v+1,n):
            if (rest>>u)&1:value+=W[v,u]
        weights[mask]=value
        if value==H-1:masks.append(mask)
    return H,U,W,masks

def run(max_n=16):
    need(1<=max_n<=18,'proved fixed-width bound only through18')
    baseline=json.loads((HERE/'tournament_receipt.json').read_text())
    old={r['n']:r for r in baseline['levels']};levels=[];records=[];controls=0;t=time.time()
    for n in range(1,max_n+1):
        words=representatives(n);hist=Counter();parts=set();count=0;presentations=0
        for wi,w in enumerate(words):
            rows=matrix(w);H,U,W,masks=tables(np.array(rows,dtype=np.int64));H=int(H)
            if n<=7 or wi in (0,len(words)//2,len(words)-1):
                h,u,ww,*_=tc.path_tables(rows)
                need(H==h and U.tolist()==u and W.tolist()==ww,'compiled versus Python full table')
                controls+=1
            hist[H]+=1;auto=aut_shifts(w);seen=set();presentations+=len(masks)
            for mask in masks:
                cm=min(move(mask,k,n) for k in auto)
                if cm in seen:continue
                seen.add(cm);count+=1
                value=sum(int(U[v]) for v in range(n) if cm>>v&1);odd=value;j=0
                while odd%2==0:odd//=2;j+=1
                parts.add(odd)
                if n>max(old):
                    oo=orbit(w)
                    records.append({'n':n,'word':''.join(map(str,w)),'rows':rows,'H':H,'selected':cm,'U':value,'odd_part':odd,'v2_U':j,'parabola_core':any(changes(v)<=3 for v in oo),'cubic_core':any(changes(v)<=4 for v in oo),'aut_order':len(auto)})
        level={'n':n,'round_isomorphism_types':len(words),'unit_marked_isomorphism_types':count,'unit_mark_presentations_on_unmarked_reps':presentations,'H_histogram':dict(sorted(hist.items())),'odd_seed_parts':sorted(parts)}
        if n in old:
            for key in level:
                expected=old[n][key]
                if key=='H_histogram':expected={int(k):v for k,v in expected.items()}
                need(level[key]==expected,'complete previous level disagreement '+str((n,key)))
        levels.append(level)
        print(json.dumps({'n':n,'classes':len(words),'seeds':count,'elapsed':round(time.time()-t,2)}),flush=True)
    out={'scope':f'All round isomorphism classes and every nonempty selected set through order{max_n}; actual stabilizer quotient. Complete agreement with original Python receipt through order14.','arithmetic_bound':'n<=18; all path and cover accumulators bounded by2*(n+1)*n!<2**63; no float arithmetic','full_array_controls':controls,'levels':levels,'additional_seed_records':records}
    tag='.optimized' if not __debug__ else ''
    (HERE/('round_extension'+tag+'.json')).write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 16)
