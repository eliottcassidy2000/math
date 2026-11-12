"""Search one-vertex extensions of the actual even-mode source cores.

Every orientation is screened by the exact attachment quadratic form.
Positive candidates require equality of the entire Hamilton skeleton.
This is a restricted physical repair search, not a completeness claim.
"""
from common import *
from hamilton_clocks import skeleton,trim
from numba import njit
from collections import Counter
from math import gcd
import time

@njit(cache=True)
def end_counts(rows):
    n=len(rows);N=1<<n;d=np.zeros((N,n),dtype=np.int64)
    for v in range(n):d[1<<v,v]=1
    for m in range(1,N):
        for v in range(n):
            if not m>>v&1:continue
            rest=m^(1<<v)
            for w in range(n):
                if rows[w]>>v&1:d[m,v]+=d[rest,w]
    return d[N-1].copy()

@njit(cache=True)
def attachments(H,U,W,E):
    n=len(U);out=np.zeros(1<<n,dtype=np.int64);linear=np.zeros(1<<n,dtype=np.int64);tax=np.zeros(1<<n,dtype=np.int64);out[0]=H
    for m in range(1,1<<n):
        b=m&-m;v=0
        while 1<<v!=b:v+=1
        rest=m^b;linear[m]=linear[rest]+U[v]-E[v];tax[m]=tax[rest]
        for w in range(v+1,n):
            if rest>>w&1:tax[m]+=W[v,w]
        out[m]=H+linear[m]-tax[m]
    return out

def add(rows,O):
    n=len(rows)
    return [r|(0 if O>>v&1 else 1<<n) for v,r in enumerate(rows)]+[O]

def run():
    records=read('even_valuations.json.gz')['records'];out=[];stats=Counter();start=time.time();cache={}
    factor_graphs={1:[0]}
    from hamilton_clocks import graphs
    for rows in graphs(8):
        h=int(tables(np.array(rows,dtype=np.int64))[0])
        if h not in factor_graphs:factor_graphs[h]=rows
    for idx,r in enumerate(records):
        rows=r['rows'];p=r['port'];e=r['exponential_coefficients'];rowkey=tuple(rows)
        if rowkey not in cache:
            a=np.array(rows,dtype=np.int64);h,u,w,_=tables(a);E=end_counts(a);counts=attachments(h,u,w,E);cache[rowkey]=counts;stats['source_graphs']+=1;stats['orientation_presentations']+=len(counts)
        counts=cache[rowkey];matches=[]
        for d in range(9):
            if any(c%2**d for c in e):continue
            target=[c//2**d for c in e];target[0]+=1;content=abs(gcd(*target))
            for G in [g for g in factor_graphs if content%g==0]:
                wanted=trim([c//G for c in target]);h1=value(wanted,1)
                if h1<=0 or h1%2==0:continue
                for O in np.flatnonzero(counts==h1):
                    O=int(O);rr=add(rows,O);u=[int(x) for x in skeleton(np.array(rr,dtype=np.int64),p)];stats['whole_kernel_tests']+=1
                    if trim(expcoeff(u))==wanted:
                        matches.append({'division_exponent':d,'factor':G,'factor_rows':factor_graphs[G],'rows':rr,'port':p,'new_vertex_out_mask':O,'H_coefficients':list(wanted),'Stirling_coefficients':u})
        if matches:out.append({'birth_id':r['birth_id'],'port':p,'source_rows':rows,'source_selected':r['selected'],'type':r['type'],'repairs':matches})
        if idx%50==0:print(json.dumps({'at':idx,'axes_with_repairs':len(out),'stats':dict(stats),'seconds':round(time.time()-start,1)}),flush=True)
    result={'scope':'All478 even higher-mode source axes, every one-new-vertex orientation on each actual source graph, fixed inherited clock port, division exponents0..8 when the whole exponential coefficients are integral, and every realizable factor through8 dividing the profile content. Every positive has full exact H-skeleton identity. No claim that all lower heads are repaired.','statistics':dict(stats),'axes_with_repairs':len(out),'records':out}
    save('extension_heads.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
