"""Complete homogeneous-inflation test below3H when all pair-run costs are positive."""
from pathlib import Path
import sys,json,gzip
import numpy as np
from numba import njit
HERE=Path(__file__).resolve().parent
from clock_probe import f,birth
from rooted_growth import need,save,endpoints

@njit
def end_table(rows):
    n=len(rows);d=np.zeros((1<<n,n),dtype=np.int64)
    for i in range(n):d[1<<i,i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask>>v&1):continue
            rest=mask^(1<<v)
            if rest:
                for u in range(n):
                    if rows[u]>>v&1:d[mask,v]+=d[rest,u]
    return d

@njit
def clone_counts(rows):
    n=len(rows);full=(1<<n)-1
    rev=np.array([full^(1<<i)^rows[i] for i in range(n)],dtype=np.int64)
    end=end_table(rows);start=end_table(rev);q=np.zeros((n,n),dtype=np.int64)
    for mask in range(1,full):
        for i in range(n):
            if end[mask,i]:
                for j in range(n):q[i,j]+=end[mask,i]*start[full^mask,j]
    counts=[]
    for v in range(n):
        mask=rows[v];h=0
        for i in range(n):
            h+=start[full,i] if mask>>i&1 else end[full,i]
            if not(mask>>i&1):
                for j in range(n):
                    if mask>>j&1:h+=q[i,j]
        counts.append(h)
    return end[full].sum(),counts

def inflate(T,sizes):
    owners=[i for i,m in enumerate(sizes) for _ in range(m)];out=[0]*len(owners)
    for i,a in enumerate(owners):
        for j in range(i+1,len(owners)):
            b=owners[j];u,v=(i,j) if a==b or T[a]>>b&1 else (j,i);out[u]|=1<<v
    return tuple(out)

def test(T,target):
    h,counts=clone_counts(np.array(T,dtype=np.int64));h=int(h);counts=list(map(int,counts))
    need(h<target<3*h,'substitution target window')
    costs=[(v-h)//2 for v in counts];need(all(v>=0 for v in costs),'nonnegative pair-run cost')
    budget=target-h;options=[]
    for a in costs:
        if a==0:
            options.append([(1,0)]);continue # singleton SCC: transitive padding has no effect
        m=1;values=[]
        while a*((1<<m)-2)<=budget:
            values.append((m,a*((1<<m)-2)));m+=1
        options.append(values)
    candidates=0;hits=[];spectrum=[]
    def rec(i,left,sizes):
        nonlocal candidates
        if i==len(T):
            if all(s==1 for s in sizes):return
            U=inflate(T,sizes);need(len(U)<=20,'exact int64 audit order bound')
            value=int(end_table(np.array(U,dtype=np.int64))[-1].sum());candidates+=1
            spectrum.append(dict(sizes=sizes,H=value))
            if value==target:hits.append(dict(sizes=sizes,rows=U,H=value))
            return
        for m,cost in options[i]:
            if cost<=left:rec(i+1,left-cost,sizes+[m])
    rec(0,budget,[])
    return dict(rows=T,H=h,target=target,clone2_counts=counts,pair_run_costs=costs,candidates=candidates,spectrum=spectrum,hits=hits)

def run():
    old=json.loads(gzip.decompress((HERE/'clock_probe.json.gz').read_bytes()))
    d=json.loads(gzip.decompress((HERE/'next_target.json.gz').read_bytes()))
    records=[]
    from component_growth import components
    from rooted_growth import cc
    root_data=json.loads(gzip.decompress((HERE/'rooted_growth.json.gz').read_bytes()))
    small=sorted({(r['n'],r['ordinary_code']) for r in root_data['records'] if r['n']<=7});zero_controls=0
    for n,code in small:
        T=tuple(cc.decode(n,code));h,vals=clone_counts(np.array(T,dtype=np.int64));blocks=components(T)
        singleton={b[0] for b in blocks if len(b)==1}
        need({i for i,value in enumerate(vals) if value==h}==singleton,'zero clone cost iff singleton SCC')
        zero_controls+=n
    for r in d['hits']:
        if r['H'] not in (881,1691,3311,6551,13031,25991):continue
        need(next(b['complete'] for b in d['bounds'] if b['H']==r['H']),'all-order source head covered')
        T=birth(f.a.unrank(r['order'],r['rank']),0);value=test(T,2*r['H']-71)
        need(value['H']==r['H'],'compiler versus independent subset DP')
        if r['H']==881:
            prior=[s for s in old['records'] if s['parent_order']==r['order'] and s['parent_rank']==r['rank']]
            need(value['clone2_counts']==[x['H2'] for x in prior],'independent literal clone table')
        records.append(dict(parent_order=r['order'],parent_rank=r['rank'],**value))
        print(r['H'],r['order'],r['rank'],'eligible',sum(x<=value['target'] for x in value['clone2_counts']),'candidates',value['candidates'],'hits',len(value['hits']),flush=True)
    save('inflation_descent',dict(scope='Every strong unit-ready source at the first six405-ray values, complete by READY plus Busch; arbitrary homogeneous substitution membership for the next405-ray value; finite transitive size vectors exhaust the sub3H window',small_ordinary_classes=len(small),zero_cost_vertex_controls=zero_controls,records=records))
if __name__=='__main__':run()
