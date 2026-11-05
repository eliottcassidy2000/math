#!/usr/bin/env python3
import sys,json,math
from pathlib import Path
from functools import lru_cache
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from verify import *
from verify_unit_offset_rigidity import graph
from verify_three_module_clock import strong
SRC=Path(__file__).resolve().parents[2]/'session197-prime-atoms/positive-five-clock'

def run_poly(a,r):
    outside=[v for v in range(len(a)) if v!=r];n=len(outside);d={}
    for i,v in enumerate(outside):d[1<<i,i]=[1,int(bool(a[r]>>v&1))]
    for mask in range(1,1<<n):
        for i,u in enumerate(outside):
            p=d.get((mask,i))
            if p is None:continue
            for j,v in enumerate(outside):
                if mask>>j&1:continue
                edge=[int(bool(a[u]>>v&1)),int(bool((a[u]>>r&1) and (a[r]>>v&1)))]
                if any(edge):
                    key=mask|1<<j,j;d[key]=add(d.get(key,[0]),mul(p,edge))
    out=[0]
    for i,v in enumerate(outside):out=add(out,mul(d.get(((1<<n)-1,i),[0]),[1,int(bool(a[v]>>r&1))]))
    return trim(out)

@lru_cache(None)
def prime(q):
    return q>=2 and all(q%d for d in range(2,math.isqrt(q)+1))

def run():
    w=json.loads((SRC/'witness.json').read_text());U=(8,1,3,6);P=R=transitive(2);V=(0,)
    source=graph((0,),P,U,R,V);runs=run_poly(source,0)
    require(runs[1:]==[3461,3456,1152],'independent edge-polynomial run DP')
    require(all(c==0 for c in w['root_run_coefficients'][len(runs):]),'source higher runs')
    require(runs==w['root_run_coefficients'][:len(runs)],'complete source run vector')
    physical=[]
    for m in range(1,5):
        a=graph(transitive(m),P,U,R,V);H=value(packing_polynomial(a),2)
        require(list(a)==w['physical_controls'][m-1]['canonical_rows'],'independent block graph')
        require(H==1152*3**m+5 and strong(a),'independent odd-cycle count')
        physical.append((m,len(a),H))
    ports={8,9,10,11};birth=list(source);edge_checks=0
    for m in range(1,9):
        a=graph(transitive(m),P,U,R,V);order=list(range(12+m-2,11,-1))+[0]+list(range(1,12))
        require(len(order)==len(a),'birth permutation size')
        for i in range(len(a)):
            for j in range(len(a)):
                require((a[i]>>j&1)==(birth[order[i]]>>order[j]&1),'retained four-port all-edge mapping');edge_checks+=1
        n=len(birth);birth.append(sum(1<<v for v in range(n) if v not in ports))
        for v in ports:birth[v]|=1<<n
    prime_checks=[]
    for q in (3461,13331):
        ps=[p for p in range(2,math.isqrt(q)+1) if prime(p)]
        require(all(q%p for p in ps),'independent prime-only division')
        prime_checks.append((q,math.isqrt(q),len(ps)))
    profiles=set();h5={};graphs=0
    for n in range(1,7):
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask);hs=h_subsets(a);h=hs[-1];full=(1<<n)-1;b=sum(hs[s]*hs[full^s] for s in range(full+1))
            profiles.add((n,h,b));graphs+=1
            if h==5:
                h5[n]=h5.get(n,0)+1;require(b==32*2**(n-4),'independent H5 profile inventory')
    require(h5=={4:24,5:240,6:2160},'H5 labelled core/padding inventory')
    checks=pure=0
    for pp,uu,rr,vv in product(sorted(profiles),repeat=4):
        np,hp,bp=pp;nu,hu,bu=uu;nr,hr,br=rr;nv,hv,bv=vv
        A=2*bu*bv*(bp-hp)*(br-hr)
        B=bv*(br-hr)*(4*hp-bp)*hu+bu*(bp-hp)*(4*hr-br)*hv
        checks+=1
        if B:continue
        C=hu*hv*(bp-3*hp)*(br-3*hr)
        if C>=5 and prime(C):
            require(pp==(2,1,4) and rr==(2,1,4),'prime-offset leaf descent')
            require((hu==C and hv==1) or (hv==C and hu==1),'prime-offset seed descent');pure+=1
            if C==5:
                k=A//1152;require(A%1152==0 and k>0 and k&(k-1)==0,'exact +5 amplitude form')
    print(json.dumps({'status':'PASS','independent_run_polynomial':runs,'odd_cycle_physical_counts':physical,'retained_port_edge_checks':edge_checks,'prime_only_certificates':prime_checks,'labelled_graphs':graphs,'profile_types':len(profiles),'H5_inventory':h5,'coefficient_cases':checks,'pure_prime_offset_cases':pure},sort_keys=True,indent=2))
if __name__=='__main__':run()
