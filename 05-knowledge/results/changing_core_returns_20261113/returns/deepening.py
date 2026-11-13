"""Reproduce the two deeper unit-ready hits, scoped miss, and defect-three repair."""
from pathlib import Path
import argparse,json,time
import inverse as inv
HERE=Path(__file__).resolve().parent

def ready_record(row,n,k):
    T,y,V,x,w=inv.physical(row,n)
    state=inv.f.advance(row['state'],w);N=405*(1<<k)+71
    inv.need(state[0]+state[4]-1==N,'deeper exact full count')
    inv.need(inv.f.a.strong(T),'deeper physical strongness')
    inv.need(sum(bool(T[v]>>y&1) for v in range(len(T)))==1,'physical indegree one')
    actual=inv.f.direct_state(inv.literal_core(row))
    inv.need(actual==row['state'],'selected core independent bins')
    inv.f.a.tables.cache_clear()
    if len(T)<=18:
        inv.need(inv.direct(T)==N and inv.direct(V,x)==1,'deeper literal full count/start')
    return dict(k=k,H=N,core=row,tail=n,word=''.join(map(str,w)),rows=T,root=y,predecessor=x,order=len(T),state=(N,state[0]-1,state[0],N+1),literal_checked=len(T)<=18)

def main(max_core):
    selected=[]
    for k in (8,21):
        for row in inv.cores(12):
            h,e,s,p,f,t=row['state'];G=s+p
            ns,_=inv.inverse(405*(1<<k)+71,G,h,e,G-e-f-1)
            if ns:
                selected.append(ready_record(row,min(ns),k));break
        else:raise RuntimeError('expected deeper hit missing')
    N=53084231;tested=1;miss=True;t0=time.time()
    ns,_=inv.inverse(N,3,1,1,0);inv.need(not ns,'singleton miss')
    for n in range(3,max_core+1):
        for rank in range(inv.f.a.class_counts(n)[n]):
            h,e,s,p,f,t=inv.f.compile_rank(n,rank);G=s+p
            ns,_=inv.inverse(N,G,h,e,G-e-f-1);tested+=1
            if ns:miss=False;raise RuntimeError(('unexpected hit',n,rank,ns))
        print('complete unit core order',n,'classes',tested,'seconds',round(time.time()-t0,2),flush=True)
    # A five-core with three start paths; an empty inner funnel makes the
    # six-core required by the successful outer inverse.
    base=(16,9,3,5,14);base_root=1
    p,e,s=inv.f.a.bin_values(base,base_root,3)
    base_state=(p[1],e[1],s[2],p[2],e[2],s[3]);b=s[1]
    core,croot=inv.f.a.funnel(base,(base_root,),())
    core_state=inv.f.advance(base_state,())
    p,e,s=inv.f.a.bin_values(core,croot,3)
    actual=(p[1],e[1],s[2],p[2],e[2],s[3])
    inv.need(core_state==actual==(31,12,34,328,89,183),'literal six-core full profile')
    inv.need(b==s[1]==inv.direct(core,croot)==3,'literal three-start propagation')
    n=145623;w=tuple((n>>i)&1 for i in range(n.bit_length()))
    parent,x=inv.f.a.funnel(core,(croot,),w);T,y=inv.f.a.birth(parent,x)
    state=inv.f.advance(core_state,w);h2,F=state[0],state[4]
    inv.need(h2+F-b==N,'repaired complete H')
    rooted=(N,h2-b,h2,N+b)
    inv.need(rooted==(53084231,5319648,5319651,53084234),'repaired marked state')
    inv.need(inv.f.a.strong(T),'actual repaired graph strong')
    inv.need([v for v in range(len(T)) if T[v]>>y&1]==[x],'unique actual predecessor')
    # A unit-start core with the same scalar response would need E2=91,
    # since the observer retains E2+start, not E2 by itself.
    shadows=[]
    for row in inv.cores(7):
        h,e,s,p,f,t=row['state']
        if (h,e,s+p,f+1)==(31,12,362,92):shadows.append(row)
    inv.need(not shadows,'no unit-start observer shadow through order7')
    repair=dict(k=17,H=N,base_rows=base,base_root=base_root,base_state=base_state,core_rows=core,core_root=croot,core_state=core_state,start_multiplicity=b,tail=n,word=''.join(map(str,w)),parent_state=state,rows=T,root=y,predecessor=x,order=len(T),state=rooted,readiness_defect=b,literal_full_graph_checked=False)
    data=dict(status='PASS',deeper_unit_ready=selected,unit_core_miss=dict(k=17,target=N,max_core_order=max_core,rooted_classes=tested,arbitrary_tail_length=True),defect_three_repair=repair,unit_observer_shadow=dict(max_order=7,observer=(31,12,362,92),hits=shadows))
    suffix='' if __debug__ else '.optimized'
    (HERE/('deepening_receipt'+suffix+'.json')).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS deeper witnesses, scoped unit miss, physical defect3 repair')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--max-core',type=int,default=17);main(p.parse_args().max_core)
