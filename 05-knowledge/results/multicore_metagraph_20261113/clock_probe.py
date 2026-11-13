"""Test every non-port vertex of the bounded ready-target parent bank."""
from pathlib import Path
import sys,json,gzip
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'berggren_target_return_20261111'/'reset'))
import funnel_six_state as f
sys.path.insert(0,str(HERE.parent/'endpoint_coverage_compiler_20261111'))
import inflation_bank as ib
from rooted_growth import need,save,endpoints

def birth(T,r):
    n=len(T)
    return tuple(row|((1<<n) if i==r else 0) for i,row in enumerate(T))+(((1<<n)-1)^(1<<r),)

def run():
    data=json.loads(gzip.decompress((HERE/'next_target.json.gz').read_bytes()))
    results=[]
    for rec in data['hits']:
        if rec['H']!=881:continue
        parent=f.a.unrank(rec['order'],rec['rank']);T=birth(parent,0);n=len(T)
        for v in range(n):
            labels,U=ib.inflate(T,v,2);r=labels.index(n-1)
            end=endpoints(U);h=sum(end)
            reverse=tuple(((1<<len(U))-1)^(1<<i)^row for i,row in enumerate(U))
            start=endpoints(reverse)[r]
            delete=f.a.induced(tuple(U),[i for i in range(len(U)) if i!=r])
            delta=sum(endpoints(delete))-start
            row=dict(parent_order=rec['order'],parent_rank=rec['rank'],vertex=v,rows=T,H2=h,defect2=delta)
            if h==1691 and delta==1:
                order,tab=ib.visit_ends(T,v)
                coeff=[0]+[sum(tab[-1][j]) for j in range(1,n+1)]
                row['full_run_coefficients']=coeff
                row['uniform_target']=coeff[1:3]==[881,405] and not any(coeff[3:])
            results.append(row)
    eligible=[row for row in results if row['H2']<=1691]
    need(len(eligible)==1 and eligible[0]['H2']==1559,'unique non-overshooting clone')
    rec=eligible[0];T=tuple(rec['rows']);v=rec['vertex']
    labels,U=ib.inflate(T,v,3);h3=sum(endpoints(U));need(h3>1691,'third clone overshoots')
    need(max(x['order'] for x in data['levels'])>=12,'Busch bound for all ready881 parents covered')
    save('clock_probe',dict(scope='Every vertex of every strong unit-ready881 class; completeness uses the published strong-order bound and READY classification; weak trailing transitive padding is handled in the proof',records=results,only_nonovershooting_clone=rec,third_clone_H=h3,lower_bound_H3=881+3*(1559-881)))
    print([(r['parent_order'],r['parent_rank'],r['vertex'],r['H2'],r['defect2']) for r in results]);print('hits',sum(r.get('uniform_target',False) for r in results))
if __name__=='__main__':run()
