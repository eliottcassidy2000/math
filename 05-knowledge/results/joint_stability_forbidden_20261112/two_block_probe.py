"""Exhaustive two-block obstruction probe above ordinary order-nine quotients.

Every strong order-eleven marked tournament with two different unselected
twin blocks can be obtained by doubling two vertices of an order-nine
tournament (other block members remain in the quotient). If both blocks are
stable, each singly doubled source is unit. This necessary L-shaped filter
cuts the search before full axis certification. No assumption that the
twice-contracted corner is unit is made.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import gzip,csv,json,time,sys
import numpy as np
from census_nine import tables,strong,need,defect,axis_values,save,ib
from joint_nine import cones,quotient
from tensor import kernel,expo
HERE=Path(__file__).resolve().parent

def decode(code,n=9):
    rows=[0]*n;b=n*(n-1)//2-1
    for i in range(n):
        for j in range(i+1,n):
            if code>>b&1:rows[i]|=1<<j
            else:rows[j]|=1<<i
            b-=1
    return rows

def run():
    source=HERE.parent/'blue_quotient_20261020b/n9_classes.tsv.gz';start=time.time();stats=Counter();records=[];failures=[];false_corners=[]
    with gzip.open(source,'rt') as f:
        for ix,item in enumerate(csv.DictReader(f,delimiter='\t'),1):
            rows=decode(int(item['code']))
            if not strong(rows):continue
            candidates={};cache={}
            for p in range(9):
                labels,rr=ib.inflate(rows,p,2);h,u,w,masks=tables(np.array(rr,dtype=np.int64));stats['single_clone_tables']+=1
                cache[p]=(labels,rr,int(h),w)
                for mask in masks:
                    if mask>>8:continue # both new clone labels are unselected
                    O=sum(1<<labels[i] for i in range(8) if mask>>i&1)
                    candidates.setdefault(O,[]).append(p)
            for O,ps in candidates.items():
                _,_,components=quotient(rows,O)
                component_of={v:i for i,c in enumerate(components) for v in c}
                for p,q in combinations(ps,2):
                    if O>>p&1 or O>>q&1:continue
                    # Exclude two representatives in the SAME maximal block.
                    if component_of[p]==component_of[q]:continue
                    stats['L_shape_pairs']+=1
                    selected=[i for i in range(9) if O>>i&1]
                    rr,oo,labels=cones.expand(rows,selected,[p,q],[2,2]);h,u,w,_=tables(np.array(rr,dtype=np.int64));stats['double_clone_tables']+=1
                    if defect(h,w,oo)!=1:continue
                    stats['unit_double_corners']+=1
                    passed=True;cert=[]
                    # Certify each direction at the singly doubled source;
                    # its positive-length family includes the double corner.
                    for fixed,free in [(p,q),(q,p)]:
                        labs,base,bh,bw=cache[fixed];v=labs.index(free);bm=sum(1<<i for i,a in enumerate(labs) if O>>a&1)
                        R=1+min(base[v].bit_count(),len(base)-1-base[v].bit_count())
                        vals=axis_values(base,bm,v,R);stats['axis_tables']+=R
                        cert.append({'fixed_doubled':fixed,'free':free,'R':R,'values':vals})
                        if any(d!=1 for d in vals):passed=False;break
                    if not passed:
                        false_corners.append({'quotient_code':int(item['code']),'rows':rows,'selected':O,'ports':[p,q],'axis_certificates':cert})
                        continue
                    stats['certified_stable_block_pairs']+=1
                    K=kernel(rows,O,[p,q]);joint=K=={(1,1):1}
                    entry={'quotient_code':int(item['code']),'rows':rows,'selected':O,'ports':[p,q],'axis_certificates':cert,'joint_stable':joint,'stirling_terms':[[*a,c] for a,c in sorted(K.items())],'double_corner_H':int(h),'double_corner_U':sum(int(u[i]) for i in oo)}
                    records.append(entry)
                    if not joint:
                        failures.append(entry);print('CONTRACTIBLE_OBSTRUCTION '+json.dumps(entry),flush=True)
            if ix%5000==0:print(json.dumps({'quotients':ix,'seconds':round(time.time()-start,1),'stats':dict(stats),'failures':len(failures)}),flush=True)
    out={'scope':'Every ordinary strong order-nine quotient, all selected sets, two different unselected blocks doubled. Equivalent complete search for stable port pairs in different contractible blocks at order eleven; triples and larger simultaneous stability are separate. Presentations are not quotiented by automorphisms.','ordinary_quotients':ix,'statistics':dict(stats),'false_unit_corner_screens':false_corners,'failures':failures,'records':records}
    save('two_block_probe.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k not in ['records','failures']}),flush=True)
if __name__=='__main__':run()
