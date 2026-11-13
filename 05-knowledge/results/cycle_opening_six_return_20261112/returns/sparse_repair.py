"""Independent module audit and radius-three parent repair after ready advance."""
from pathlib import Path
from itertools import combinations
from collections import defaultdict, Counter
import json, time
import constructor_language as c
HERE=Path(__file__).resolve().parent
P=(78,12,24,16,35,79,30)

def inflate(rows,v,k):
    types=list(range(len(rows)))+[v]*(k-1)
    out=[0]*len(types)
    for i,j in combinations(range(len(types)),2):
        win=True if types[i]==types[j] else bool(rows[types[i]]>>types[j]&1)
        a,b=(i,j) if win else (j,i)
        out[a]|=1<<b
    return tuple(out)

def module_coeff(rows,module,start=None):
    """Collapsed path runs, coefficient of j! S(k,j), independently enumerated."""
    ids=[i for i in range(len(rows)) if i!=module]; m=len(ids); full=(1<<m)-1
    # State (outside mask,last,number of visits to repeated module).
    bank=[defaultdict(int) for _ in range(full+1)]
    if start is None or start==module:bank[0][(m,1)]=1
    for i,v in enumerate(ids):
        if start is None or start==v:bank[1<<i][(i,0)]=1
    for mask in range(full+1):
        # Last outside vertex can append a fresh module run at this same mask.
        for (last,j),number in list(bank[mask].items()):
            if last<m and rows[ids[last]]>>module&1:
                bank[mask][(m,j+1)]+=number
        # Every outside append increases mask; no adjacent repeated-module runs.
        for (last,j),number in bank[mask].items():
            v=module if last==m else ids[last]
            for i,u in enumerate(ids):
                if not(mask>>i&1) and rows[v]>>u&1:
                    bank[mask|1<<i][(i,j)]+=number
    out=Counter()
    for (last,j),number in bank[full].items():
        if j:out[j]+=number
    return dict(sorted(out.items()))

def exp_coeff(run):
    from math import comb
    out=Counter()
    for j,cnt in run.items():
        for i in range(1,j+1):out[i]+=cnt*(-1)**(j-i)*comb(j,i)
    return {j:v for j,v in sorted(out.items()) if v}

def flip(rows,edges):
    out=list(rows)
    for i,j in edges:out[i]^=1<<j;out[j]^=1<<i
    return tuple(out)

def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument("--radius",type=int,default=3)
    args=parser.parse_args()
    audit=[]
    for port in (0,4):
        child=c.birth(P,port)
        hc=exp_coeff(module_coeff(child,1))
        dc=exp_coeff(module_coeff(P,1,start=port))
        c.need(hc=={1:7,2:45},'exact all-height child H')
        c.need(dc==({1:1,2:6} if port==0 else {1:1}),'exact all-height defect')
        audit.append(dict(port=port,H=hc,defect=dc,H_runs=module_coeff(child,1),defect_runs=module_coeff(P,1,start=port)))
    literal=[]
    for k in range(1,6):
        par=inflate(P,1,k)
        for port in (0,4):
            child=c.birth(par,port);v=c.state(child,len(child)-1)
            c.need(v[0]+v[2]==45*2**k+7 and v[1]==(6*2**k+1 if port==0 else 1),'literal inflated neutral switch')
            literal.append(dict(k=k,port=port,state=v))
    # Start with ready port4, add two successive sole-predecessor births.
    ready=c.birth(P,4); X=c.birth(ready,7)
    edges=list(combinations(range(len(X)),2)); hits=[];stats=Counter();started=time.time()
    for radius in range(args.radius+1):
        for chosen in combinations(edges,radius):
            alt=flip(X,chosen)
            h,s,F,birthH=c.probe.profiles(alt)
            stats[radius]+=1
            for port,H in enumerate(birthH):
                if H!=293:continue
                # This scope includes every potential new-root defect, not just 1/2.
                if s[port]>2:continue
                child=c.birth(alt,port)
                hc=exp_coeff(module_coeff(child,1));dc=exp_coeff(module_coeff(alt,1,start=port))
                record=dict(radius=radius,flips=chosen,port=port,defect=s[port],H_coeff=hc,defect_coeff=dc,all_height_H=hc=={1:23,2:135},rows=child,strong=c.strong(child))
                # Replacing a chosen entry inside the module does not commute with inflation.
                record['module_port_compatible']=port!=1
                hits.append(record)
        print('radius',radius,'parents',stats[radius],'hits',len(hits),'seconds',round(time.time()-started,2),flush=True)
    out=dict(status='FINITE-EXACT bounded sparse repair scope plus exact module-polynomial audit',max_radius=args.radius,audit=audit,literal=literal,scope='all parent arc subsets of radius <=requested maximum in X=P+z+x (9 vertices), every new sole predecessor, child total H=293, defect<=2; base module vertex1',stats=dict(stats),hits=hits)
    (HERE/('sparse_repair.json' if args.radius==3 else 'sparse_repair_radius'+str(args.radius)+'.json')).write_text(json.dumps(out,indent=2)+'\n')
    print('PASS independent all-height neutral family; repair hits',len(hits))
if __name__=='__main__':main()
