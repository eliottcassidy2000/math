"""Selected-defect engine and complete stable-port certification at n=9,10.

This engine sums selected boundary vectors directly, avoiding the full W
matrix. Its subset loops use only present neighbours. It is independently
checked against the full table engine and the complete order-nine census.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import sys,json,time,gzip,hashlib
import numpy as np
from numba import njit
from census_nine import read,save,need,tables,ib,defect
from joint_nine import quotient
from cyclic_port_probe import triangles
HERE=Path(__file__).resolve().parent

@njit(cache=True)
def selected_defects(rows,selections):
    n=len(rows);N=1<<n;full=N-1
    e=np.zeros((N,n),dtype=np.int64);s=np.zeros((N,n),dtype=np.int64)
    indices=np.zeros(N,dtype=np.int64);incoming=np.zeros(n,dtype=np.int64)
    for v in range(n):
        indices[1<<v]=v
        for u in range(n):
            if rows[u]>>v&1:incoming[v]|=1<<u
        e[1<<v,v]=1;s[1<<v,v]=1
    for mask in range(1,N):
        if not mask&(mask-1):continue
        todo=mask
        while todo:
            bit=todo&-todo;todo-=bit;v=indices[bit];rest=mask^bit
            pred=rest&incoming[v]
            while pred:
                bb=pred&-pred;pred-=bb;e[mask,v]+=e[rest,indices[bb]]
            succ=rest&rows[v]
            while succ:
                bb=succ&-succ;succ-=bb;s[mask,v]+=s[rest,indices[bb]]
    H=0
    for v in range(n):H+=e[full,v]
    out=np.empty(len(selections),dtype=np.int64)
    for j,O in enumerate(selections):
        tax=0
        for mask in range(1,full):
            a=0;b=0;todo=mask&O
            while todo:
                bb=todo&-todo;todo-=bb;a+=e[mask,indices[bb]]
            if not a:continue
            rest=full^mask;todo=rest&O
            while todo:
                bb=todo&-todo;todo-=bb;b+=s[rest,indices[bb]]
            tax+=a*b
        out[j]=H-tax
    return out

def certify(rows,marks):
    n=len(rows);stable={m:[] for m in marks};cert={m:{} for m in marks};false=[];calls=0
    for v in range(n):
        active=[m for m in marks if not m>>v&1]
        R=1+min(rows[v].bit_count(),n-1-rows[v].bit_count());values={m:[1] for m in active}
        for size in range(2,R+1):
            if not active:break
            labels,rr=ib.inflate(rows,v,size)
            sm=[sum(1<<i for i,q in enumerate(labels) if m>>q&1) for m in active]
            vals=selected_defects(np.array(rr,dtype=np.int64),np.array(sm,dtype=np.int64));calls+=1;remaining=[]
            for m,d in zip(active,vals):
                d=int(d);values[m].append(d)
                if d==1:remaining.append(m)
                elif size>2:false.append({'selected':m,'port':v,'R':R,'first_failed_length':size,'values':values[m]})
            active=remaining
        for m in active:stable[m].append(v);cert[m][v]={'R':R,'values':values[m]}
    return stable,cert,false,calls

def run(n):
    need(n in [9,10],'declared complete input orders')
    stem='nine' if n==9 else 'ten';source=read('unit_'+stem+'.json.gz');out=[];hist=Counter();false=[];calls=0;start=time.time();controls=0
    for i,r in enumerate(source['records']):
        rows=r['rows'];marks={m['selected']:m['U'] for m in r['marks']}
        P,C,F,num=certify(rows,list(marks));calls+=num
        for entry in F:false.append({'code':r['code'],'rows':rows,**entry})
        for m,ports in P.items():
            hist[len(ports)]+=1
            if not ports:continue
            _,_,comps=quotient(rows,m);eligible=[c for c in comps if len(c)>1 and set(c)&set(ports)]
            out.append({'code':r['code'],'rows':rows,'selected':m,'H':r['H'],'U':marks[m],'stable_ports':ports,'cyclic_triples':triangles(rows,ports),'eligible_blocks':eligible,'axis_certificates':C[m]})
        if i<20 or i%10000==0:
            h,u,w,_=tables(np.array(rows,dtype=np.int64));vals=selected_defects(np.array(rows,dtype=np.int64),np.array(list(marks),dtype=np.int64));need(all(int(x)==1 for x in vals),'full table positive control')
            allmarks=[0,(1<<n)-1,1,3];vals=selected_defects(np.array(rows,dtype=np.int64),np.array(allmarks,dtype=np.int64))
            for m,d in zip(allmarks,vals):need(int(d)==defect(h,w,[j for j in range(n) if m>>j&1]),'full table hostile control')
            controls+=len(allmarks)+len(marks)
        if i%10000==0:print(json.dumps({'order':n,'at':i,'total':len(source['records']),'stable_states':len(out),'calls':calls,'seconds':round(time.time()-start,1)}),flush=True)
    if n==9:
        old=read('stable_nine.json.gz');need(hist==Counter({int(k):v for k,v in old['stable_port_histogram'].items()}),'complete independent histogram')
        need([(r['code'],r['selected'],r['stable_ports']) for r in out]==[(r['code'],r['selected'],r['stable_ports']) for r in old['records']],'complete independent stable-state equality')
    result={'scope':f'Complete ordinary strong marked unit seed stable-port census at order {n}. Each axis certified through the selected-endpoint degree bound. Selected-boundary engine independently validated against full W tables; no converse merging.','input_sha256':hashlib.sha256((HERE/('unit_'+stem+'.json.gz')).read_bytes()).hexdigest(),'stable_states':len(out),'births':sum(not r['eligible_blocks'] for r in out),'descendants':sum(bool(r['eligible_blocks']) for r in out),'stable_port_histogram':dict(hist),'cyclic_stable_sets':sum(bool(r['cyclic_triples']) for r in out),'false_two_clone_screens':false,'literal_calls':calls,'full_table_controls':controls,'records':out}
    save('stable_'+stem+'_selected.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k not in ['records','false_two_clone_screens']}),flush=True)
if __name__=='__main__':run(int(sys.argv[1]))
