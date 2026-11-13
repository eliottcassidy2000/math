"""Exact scalar-tail inversion for physical unique-start funnel cores."""
from pathlib import Path
from functools import lru_cache
import sys,json,time,argparse,math
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'berggren_target_return_20261111'/'reset'))
import funnel_six_state as f

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def g(n):
    value=1
    while n:
        half=n>>1
        if not n&1:value+=half+(half&1)
        n=half
    return value

def inverse(target,A,B,C,D,leaf=4):
    """All n>=0 with A*n+B*g(n)+C*parity(n)+D=target; A,B>0,C>=0."""
    need(A>0 and B>0 and C>=0,'inverse cone')
    stack=[(A,C,D,0,1)];hits=[];nodes=leaves=0
    while stack:
        a,c,d,base,scale=stack.pop();nodes+=1
        lo=max(0,(target-B-c-d+a+B-1)//(a+B));hi=(target-B-d)//a
        if lo>hi:continue
        if hi-lo<=leaf:
            leaves+=hi-lo+1
            for n in range(lo,hi+1):
                if a*n+B*g(n)+c*(n&1)+d==target:hits.append(base+scale*n)
            continue
        stack.append((2*a,0,a+c+d,base+scale,2*scale))
        stack.append((2*a+B,B,d,base,2*scale))
    return sorted(hits),dict(nodes=nodes,candidates=leaves)

def address(n,rank):
    if n==1:return []
    m,j,w=f.decode(n,rank)
    return address(m,j)+[''.join(map(str,w))]

def cores(order):
    rows=[dict(order=1,rank=0,state=f.SINGLE,words=[])]
    seen={(1,1,3,1)}
    for n in range(3,order+1):
        for rank in range(f.a.class_counts(n)[n]):
            c=f.compile_rank(n,rank);h,e,s,p,ee,t=c;key=(h,e,s+p,ee)
            if key in seen:continue
            seen.add(key);rows.append(dict(order=n,rank=rank,state=c,words=address(n,rank)))
    return rows

def direct(T,start=None):
    n=len(T);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):
        if start is None or start==v:dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            val=dp[mask][v]
            if not val:continue
            rem=T[v]&~mask
            while rem:
                bit=rem&-rem;rem-=bit;w=bit.bit_length()-1
                dp[mask|bit][w]+=val
    return sum(dp[-1])

def literal_core(row):
    return (0,) if row['order']==1 else f.a.unrank(row['order'],row['rank'])

def physical(row,n):
    U=literal_core(row);word=tuple((n>>j)&1 for j in range(n.bit_length()))
    V,x=f.a.funnel(U,(0,),word);T,z=f.a.birth(V,x)
    return T,z,V,x,word

def verify_controls():
    checks=0
    for n in range(1,1024):
        w=tuple((n>>j)&1 for j in range(n.bit_length()))
        need(g(n)==f.a.tail_coefficients(w,1)[2],'g is actual B1');checks+=1
    for A,B,C,D in [(3,1,1,0),(15,3,1,10),(7,5,3,-11),(1,2,0,0),(4,3,5,17)]:
        for target in range(-5,501):
            expected=[n for n in range(max(0,(target-B-D)//A)+1) if A*n+B*g(n)+C*(n&1)+D==target]
            actual,_=inverse(target,A,B,C,D)
            need(actual==expected,('complete inverse',A,B,C,D,target));checks+=1
    return checks

def run(order,limit):
    cs=cores(order);out=[];totalnodes=trials=0;t0=time.time()
    for k in range(1,limit+1):
        target=405*(1<<k)+71;found=None;knodes=0;kcores=0
        for row in cs:
            h,e,s,p,ee,t=row['state'];G=s+p
            ns,meta=inverse(target,G,h,e,G-e-ee-1);knodes+=meta['nodes'];kcores+=1
            if row['order']==1:ns=[n for n in ns if n>0]
            if not ns:continue
            n=min(ns,key=lambda n:(n.bit_length(),n));word=tuple((n>>j)&1 for j in range(n.bit_length()))
            state=f.advance(row['state'],word)
            need(state[0]+state[4]-1==target,'six-state independent target')
            T,z,V,x,w=physical(row,n)
            rec=dict(k=k,H=target,core_order=row['order'],core_rank=row['rank'],core_words=row['words'],tail_address=n,tail_word=''.join(map(str,w)),parent_state=state,rows=T,root=z,predecessor=x,order=len(T),cores_tried=kcores,inverse_nodes=knodes)
            if len(T)<=18:
                need(direct(T)==target,'literal full count');need(direct(V,x)==1,'literal unique-start parent')
                keep=[i for i in range(len(T)) if i!=z];delT=f.a.induced(T,keep)
                need(direct(delT)-direct(T,z)==1,'literal defect1');rec['literal_checked']=True
            else:rec['literal_checked']=False
            found=rec;break
        totalnodes+=knodes;trials+=kcores
        if found is None:found=dict(k=k,H=target,scope=f'no tail over the declared core bank through order{order}',cores_tried=kcores,inverse_nodes=knodes)
        out.append(found)
        print('k',k,'H',target,'source',None if 'rows' not in found else (found['core_words'],found['tail_word']),'nodes',knodes,'seconds',round(time.time()-t0,1),flush=True)
    return dict(core_order=order,distinct_observers=len(cs),targets=out,total_inverse_nodes=totalnodes,core_trials=trials)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--core-order',type=int,default=10);p.add_argument('--k',type=int,default=30);p.add_argument('--controls-only',action='store_true');a=p.parse_args()
    checks=verify_controls();data=dict(status='PASS',exact_inverse_controls=checks)
    if not a.controls_only:data['search']=run(a.core_order,a.k)
    suffix='' if __debug__ else '.optimized'
    (HERE/('inverse_receipt'+suffix+'.json')).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS',checks,'inverse controls')
