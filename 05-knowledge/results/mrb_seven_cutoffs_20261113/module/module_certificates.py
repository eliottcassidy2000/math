"""Positive run support, three-size certificates, and bounded ready-family search."""
from pathlib import Path
from itertools import product,permutations
import sys,json,time,argparse
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'berggren_target_return_20261111'/'reset'))
import funnel_six_state as f

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def coefficients(T,v,start=None,limit=None):
    """Coefficient of ordered j nonempty module paths; optional external start."""
    fixed=[i for i in range(len(T)) if i!=v];m=len(fixed)
    J=min(m+1,limit) if limit else m+1
    full=(1<<m)-1;dp=[[[0]*(J+1) for _ in range(m+1)] for _ in range(full+1)]
    if start is None:dp[0][m][1]=1
    for i,u in enumerate(fixed):
        if start is None or start==u:dp[1<<i][i][0]=1
    for mask in range(full+1):
        for j in range(J+1):
            for last in range(m+1):
                value=dp[mask][last][j]
                if not value:continue
                u=v if last==m else fixed[last]
                if last!=m and j<J and T[u]>>v&1:dp[mask][m][j+1]+=value
                for i,w in enumerate(fixed):
                    if not mask>>i&1 and T[u]>>w&1:dp[mask|1<<i][i][j]+=value
    return tuple(sum(dp[full][last][j] for last in range(m+1)) for j in range(1,J+1))

def inflate(T,v,k):
    owner=[u for u in range(len(T)) for _ in range(k if u==v else 1)]
    rows=tuple(sum((int(i<j) if u==w else (T[u]>>w&1))<<j for j,w in enumerate(owner) if i!=j) for i,u in enumerate(owner))
    return rows,{u:owner.index(u) for u in range(len(T))}

def count(T,start=None):
    n=len(T);dp=[[0]*n for _ in range(1<<n)]
    for u in range(n):
        if start is None or u==start:dp[1<<u][u]=1
    for mask in range(1,1<<n):
        for last in range(n):
            val=dp[mask][last]
            if not val:continue
            rem=T[last]&~mask
            while rem:
                bit=rem&-rem;w=bit.bit_length()-1;rem-=bit
                dp[mask|bit][w]+=val
    return sum(dp[-1])

def kj(k,j):
    from math import comb
    return sum((-1)**(j-i)*comb(j,i)*i**k for i in range(j+1))

def controls(order):
    contexts=observations=literal=0;hostile=None;support={}
    for n in range(1,order+1):
        for bits in range(1<<((n-1)*(n-2)//2)):
            T=f.a.tiling(n,bits)
            for v in range(n):
                contexts+=1
                for start in [None]+[u for u in range(n) if u!=v]:
                    cc=coefficients(T,v,start);observations+=1
                    active=[i+1 for i,x in enumerate(cc) if x]
                    need(not active or active==list(range(1,max(active)+1)),('run support',T,v,start,cc))
                    support[str(max(active,default=0))]=support.get(str(max(active,default=0)),0)+1
                    vals=[]
                    for k in (1,2,3):
                        U,ix=inflate(T,v,k);h=count(U,None if start is None else ix[start]);literal+=1
                        need(h==sum(c*kj(k,j) for j,c in enumerate(cc,1)),('skeleton/literal',T,v,start,k))
                        vals.append(h)
                    E=vals[2]-3*vals[1]+2*vals[0]
                    need(E==6*(cc[2] if len(cc)>2 else 0),'third finite difference')
                    need((E==0)==(not any(cc[2:])),('three-size iff',T,v,start))
                    if E and hostile is None and start is None:hostile=dict(T=T,module=v,coefficients=cc,values=vals,E=E)
    # Independent literal permutations on the sharp five-vertex hostile.
    if hostile:
        T=hostile['T'];v=hostile['module']
        U,_=inflate(T,v,3)
        hp=sum(all(U[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(len(U))))
        need(hp==hostile['values'][2],'literal permutations hostile')
    return dict(scope=f'All Hamilton-frame presentations through order{order}, every module vertex, full H and every external starting vertex',contexts=contexts,observations=observations,literal_counts=literal,support_histogram=support,hostile=hostile)

def address(n,rank):
    if n==1:return []
    m,j,w=f.decode(n,rank)
    return address(m,j)+[''.join(map(str,w))]

def search(order):
    targets={97:(45,7),293:(135,23),881:(405,71),2645:(1215,215),7937:(3645,647)}
    banks={1:[f.SINGLE]};tested=0;matches=[];candidates=[];census=[];t0=time.time()
    for n in range(3,order+1):
        states=[]
        for m in [1]+list(range(3,n)):
            t=n-m-1
            words=[()] if t==0 else [tuple((j>>i)&1 for i in range(t-1))+(1,) for j in range(1<<(t-1))]
            for w in words:states.extend(f.advance(c,w) for c in banks[m])
        need(len(states)==f.a.class_counts(n)[n],'complete unique-start class rank order')
        banks[n]=states;here=[]
        for rank,c in enumerate(states):
            h=c[0]+c[4]-1
            if h not in targets:continue
            core=f.a.unrank(n,rank);T,z=f.a.birth(core,0);a,b=targets[h]
            need(count(T)==h,'direct ready target count')
            good=[];profiles=[]
            for v in range(1,n):
                cc=coefficients(T,v,limit=3);tested+=1
                profiles.append([v,*cc])
                if cc[:3]!=(h,a,0):continue
                ss=coefficients(core,v,start=0,limit=3)
                need(ss==(1,0,0),'ready start survives all-size certificate')
                good.append(v);matches.append(dict(target=h,A=a,B=b,core_order=n,rank=rank,words=address(n,rank),module=v,ready=T,root=z,pred=0,H_runs3=cc,start_runs3=ss))
            candidates.append(dict(target=h,order=n+1,rank=rank,words=address(n,rank),profiles3=profiles,matching_vertices=good))
            here.append(h)
        census.append(dict(core_order=n,classes=len(states),min_ready=min(c[0]+c[4]-1 for c in states),target_candidates=len(here)))
        print('order',n,'classes',len(states),'target_candidates',len(here),'module_tests',tested,'hits',[(q['target'],q['module']) for q in matches],'seconds',round(time.time()-t0,1),flush=True)
    return dict(scope=f'Every strong unique-start rooted core through order{order}, followed by one indegree-one root birth; every vertex except that root and its sole predecessor is the uniform TT_k module; target rays45/135/405/1215/3645*2^k+7/23/71/215/647',census=census,tested_modules=tested,matches=matches,candidates=candidates)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--controls',type=int,default=5);p.add_argument('--search',type=int,default=17);a=p.parse_args()
    result=dict(controls=controls(a.controls))
    if a.search:result['search']=search(a.search)
    (HERE/('verification.json' if __debug__ else 'verification.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS',result['controls']['observations'],'observations',flush=True)
