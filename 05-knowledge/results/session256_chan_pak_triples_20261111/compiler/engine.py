"""Positive two-flank compiler and terminating modular Pareto optimizer.
No packages, no network. Run this file to reproduce results.json.
"""
from collections import deque
from functools import lru_cache
from itertools import permutations, product
import json
from pathlib import Path
CORES=((0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11))
SEED=(1,1,0,0,0,0,0,0,1)

def step(v,s,t):
    ga,gb,a,b,ab,p,q,pq,u=v
    return (ga+(1-s)*(a+p),gb+(1-t)*(b+q),2*a+s*u,
            2*b+t*u,4*ab+2*t*a+2*s*b+s*t*u,s*u,t*u,s*t*u,u)

def matrix(s,t):
    return tuple(tuple(step(tuple(int(i==j) for i in range(9)),s,t)[k]
                       for j in range(9)) for k in range(9))

@lru_cache(None)
def g(n):
    if n<2:return 1
    q=n//2
    return g(q)+(q+(q&1) if not n&1 else 0)

def state(a,b):return (g(a),g(b),a,b,a*b,a&1,b&1,(a&1)*(b&1),1)

def H(rows):
    n=len(rows)
    if not n:return 1
    dp={ (1<<i,i):1 for i in range(n)}
    for mask in range(1,1<<n):
        for i in range(n):
            v=dp.get((mask,i),0)
            if not v:continue
            nxt=rows[i]&~mask
            while nxt:
                bit=nxt&-nxt;nxt-=bit;j=bit.bit_length()-1
                dp[mask|bit,j]=dp.get((mask|bit,j),0)+v
    return sum(dp.get(((1<<n)-1,i),0) for i in range(n))

def induced(rows,mask):
    ids=[i for i in range(len(rows)) if mask>>i&1]
    return tuple(sum((1<<j) for j,k in enumerate(ids) if rows[i]>>k&1) for i in ids)

def signature(rows,cut):
    n=len(rows);full=(1<<n)-1;out=[0]*8
    for l in range(1,n+1):
        for path in permutations(range(n),l):
            if not all(rows[a]>>b&1 for a,b in zip(path,path[1:])):continue
            w=H(induced(rows,full^sum(1<<i for i in path)))
            flags=(True,not(cut>>path[-1]&1),bool(cut>>path[0]&1),
                   bool(cut>>path[0]&1) and not(cut>>path[-1]&1))
            for i,yes in enumerate(flags):
                if yes:out[2*i]+=w;out[2*i+1]+=w*((-1)**l)
    return (H(rows),*out)

def value(sig,v):
    d,kp,km,ip,im,op,om,jp,jm=sig
    ga,gb,a,b,ab,p,q,pq,u=v
    return d*(ga+gb-u)+(d+kp)*ab+(d+km)*pq+ip*a-im*p+op*b-om*q+(jp+jm)*u

def graph(rows,cut,a,b):
    left=[] if not a else list(map(int,bin(a)[2:]))
    right=[] if not b else [1-int(x) for x in bin(b)[2:][::-1]]
    na,nc,nb=len(left),len(rows),len(right);n=na+nc+nb
    ans=[0]*(n+1)
    for i in range(n):
        for j in range(i+1,n):
            win=(rows[i-na]>>(j-na)&1) if na<=i<j<na+nc else 1
            ans[i if win else j]|=1<<(j if win else i)
    bits=left+[(cut>>i)&1 for i in range(nc)]+right
    for i,z in enumerate(bits):ans[n if z else i]|=1<<(i if z else n)
    return tuple(ans)

def dominated(v,w):return all(a<=b for a,b in zip(v,w))

def optimize(modulus, controller=None):
    """Return complete nondominated reachable states, with address words.
    Controller = (initial, transition(state,s,t), accepting(state)); default all.
    Modulus must be positive even, preserving exact Boolean coordinates.
    """
    if modulus<2 or modulus%2:raise ValueError('use a positive even modulus')
    if controller is None:controller=(0,lambda c,s,t:0,lambda c:True)
    start,transition,accept=controller
    def key(c,v):return c,tuple(x%modulus for x in v)
    bank={key(start,SEED):{SEED:''}};queue=deque([(start,SEED)])
    additions=1;max_live=1
    while queue:
        c,v=queue.popleft();k=key(c,v)
        if v not in bank.get(k,{}):continue
        word=bank[k][v]
        for s,t in product(range(2),repeat=2):
            cc=transition(c,s,t)
            if cc is None:continue
            w=step(v,s,t);kk=key(cc,w);bucket=bank.setdefault(kk,{})
            if any(dominated(old,w) for old in bucket):continue
            for old in list(bucket):
                if dominated(w,old):del bucket[old]
            bucket[w]=word+str(2*s+t);queue.append((cc,w));additions+=1
        max_live=max(max_live,sum(map(len,bank.values())))
    accepted=[(c,v,w) for (c,res),items in bank.items() if accept(c) for v,w in items.items()]
    return accepted,{'modulus':modulus,'buckets':len(bank),'retained':sum(map(len,bank.values())),
                    'accepted':len(accepted),'additions':additions,'max_live':max_live,
                    'max_word_length':max(map(lambda x:len(x[2]),accepted),default=0)}

def run():
    matrices={str(2*s+t):matrix(s,t) for s,t in product(range(2),repeat=2)}
    controls=0
    for a,b,s,t in product(range(32),range(32),range(2),range(2)):
        assert step(state(a,b),s,t)==state(2*a+s,2*b+t);controls+=1
    sigs=[(ci,cut,signature(rows,cut)) for ci,rows in enumerate(CORES) for cut in range(1<<len(rows))]
    direct=0
    for ci,cut,sig in sigs:
        for a,b in product(range(4),repeat=2):
            assert H(graph(CORES[ci],cut,a,b))==value(sig,state(a,b));direct+=1
    experiments=[]
    for mod in (2,4,6,8):
        accepted,stats=optimize(mod)
        best={}
        for ci,cut,sig in sigs:
            for c,v,word in accepted:
                h=value(sig,v);r=h%mod
                rec=(h,len(word),ci,cut,v[2],v[3],word)
                if r not in best or rec<best[r]:best[r]=rec
        records=[]
        for r,rec in sorted(best.items()):
            h,_,ci,cut,a,b,word=rec;rows=graph(CORES[ci],cut,a,b)
            assert H(rows)==h
            records.append({'residue':r,'minimum_H':h,'core':ci,'cut':cut,'A':a,'B':b,'word':word,'rows':rows})
        experiments.append({'stats':stats,'minimum_by_residue':records})
    # Regular driver: both flanks nonzero, and each address has an even number of 1 bits.
    def tr(c,s,t):
        seen_a,seen_b,pa,pb=c
        return (seen_a or s,seen_b or t,pa^s,pb^t)
    controller=((0,0,0,0),tr,lambda c:c==(1,1,0,0))
    accepted,stats=optimize(4,controller)
    best={}
    for ci,cut,sig in sigs:
        for c,v,word in accepted:
            h=value(sig,v);r=h%4;rec=(h,ci,cut,v[2],v[3],word)
            if r not in best or rec<best[r]:best[r]=rec
    regular=[]
    for r,(h,ci,cut,a,b,word) in sorted(best.items()):
        rows=graph(CORES[ci],cut,a,b);assert H(rows)==h
        regular.append({'residue':r,'minimum_H':h,'core':ci,'cut':cut,'A':a,'B':b,'word':word,'rows':rows})
    # Synchronized primitive-Pythagorean parameter controller modulo 5.
    # Digits 0,1,2 mean bit pairs 00,01,10. First digit nonzero;
    # require both addresses positive and parameter residues (1,0).
    def ppt_tr(c,s,t):
        m,n,started,sa,sb=c
        digit=2*s+t
        if digit==3 or (not started and digit==0):return None
        mm,nn=((2*m-n,m),(2*m+n,m),(m+2*n,n))[digit]
        return mm%5,nn%5,1,sa or s,sb or t
    ppt_ctl=((2,1,0,0,0),ppt_tr,lambda c:c==(1,0,1,1,1))
    accepted,ppt_stats=optimize(4,ppt_ctl)
    ppt_best={}
    for ci,cut,sig in sigs:
        for c,v,word in accepted:
            h=value(sig,v);r=h%4;rec=(h,ci,cut,v[2],v[3],word)
            if r not in ppt_best or rec<ppt_best[r]:ppt_best[r]=rec
    ppt_records=[]
    for r,(h,ci,cut,a,b,word) in sorted(ppt_best.items()):
        m,n=2,1
        for digit in map(int,word):m,n=((2*m-n,m),(2*m+n,m),(m+2*n,n))[digit]
        triple=(m*m-n*n,2*m*n,m*m+n*n)
        rows=graph(CORES[ci],cut,a,b);assert H(rows)==h
        assert (m%5,n%5)==(1,0) and triple[0]**2+triple[1]**2==triple[2]**2
        ppt_records.append({'residue':r,'minimum_H':h,'core':ci,'cut':cut,'A':a,'B':b,
                            'word':word,'parameters':[m,n],'triple':triple,'rows':rows})
    result={'digit_controls':controls,'direct_graph_controls':direct,'matrices':matrices,
            'pythagorean_driver':{'stats':ppt_stats,'minimum_by_residue':ppt_records},
            'experiments':experiments,'regular_driver':{'stats':stats,'minimum_by_residue':regular}}
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='matrices'},indent=2))
if __name__=='__main__':run()
