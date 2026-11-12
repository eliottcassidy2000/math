#!/usr/bin/env python3
"""Exact controls for marked batch extension normal forms. No external packages."""
from itertools import combinations, permutations, product
from collections import Counter
from fractions import Fraction
from math import comb, factorial
import json, random

def need(c,m):
    if not c: raise AssertionError(m)
def rows_from_code(n,code):
    rows=[0]*n
    for t,(i,j) in enumerate(combinations(range(n),2)):
        a,b=(i,j) if (code>>t)&1 else (j,i)
        rows[a]|=1<<b
    return tuple(rows)
def arc(t,i,j):return (t[i]>>j)&1
def chains(t):
    n=len(t);out={};inc={}
    for i,j in combinations(range(n),2):
        if ((t[i]^t[j]) & ~((1<<i)|(1<<j)))==0:
            a,b=(i,j) if arc(t,i,j) else (j,i)
            need(a not in out and b not in inc,'path degrees')
            out[a]=b;inc[b]=a
    ans=[]
    for v in range(n):
        if v not in inc:
            p=[v]
            while p[-1] in out:p.append(out[p[-1]])
            ans.append(p)
    need(sum(map(len,ans))==n,'path forest')
    return ans
def induced(t,vs):return tuple(sum(arc(t,v,w)<<j for j,w in enumerate(vs)) for v in vs)
def skeleton(t,n):
    old=induced(t,list(range(n)));ch=chains(old);groups=[];words=[]
    for p in ch:
        runs=[];last=None
        for v in p:
            s=tuple(arc(t,x,v) for x in range(n,len(t)))
            if s!=last:runs.append([v]);last=s
            else:runs[-1].append(v)
        groups+=runs;words.append([tuple(arc(t,x,g[0]) for x in range(n,len(t))) for g in runs])
    groups += [[x] for x in range(n,len(t))]
    sk=induced(t,[g[0] for g in groups])
    return sk,groups,ch,words

def check_normal(t,n):
    sk,groups,ch,words=skeleton(t,n)
    for i,g in enumerate(groups):
        for a,b in combinations(range(len(g)),2):need(arc(t,g[a],g[b])==1,'run is ordered')
        for j,h in enumerate(groups):
            if i!=j:
                for v in g:
                    for w in h:need(arc(t,v,w)==arc(sk,i,j),'exact expansion')
    need(len(chains(t))==len(chains(sk)),'core rank invariant')
    for cp in chains(t):
        origins={i for i,p in enumerate(ch) if any(v in p for v in cp)}
        need(len(origins)<=1,'cannot bridge distinct parent chains')
    return len(chains(t))

def ensemble(n,k):
    N=n+k;hist=Counter();oldsurv=Counter();E=Fraction(0);cases=0
    for code in range(1<<comb(N,2)):
        t=rows_from_code(N,code);old=induced(t,list(range(n)));b=len(chains(old));bp=check_normal(t,n)
        oldedges=[(p[i],p[i+1]) for p in chains(old) for i in range(len(p)-1)]
        surv=sum(all(arc(t,x,u)==arc(t,x,v) for x in range(n,N)) for u,v in oldedges)
        hist[b,bp]+=1;oldsurv[b,surv]+=1;cases+=1
    byb=Counter()
    for (b,bp),c in hist.items():byb[b]+=c
    for b,c in byb.items():
        mean=sum(Fraction(bp*v,c) for (bb,bp),v in hist.items() if bb==b)
        expected=Fraction(N)-Fraction(n-b,2**k)-Fraction(n*k+comb(k,2),2**(N-2))
        need(mean==expected,'exact batch mean')
        m=n-b
        for s in range(m+1):
            prob=Fraction(comb(m,s)*(2**k-1)**(m-s),2**(k*m))
            need(Fraction(oldsurv[b,s],c)==prob,'binomial old survival')
        # Coupling TV upper bound checked from distributions.
        supp=set(bp for bb,bp in hist if bb==b)|set(range(N-m,N+1))
        tv=sum(abs(Fraction(hist[b,bp],c)-(Fraction(comb(m,N-bp)*(2**k-1)**(m-N+bp),2**(k*m)) if 0<=N-bp<=m else 0)) for bp in supp)/2
        need(tv<=Fraction(n*k+comb(k,2),2**(N-2)),'TV bound')
    return {'n':n,'k':k,'labelled_children':cases,'parent_core_ranks':sorted(byb)}

def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def make_template(Q,R,words):
    r=sum(map(len,words));k=len(R); sk=[0]*(r+k); orig=[];sig=[]
    for i,w in enumerate(words):
        for x in w:orig.append(i);sig.append(x)
    for a,b in combinations(range(r+k),2):
        if b<r:v=1 if orig[a]==orig[b] else arc(Q,orig[a],orig[b])
        elif a<r:v=1-sig[a][b-r]
        else:v=arc(R,a-r,b-r)
        u,vv=(a,b) if v else (b,a);sk[u]|=1<<vv
    return tuple(sk)
def autos(t):
    return [p for p in permutations(range(len(t))) if all(arc(t,i,j)==arc(t,p[i],p[j]) for i,j in combinations(range(len(t)),2))]
def template_group(Q,R,words):
    offsets=[];a=0
    for w in words:offsets.append(a);a+=len(w)
    G=[]
    for p in autos(Q):
        for q in autos(R):
            good=True
            for i,w in enumerate(words):
                mapped=[tuple(s[q.index(j)] for j in range(len(R))) for s in w]
                if mapped!=words[p[i]]:good=False;break
            if good:G.append(tuple(offsets[p[i]]+j for i,w in enumerate(words) for j in range(len(w))))
    return G

def compositions(n,r):
    if r==1:
        if n>=1:yield(n,)
        return
    for a in range(1,n-r+2):
        for rest in compositions(n-a,r-1):yield(a,)+rest

def cycles(p):
    unseen=set(range(len(p)));out=[]
    while unseen:
        x=min(unseen);c=[]
        while x in unseen:c.append(x);unseen.remove(x);x=p[x]
        out.append(c)
    return out

def expand(sk,weights):
    n=sum(weights);groups=[];a=0
    for w in weights:groups.append(list(range(a,a+w)));a+=w
    rows=[0]*n
    for i,g in enumerate(groups):
        for a,b in combinations(g,2):rows[a]|=1<<b
        for j,h in enumerate(groups):
            if i!=j and arc(sk,i,j):
                for a in g:
                    for b in h:rows[a]|=1<<b
    return tuple(rows)
def marked_code(t,n):
    best=None;N=len(t)
    for a in permutations(range(n)):
        for b in permutations(range(n,N)):
            p=a+b;c=tuple(arc(t,p[i],p[j]) for i,j in combinations(range(N),2))
            if best is None or c<best:best=c
    return best

def family_control(Q,R,words,max_old):
    sk=make_template(Q,R,words);G=template_group(Q,R,words);r=sum(map(len,words));k=len(R);coeff={}
    for n in range(r,max_old+1):
        vectors=list(compositions(n,r));orbits=set();graphs=set()
        for w in vectors:
            orbits.add(min(tuple(w[p[i]] for i in range(r)) for p in G))
            t=expand(sk,w+(1,)*k);need(len(chains(t))==len(chains(sk)),'weighted rank')
            graphs.add(marked_code(t,n))
        fixed=sum(sum(all(w[i]==w[p[i]] for i in range(r)) for w in vectors) for p in G)
        need(fixed%len(G)==0,'Burnside integer')
        need(len(orbits)==len(graphs)==fixed//len(G),'marked coefficient bijection')
        coeff[n+k]=len(orbits)
    return {'runs':r,'new':k,'group_order':len(G),'kernel_possible':len(set(G))<len(G),'coefficients':coeff}

def ham(t):
    n=len(t);dp={}
    for v in range(n):dp[1<<v,v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            a=dp.get((mask,v),0)
            if a:
                for w in range(n):
                    if not (mask>>w)&1 and arc(t,v,w):dp[mask|(1<<w),w]=dp.get((mask|(1<<w),w),0)+a
    return sum(dp.get(((1<<n)-1,v),0) for v in range(n))
def stirling(n,k):
    a=[0]*(n+1);a[0]=1
    for j in range(1,n+1):
        b=[0]*(n+1)
        for r in range(1,j+1):b[r]=a[r-1]+r*a[r]
        a=b
    return a[k]
def content_ham(q,weights):
    r=len(q);dp={};ans=0
    for i in range(r):
        c=tuple(int(j==i) for j in range(r));dp[c,i]=1
    for total in range(1,sum(weights)+1):
        nxt={}
        for (c,last),a in dp.items():
            if all(c):
                w=a
                for n,k in zip(weights,c):w*=factorial(k)*stirling(n,k)
                ans+=w
            for j in range(r):
                if arc(q,last,j) and c[j]<weights[j]:
                    d=tuple(v+int(i==j) for i,v in enumerate(c));nxt[d,j]=nxt.get((d,j),0)+a
        dp=nxt
    return ans

def main():
    exact=[ensemble(n,k) for n,k in [(1,1),(2,1),(3,1),(4,1),(1,4),(2,3),(3,2),(4,2),(3,3)]]
    rng=random.Random(284)
    for N in range(7,13):
        for _ in range(30):
            t=rows_from_code(N,rng.getrandbits(comb(N,2)));check_normal(t,rng.randrange(1,N))
    c3=(2,4,1)
    fam=[family_control(c3,(0,),[[(0,)],[(0,)],[(0,)]],6),family_control((0,),(0,),[[(0,),(1,),(0,)]],6),family_control(c3,c3,[[(1,0,0)],[(0,1,0)],[(0,0,1)]],4)]
    # External full-signature births, with Gray ordered old witnesses.
    births=[]
    for k in range(1,9):
        old=1<<k;words=[[tuple(((i^(i>>1))>>j)&1 for j in range(k)) for i in range(old)]]
        sk=make_template((0,),transitive(k),words)
        need(len(set(words[0]))==old,'all external signatures')
        need(all(sum(a!=b for a,b in zip(x,y))==1 for x,y in zip(words[0],words[0][1:])),'Gray one-change')
        need(len(sk)==old+k,'sharp birth')
        births.append({'boundary':k,'first_order':old+k,'core_rank':len(chains(sk))})
    path_controls=0
    for q in [(2,4,1),transitive(3),make_template((0,),(0,),[[(0,),(1,),(0,)]])]:
        for w in product((1,2),repeat=len(q)):
            need(content_ham(q,w)==ham(expand(q,w)),'walk/Stirling Hamiltonian semantics')
            path_controls+=1
    # Hostile: old runs can merge through a newborn.
    t=(6,0,2) # old0 -> new2 -> old1, old0 -> old1
    sk,g,_,_=skeleton(t,2)
    need(len(g)==3 and len(chains(t))==1,'newborn bridges old runs')
    return {'status':'PASS','exact_ensembles':exact,'random_normal_forms':180,'hamiltonian_substitution_controls':path_controls,'families':fam,'external_births':births,'hostile':'two old runs joined by newborn middle vertex','scope':'labelled-uniform ensemble; marked batch isomorphisms, not unmarked A000568'}
if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
