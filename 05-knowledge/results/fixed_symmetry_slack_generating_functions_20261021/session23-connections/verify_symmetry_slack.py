#!/usr/bin/env python3
"""Independent finite controls for maximal-W compiler. Exact exhaustive unlabelled n<=6."""
from itertools import permutations, product
from functools import lru_cache
from collections import Counter, defaultdict
import json,sys
MAX_N=int(sys.argv[1]) if len(sys.argv)>1 else 6
if MAX_N not in (6,7): raise ValueError("Choose exhaustive maximum order6 or7")

CHECKS=Counter()
def ck(x,name):
    CHECKS[name]+=1
    if not x: raise RuntimeError(name)

def matrix(n,mask):
    a=[[0]*n for _ in range(n)]
    k=0
    for i in range(n):
        for j in range(i+1,n):
            a[i][j]=(mask>>k)&1;a[j][i]=1-a[i][j];k+=1
    return tuple(map(tuple,a))

def encode(a,p):
    z=0;k=0
    for i in range(len(p)):
        for j in range(i+1,len(p)):
            z|=a[p[i]][p[j]]<<k;k+=1
    return z

@lru_cache(None)
def canonical(a):
    n=len(a);by=defaultdict(list)
    for i in range(n):by[sum(a[i])].append(i)
    groups=[by[d] for d in sorted(by)]
    return min(encode(a,sum((tuple(b) for b in bs),())) for bs in product(*(permutations(g) for g in groups))) if n else 0

@lru_cache(None)
def autos(a):
    n=len(a);degs=list(map(sum,a));gs=defaultdict(list)
    for i,d in enumerate(degs):gs[d].append(i)
    groups=list(gs.values());ans=[]
    for bs in product(*(permutations(g) for g in groups)):
        p=list(range(n))
        for old,new in zip(groups,bs):
            for i,j in zip(old,new):p[i]=j
        if all(a[i][j]==a[p[i]][p[j]] for i in range(n) for j in range(i+1,n)):
            ans.append(tuple(p))
    return ans

@lru_cache(None)
def anti_maps(a):
    n=len(a);ds=list(map(sum,a));by=defaultdict(list)
    for i,d in enumerate(ds):by[d].append(i)
    pairs=[(vs,by.get(n-1-d,[])) for d,vs in by.items()]
    if any(len(x)!=len(y) for x,y in pairs):return []
    out=[]
    for bs in product(*(permutations(ys) for xs,ys in pairs)):
        p=list(range(n))
        for (xs,ys),zs in zip(pairs,bs):
            for x,z in zip(xs,zs):p[x]=z
        if all(a[p[i]][p[j]]==1-a[i][j] for i in range(n) for j in range(i+1,n)):out.append(tuple(p))
    return out

def antis(a):
    n=len(a);b=tuple(tuple(0 if i==j else 1-a[i][j] for j in range(n)) for i in range(n))
    return canonical(a)==canonical(b)

@lru_cache(None)
def feedback(a):
    n=len(a);dp=[0]*(1<<n)
    for s in range(1,1<<n):
        dp[s]=min(dp[s^(1<<v)]+sum(a[v][u] for u in range(n) if s>>u&1 and u!=v) for v in range(n) if s>>v&1)
    return dp[-1]

def induced(a,ss):return tuple(tuple(a[i][j] for j in ss) for i in ss)

def sccs(a,ss):
    ss=list(ss);reach={i:{i}|{j for j in ss if a[i][j]} for i in ss}
    for k in ss:
        for i in ss:
            if k in reach[i]:reach[i]|=reach[k]
    comps=[];seen=set()
    for i in ss:
        if i not in seen:
            c={j for j in ss if j in reach[i] and i in reach[j]};seen|=c;comps.append(c)
    # score among components determines condensation order
    rank={min(c):sum(a[next(iter(c))][next(iter(d))] for d in comps if c!=d) for c in comps}
    comps.sort(key=lambda c:-rank[min(c)])
    return comps

def word(a,ss):
    cs=sccs(a,ss)
    return tuple(len(c) for c in cs) if all(len(c) in (1,3) for c in cs) else None

def module(a,ss):
    ss=set(ss)
    return all(len({a[v][u] for u in ss})<=1 for v in range(len(a)) if v not in ss)

def maximal_w(a):
    n=len(a);mods=[]
    for s in range(1,1<<n):
        ss={i for i in range(n) if s>>i&1}
        if module(a,ss) and word(a,ss) is not None:mods.append(ss)
    for A in mods:
        for B in mods:
            if A&B:
                ck(module(a,A|B) and word(a,A|B) is not None,'overlap_union')
    ms=[s for s in mods if not any(s<t for t in mods)]
    ck(set().union(*ms)==set(range(n)) and sum(map(len,ms))==n,'partition')
    return sorted(ms,key=lambda s:min(s))

def canonical_test(q,flags):
    r=len(q)
    for mask in range(1,1<<r):
        I=[i for i in range(r) if mask>>i&1]
        if len(I)<2 or not module(q,I):continue
        cs=sccs(q,I)
        if all(len(c)==1 or (len(c)==3 and all(flags[i] for i in c)) for c in cs):return False
    return True

def inflate(q,ws):
    # standardized directed triangles and ordinal words
    blocks=[];off=0;a=[]
    sizes=list(map(sum,ws));N=sum(sizes);a=[[0]*N for _ in range(N)]
    for w in ws:
        ids=list(range(off,off+sum(w)));blocks.append(ids)
        pos=off
        for letter in w:
            for i in range(pos,pos+letter):
                for j in range(i+1,off+sum(w)):a[i][j]=1;a[j][i]=0
            if letter==3:a[pos][pos+2]=0;a[pos+2][pos]=1
            pos+=letter
        off+=sum(w)
    for i in range(len(ws)):
        for j in range(i+1,len(ws)):
            for u in blocks[i]:
                for v in blocks[j]:a[u][v]=q[i][j];a[v][u]=q[j][i]
    return tuple(map(tuple,a))

def allwords_upto(N,nonempty=False):
    ans=[]
    def go(w,n):
        if w or not nonempty:ans.append(w)
        for k in (1,3):
            if n+k<=N:go(w+(k,),n+k)
    go((),0);return ans

def analyze(a):
    n=len(a);G=autos(a);f=feedback(a);s=sum(any(g[v]!=v for g in G) for v in range(n));d=3*f-s
    ck(d>=0,'nonnegative_slack');blocks=maximal_w(a);r=len(blocks)
    q=tuple(tuple(0 if i==j else a[min(blocks[i])][min(blocks[j])] for j in range(r)) for i in range(r))
    ws=tuple(word(a,b) for b in blocks);ns=list(map(sum,ws));ts=[w.count(3) for w in ws];us=[w.count(1) for w in ws]
    ck(canonical_test(q,[w==(1,) for w in ws]),'finite_maximality')
    ck(canonical(inflate(q,ws))==canonical(a),'reconstruction')
    K=[g for g in autos(q) if all(ws[i]==ws[g[i]] for i in range(r))]
    J=[i for i in range(r) if any(g[i]!=i for g in K)]
    U=sum(us[i] for i in J);V=sum(ns[i] for i in J)
    costs=[]
    for p in permutations(range(r)):
        pos={v:i for i,v in enumerate(p)}
        costs.append((sum(ns[i]*ns[j] for i in range(r) for j in range(r) if q[i][j] and pos[i]>pos[j]),p))
    C=min(c for c,p in costs)
    ck(f==sum(ts)+C,'weighted_straightening')
    ck(s==3*sum(ts)+U,'support_formula')
    ck(d==3*C-U,'slack_formula');ck(V<=2*d,'moved_mass_bound');ck(C<=d,'core_feedback_bound')
    if d==0:ck(r==1 and word(a,range(n)) is not None,'zero_slack');return d,f,antis(a),None
    ck(r<=4*d+1,'core_order_bound')
    L={i for i in range(r) if ns[i]>d};ck(len(L)<=2*d+1,'long_slot_bound')
    admissible=[]
    for c,p in costs:
        pos={v:i for i,v in enumerate(p)}
        if all(pos[i]<pos[j] for i in range(r) for j in range(r) if q[i][j] and (i in L or j in L)):admissible.append((c,p))
    ck(bool(admissible) and min(c for c,p in admissible)==C,'admissible_cost')
    K0=[g for g in autos(q) if {g[i] for i in L}==L and all(ws[i]==ws[g[i]] for i in range(r) if i not in L)]
    ck(all(g[i]==i for g in K0 for i in L),'long_fixed');ck(K==K0,'exact_stabilizer')
    AC=[g for g in anti_maps(q) if {g[i] for i in L}==L and all(ws[g[i]]==ws[i][::-1] for i in range(r) if i not in L)]
    ck(antis(a)==any(all(ws[g[i]]==ws[i][::-1] for i in L) for g in AC),'converse_configuration')
    if AC:
        ck(all(tuple(g[i] for i in sorted(L))==tuple(AC[0][i] for i in sorted(L)) for g in AC),'unique_long_anti_action')
        ck(sum(AC[0][i]==i for i in L)==len(L)%2 and all(AC[0][AC[0][i]]==i for i in L),'long_anti_pairs_middle')

    # Change every long word independently, retain the full finite configuration.
    replacements=[w for w in allwords_upto(min(7,n+1),True) if sum(w)>d]
    tested=0
    for i in L:
        for w in replacements:
            nw=list(ws);nw[i]=w;b=inflate(q,nw)
            if len(b)>7:continue
            Gb=autos(b);fb=feedback(b);sb=sum(any(g[v]!=v for g in Gb) for v in range(len(b)))
            ck(3*fb-sb==d,'long_fill_slack')
            ck(len(maximal_w(b))==r,'long_fill_maximality')
            tested+=1
    return d,f,antis(a),(q,ws,L,C)

# Generate isomorphism classes by all one-vertex extensions; canonical score classes are exact.
classes={0:{0}};counts=defaultdict(Counter);sc=defaultdict(Counter);joint=defaultdict(Counter);scjoint=defaultdict(Counter);examples={}
for n in range(1,MAX_N+1):
    masks=set()
    for m in classes[n-1]:
        a0=matrix(n-1,m)
        for row in range(1<<(n-1)):
            a=[list(x)+[(row>>i)&1] for i,x in enumerate(a0)]
            a.append([1-((row>>i)&1) for i in range(n-1)]+[0])
            masks.add(canonical(tuple(map(tuple,a))))
    classes[n]=masks
    ck(len(masks)==[1,1,1,2,4,12,56,456][n],'class_count')
    for m in sorted(masks):
        a=matrix(n,m);d,f,isc,config=analyze(a);counts[d][n]+=1;sc[d][n]+=isc;joint[d][n,f]+=1;scjoint[d][n,f]+=isc;examples.setdefault(d,(n,m))
    print('n',n,'classes',len(masks),'slack',dict(sorted((d,c[n]) for d,c in counts.items() if c[n])),'self_converse',sum(c[n] for c in sc.values()),flush=True)
# Independent direct word-series controls for inherited exact layers.
words=allwords_upto(MAX_N)
for n in range(1,MAX_N+1):
    ck(counts[0][n]==sum(sum(w)==n for w in words),'delta0_series')
    e3=sum(sum(l)+2+sum(mid)+sum(r)==n for l in words for mid in words if sum(mid)>=2 for r in words)
    e4=sum(sum(l)+5+sum(r)==n for l in words for r in words)
    ck(counts[3][n]==e3,'delta3_series');ck(counts[4][n]==e4,'delta4_series')
    ck(sc[0][n]==sum(sum(w)==n and w==w[::-1] for w in words),'delta0_palindrome')
    for f in range(7):
        e3b=sum(sum(l)+2+sum(mid)+sum(r)==n and 1+l.count(3)+mid.count(3)+r.count(3)==f for l in words for mid in words if sum(mid)>=2 for r in words)
        e4b=sum(sum(l)+5+sum(r)==n and 3+l.count(3)+r.count(3)==f for l in words for r in words)
        s3b=sum(2*sum(l)+2+sum(mid)==n and 1+2*l.count(3)+mid.count(3)==f for l in words for mid in words if sum(mid)>=2 and mid==mid[::-1])
        s4b=sum(2*sum(l)+5==n and 3+2*l.count(3)==f for l in words)
        ck(joint[3][n,f]==e3b,'delta3_bivariate')
        ck(joint[4][n,f]==e4b,'delta4_bivariate')
        ck(scjoint[3][n,f]==s3b,'delta3_converse_bivariate')
        ck(scjoint[4][n,f]==s4b,'delta4_converse_bivariate')
# Actual phase-locked diagonal C3 hostile.
a=[[0]*6 for _ in range(6)]
for shift in (0,3):
    for i in range(3):a[shift+i][shift+(i+1)%3]=1
for i in range(3):
    for j in range(3):a[i][3+j]=int(i==j);a[3+j][i]=int(i!=j)
a=tuple(map(tuple,a));d,f,isc,config=analyze(a)
ck((f,d,len(autos(a)))==(4,6,3),'diagonal_c3_hostile')
print('diagonal_c3_maximal_W_blocks', [sorted(b) for b in maximal_w(a)], 'configuration',config)
ck(not module(a,{0,1,2}) and not module(a,{3,4,5}),'diagonal_orbits_not_modules')
for m in classes[6]:
    a=matrix(6,m);g=autos(a);f=feedback(a);s=sum(any(p[v]!=v for p in g) for v in range(6));d=3*f-s
    if d==6 and len(g)>1:
        print('delta6_symmetric_order6',m,'aut',len(g),'feedback',f,'scores',sorted(map(sum,a)),'maxW',[sorted(b) for b in maximal_w(a)],'self_converse',antis(a))
print('exact checks',sum(CHECKS.values()),dict(sorted(CHECKS.items())))
print('delta_histograms',json.dumps({d:dict(sorted(c.items())) for d,c in sorted(counts.items())},sort_keys=True))
print('self_converse_histograms',json.dumps({d:dict(sorted(c.items())) for d,c in sorted(sc.items())},sort_keys=True))
print('PASS')
