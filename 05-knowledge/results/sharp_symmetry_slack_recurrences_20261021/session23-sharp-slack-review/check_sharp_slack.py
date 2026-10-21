#!/usr/bin/env python3
"""Independent finite controls for sharp slack poles; explicit raises survive -O."""
from collections import Counter
from itertools import product
from functools import lru_cache
import json, random, hashlib

checks=Counter()
def ck(ok,name):
    checks[name]+=1
    if not ok: raise RuntimeError(name)

def graph(n,mask):
    rows=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            if (mask>>k)&1: rows[i]|=1<<j
            else: rows[j]|=1<<i
            k+=1
    return tuple(rows)
def edge(g,i,j):return (g[i]>>j)&1

def maps(g, anti=False):
    n=len(g);deg=[x.bit_count() for x in g]
    poss=[[j for j in range(n) if deg[j]==(n-1-deg[i] if anti else deg[i])] for i in range(n)]
    order=sorted(range(n),key=lambda i:len(poss[i]));p=[-1]*n;ans=[]
    def go(k,used):
        if k==n:ans.append(tuple(p));return
        i=order[k]
        for j in poss[i]:
            if used>>j&1:continue
            if all(edge(g,i,h) ^ anti == edge(g,j,p[h]) for h in order[:k]):
                p[i]=j;go(k+1,used|1<<j);p[i]=-1
    go(0,0);return ans

def weighted_feedback(g,weights=None):
    n=len(g);weights=weights or [1]*n;N=1<<n;dp=[0]*N;last=[0]*N
    # Appending v last pays all of its outgoing edges to earlier vertices.
    subsets=[]
    for v in range(n):
        row=[0]*N
        bits=g[v]
        for s in range(1,N):
            k=s&-s;i=k.bit_length()-1
            row[s]=row[s^k]+(weights[i] if bits&k else 0)
        subsets.append(row)
    for s in range(1,N):
        best=10**12;vbest=-1;r=s
        while r:
            bit=r&-r;v=bit.bit_length()-1;r^=bit
            x=dp[s^bit]+weights[v]*subsets[v][s^bit]
            if x<best:best=x;vbest=v
        dp[s]=best;last[s]=vbest
    order=[];s=N-1
    while s:
        v=last[s];order.append(v);s^=1<<v
    return dp[-1],tuple(reversed(order))

def module(g,mask):
    return all((g[v]&mask) in (0,mask) for v in range(len(g)) if not mask>>v&1)

def components(g,mask):
    verts=[i for i in range(len(g)) if mask>>i&1]
    reach={i:(g[i]&mask)|1<<i for i in verts}
    for k in verts:
        for i in verts:
            if reach[i]>>k&1:reach[i]|=reach[k]
    cs=[];seen=0
    for i in verts:
        if not seen>>i&1:
            c=sum(1<<j for j in verts if reach[i]>>j&1 and reach[j]>>i&1)
            cs.append(c);seen|=c
    cs.sort(key=lambda c: -(g[(c&-c).bit_length()-1]&mask).bit_count())
    # Any representative in an earlier SCC has larger induced outdegree
    # than every representative in a later SCC, so this order is exact.
    return cs

def wword(g,mask):
    cs=components(g,mask)
    return tuple(c.bit_count() for c in cs) if all(c.bit_count() in (1,3) for c in cs) else None

def max_w(g):
    mods=[s for s in range(1,1<<len(g)) if module(g,s) and wword(g,s) is not None]
    return [s for s in mods if not any(s!=t and s&t==s for t in mods)]

def induced(g,verts):
    return tuple(sum(edge(g,i,j)<<b for b,j in enumerate(verts) if i!=j) for i in verts)

def join(gs):
    n=sum(map(len,gs));out=[0]*n;off=0
    for g in gs:
        for i,row in enumerate(g):out[off+i]=(row<<off)|sum(1<<j for j in range(off+len(g),n))
        off+=len(g)
    return tuple(out)

def word_graph(w):
    return join([graph(1,0) if x==1 else (2,4,1) for x in w])

R5=tuple(sum(1<<((i+j)%5) for j in (1,2)) for i in range(5))

def wrap(g):
    out=list(join([(0,),g,(0,)]));n=len(out)
    out[0]&=~(1<<(n-1));out[n-1]|=1
    return tuple(out)

def orbits(G,n):
    seen=set();ans=[]
    for i in range(n):
        if i not in seen:
            o={g[i] for g in G};seen|=o;ans.append(o)
    return ans

# All labelled tournaments through order five, then deterministic larger hostiles.
rng=random.Random(20261021)
samples=[graph(n,m) for n in range(1,6) for m in range(1<<(n*(n-1)//2))]
samples += [graph(n,rng.getrandbits(n*(n-1)//2)) for n in (6,7,8) for _ in range(35)]
# Diagonal synchronized triangles from inherited obstruction.
diag=[0]*6
for o in (0,3):
    for i in range(3):diag[o+i]|=1<<(o+(i+1)%3)
for i in range(3):
    for j in range(3):
        if i==j:diag[i]|=1<<(3+j)
        else:diag[3+j]|=1<<i
samples.append(tuple(diag))
slack_hist=Counter()
for g in samples:
    n=len(g);G=maps(g);O=orbits(G,n);oid={v:i for i,o in enumerate(O) for v in o}
    f,_=weighted_feedback(g);s=sum(len(o) for o in O if len(o)>1);d=3*f-s;slack_hist[d]+=1
    ck(d not in (1,2,5) and d>=0,'nonempty_spectrum')
    D=sum(3*(len(o)**2-1)//8-len(o) for o in O if len(o)>1);R=sum(len(o)>=5 for o in O)
    ck(0<=D<=d and (d-D)%3==0,'charge_congruence')
    blocks=max_w(g);reps=[(m&-m).bit_length()-1 for m in blocks];q=induced(g,reps)
    ws=[wword(g,m) for m in blocks];weights=[m.bit_count() for m in blocks]
    C,pi=weighted_feedback(q,weights);ck(f==C+sum(w.count(3) for w in ws),'module_feedback')
    order=[]
    for k in pi:
        for cm in components(g,blocks[k]):
            vv=[v for v in range(n) if cm>>v&1]
            _,ip=weighted_feedback(induced(g,vv));order.extend(vv[v] for v in ip)
    pos={v:i for i,v in enumerate(order)}
    ck(sum(edge(g,u,v) and pos[u]>pos[v] for u in range(n) for v in range(n))==f,'rebuilt_order_optimal')
    ext=[(u,v) for u in range(n) for v in range(n) if edge(g,u,v) and pos[u]>pos[v] and oid[u]!=oid[v]]
    e=len(ext);Bids={i for i,o in enumerate(O) if len(o)>=5}|{oid[v] for uv in ext for v in uv}
    B=set().union(*(O[i] for i in Bids)) if Bids else set();L=[k for k in pi if weights[k]>d]
    ck(e<=(d-D)//3 and len(Bids)<=R+2*e,'external_charge_bound')
    ck(all(not (blocks[k]&sum(1<<v for v in B)) for k in L),'long_disjoint_B')
    for oi in Bids:
        for k in L:
            pp=[pos[v] for v in O[oi]];ll=[pos[v] for v in range(n) if blocks[k]>>v&1]
            ck(max(pp)<min(ll) or min(pp)>max(ll),'B_orbit_no_straddle')
    separators=[]
    for i,j in zip(L,L[1:]):
        x=max(pos[v] for v in range(n) if blocks[i]>>v&1);y=min(pos[v] for v in range(n) if blocks[j]>>v&1)
        between=[v for v in B if x<pos[v]<y];ck(bool(between),'long_separator_exists');separators.append(oid[between[0]])
    ck(len(set(separators))==len(separators),'separator_orbits_distinct')
    m=(2*d//3+1 if d%3==0 else (2*d-2)//3 if d%3==1 else (2*d-7)//3)
    ck(len(L)<=R+2*(d-D)//3+1<=m,'sharp_long_slot_bound')

# Test wrapper on every labelled tournament of order2..5, including nonsymmetric U.
wrapper_count=0
for n in range(2,6):
    for mask in range(1<<(n*(n-1)//2)):
        g=graph(n,mask);h=wrap(g);G=maps(g);H=maps(h);A=maps(g,True);B=maps(h,True)
        fg,_=weighted_feedback(g);fh,_=weighted_feedback(h)
        ck(fh==fg+1,'wrapper_feedback')
        ck(len(H)==len(G) and all(p[0]==0 and p[-1]==n+1 for p in H),'wrapper_actual_group')
        ck({tuple(p[i+1]-1 for i in range(n)) for p in H}==set(G),'wrapper_group_restrictions')
        ck(bool(A)==bool(B) and len(A)==len(B),'wrapper_converse')
        ck(all(p[0]==n+1 and p[n+1]==0 for p in B),'wrapper_anti_endpoints')
        ck(len(components(h,(1<<len(h))-1))==1,'wrapper_strong')
        wrapper_count+=1

# Exact family controls: long words >d, full actual groups, quotient weighted DP,
# and exact finite maximal-W test on intended blocks (not a premise of injectivity).
def family(a,b,variant=0):
    d=3*a+4*b;L=2*a+b+1;w=(1,)*(d+1)
    words=[w for _ in range(L)]
    if variant and d>=3:words[0]=(3,)+(1,)*(d-2) # same vertex size, not a palindrome for d>=4
    blocks=[];types=[]
    def W(i):
        k=len(blocks);blocks.append(word_graph(words[i]));types.append(('w',i));return (0,),[k]
    def base():
        gs=[];ids=[]
        for i in range(b+1):
            x,k=W(a+i);gs.append(x);ids+=k
            if i<b:
                start=len(blocks)
                blocks.extend([(0,)]*5);types.extend([('r',i)]*5);gs.append(R5);ids+=list(range(start,start+5))
        return join(gs),ids
    def rec(k):
        if k==0:return base()
        left,li=W(a-k);center,ci=rec(k-1)
        x=len(blocks);blocks.append((0,));types.append(('x',k))
        y=len(blocks);blocks.append((0,));types.append(('y',k))
        right,ri=W(a+b+k)
        return join([left,wrap(center),right]),li+[x]+ci+[y]+ri
    q,ids=rec(a)
    # reorder actual blocks to the produced quotient order
    gs=[blocks[i] for i in ids];ts=[types[i] for i in ids];sizes=list(map(len,gs));N=sum(sizes);rows=[0]*N;off=0;starts=[]
    for g in gs:
        starts.append(off)
        for i,row in enumerate(g):rows[off+i]=row<<off
        off+=len(g)
    for i in range(len(q)):
        for j in range(len(q)):
            if edge(q,i,j):
                target=((1<<sizes[j])-1)<<starts[j]
                for v in range(starts[i],starts[i]+sizes[i]):rows[v]|=target
    return tuple(rows),q,gs,ts,words

family_rows=[]
for a,b in [(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(2,1),(0,2),(1,2)]:
    d=3*a+4*b;m=2*a+b+1
    if d==0:continue
    for variant in (0,1):
        g,q,gs,ts,words=family(a,b,variant)
        C,pi=weighted_feedback(q,list(map(len,gs)));t=sum(w.count(3) for w in words)
        ck(C==a+3*b,'family_weighted_feedback')
        G=maps(g);support=sum(any(p[i]!=i for p in G) for i in range(len(g)))
        ck(len(G)==5**b*3**t and support==5*b+3*t,'family_actual_group')
        ck(3*(C+t)-support==d,'family_exact_slack')
        anti=maps(g,True)
        expected=all(words[i]==words[-1-i][::-1] for i in range(m))
        ck(bool(anti)==expected,'family_actual_self_converse')
        # Any union of displayed blocks that lies in W would defeat maximality.
        for S in range(1,1<<len(q)):
            if S.bit_count()<2 or not module(q,S):continue
            cs=components(q,S)
            isin=all(c.bit_count()==1 or (c.bit_count()==3 and all(len(gs[i])==1 for i in range(len(q)) if c>>i&1)) for c in cs)
            ck(not isin,'family_canonical_W_maximality')
        # Rebuild B in this independently optimized block order and check that
        # its ACTUAL orbit separators attain the sharp bound, not just weights.
        n=len(g);O=orbits(G,n);oid={v:i for i,o in enumerate(O) for v in o}
        starts=[];off=0
        for block in gs:starts.append(off);off+=len(block)
        order=[]
        for i in pi:order.extend(range(starts[i],starts[i]+len(gs[i])))
        pos={v:i for i,v in enumerate(order)}
        ck(sum(edge(g,u,v) and pos[u]>pos[v] for u in range(n) for v in range(n))==C+t,'family_rebuilt_order_optimal')
        ext=[(u,v) for u in range(n) for v in range(n) if edge(g,u,v) and pos[u]>pos[v] and oid[u]!=oid[v]]
        e=len(ext);D=sum(3*(len(o)**2-1)//8-len(o) for o in O if len(o)>1);R=sum(len(o)>=5 for o in O)
        Bids={i for i,o in enumerate(O) if len(o)>=5}|{oid[v] for uv in ext for v in uv}
        B=set().union(*(O[i] for i in Bids)) if Bids else set()
        L=[i for i in pi if ts[i][0]=='w'];ck(len(L)==m,'family_long_slot_count')
        ck((e,D,R,len(Bids))==(a,4*b,b,2*a+b),'family_actual_orbit_charges')
        ck(all(not(set(range(starts[i],starts[i]+len(gs[i])))&B) for i in L),'family_long_disjoint_B')
        for oi in Bids:
            pp=[pos[v] for v in O[oi]]
            for i in L:
                ll=[pos[v] for v in range(starts[i],starts[i]+len(gs[i]))]
                ck(max(pp)<min(ll) or min(pp)>max(ll),'family_B_orbit_no_straddle')
        separators=[]
        for i,j in zip(L,L[1:]):
            x=max(pos[v] for v in range(starts[i],starts[i]+len(gs[i])))
            y=min(pos[v] for v in range(starts[j],starts[j]+len(gs[j])))
            between=[v for v in B if x<pos[v]<y]
            ck(bool(between),'family_long_separator_exists');separators.append(oid[between[0]])
        ck(len(set(separators))==len(separators)==m-1,'family_separator_orbits_distinct')
        ck(m==R+2*(d-D)//3+1,'family_sharp_bound_equality')
        family_rows.append({'a':a,'b':b,'variant':variant,'d':d,'vertices':len(g),'quotient_order':len(q),'long_slots':m,'feedback':C+t,'support':support,'aut_order':len(G),'self_converse':bool(anti)})

# Finite arithmetic orbit-charge optimization over every allowed multiset D<=300.
N=300;charges=[(3*(l*l-1)//8-l,l) for l in range(5,100,2) if 3*(l*l-1)//8-l<=N]
pen=[None]*(N+1);pen[0]=0
for D in range(N+1):
    if pen[D] is None:continue
    for q,l in charges:
        if D+q<=N:
            v=pen[D]+2*q-3
            if pen[D+q] is None or v<pen[D+q]:pen[D+q]=v
for d in range(N+1):
    avail=[(2*d-pen[D])//3+1 for D in range(d+1) if pen[D] is not None and (d-D)%3==0]
    if d in (1,2,5):ck(not avail,'penalty_empty_layer');continue
    b=d%3;a=(d-4*b)//3;m=2*a+b+1
    ck(a>=0 and max(avail)==m,'penalty_maximum')

# Small-prime irreducibility certificates, checked without symbolic factor command.
def pmul(f,g,p):
    h=[0]*(len(f)+len(g)-1)
    for i,x in enumerate(f):
        for j,y in enumerate(g):h[i+j]=(h[i+j]+x*y)%p
    return h
def divrem(f,g,p):
    f=f[:]
    while f and not f[-1]:f.pop()
    while len(f)>=len(g):
        a=f[-1]*pow(g[-1],-1,p)%p;shift=len(f)-len(g)
        for j,y in enumerate(g):f[j+shift]=(f[j+shift]-a*y)%p
        while f and not f[-1]:f.pop()
    return f
def irreducible(f,p):
    for d in range(1,(len(f)-1)//2+1):
        for coeff in product(range(p),repeat=d):
            if not divrem(f,list(coeff)+[1],p):return False
    return True
r=[-1,0,1,0,0,0,1];f3=[[1,0,1],[-1,0,-1,0,1]];f5=[[-1,-2,1,1],[1,-2,-1,1]]
for p,fs in [(3,f3),(5,f5)]:
    fs=[[x%p for x in f] for f in fs]
    ck(pmul(*fs,p)==[x%p for x in r],'irreducibility_factorization')
    for f in fs:ck(irreducible(f,p),'irreducibility_factor')
ck(not ({2,4}&{3}),'incompatible_factor_degrees')

result={'status':'PASS','checks':sum(checks.values()),'counts':dict(sorted(checks.items())),'upper_bound_actual_tournaments':len(samples),'upper_bound_slack_histogram':dict(sorted(slack_hist.items())),'wrapper_actual_tournaments':wrapper_count,'families':family_rows,'penalty_max_d':N,'scope':'Finite exact controls; unbounded proof audited separately.'}
print(json.dumps(result,indent=2,sort_keys=True))
