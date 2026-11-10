#!/usr/bin/env python3
import sys as _sys
if "--fail" in _sys.argv:
    raise RuntimeError("deliberate checker failure")
"""Independent checker: symbolic vertex labels, graph deletion recurrence, DFS cycles."""
import itertools as it,json,random
from functools import lru_cache
rng=random.Random(237)
def need(v,msg):
    if not v:raise RuntimeError(msg)
def make(m,edges,mode=0):
    vertices=[('edge',a,b) for a,b in edges]
    core=[];tri=[];extra=[];arcs=set()
    for i in range(m):
        incident=[v for v in vertices if v[0]=='edge' and i in v[1:]]
        ps=[('private',i,j) for j in range(5+len(incident)%2)]
        vertices+=ps;c=incident+ps;core.append(c);tri.append(ps[:3]);extra.append(ps[3])
        arcs.update(zip(c,c[1:]+c[:1]));arcs.add((ps[2],ps[0]))
    need(all((b,a) not in arcs for a,b in arcs),'arc conflict')
    nums={v:i for i,v in enumerate(vertices)};n=len(vertices)
    rows=[0]*n;free=[]
    for a,b in it.combinations(vertices,2):
        if (a,b) in arcs:src,dst=a,b
        elif (b,a) in arcs:src,dst=b,a
        else:
            free.append((nums[a],nums[b]))
            src,dst=(a,b) if mode==0 or rng.randrange(2)==0 else (b,a)
        rows[nums[src]]|=1<<nums[dst]
    return rows,[tuple(nums[x] for x in c) for c in core],[tuple(nums[x] for x in c) for c in tri],[nums[x] for x in extra],free

def poly_graph(adj):
    @lru_cache(None)
    def rec(mask):
        if not mask:return (1,)
        v=(mask&-mask).bit_length()-1
        no=rec(mask&~(1<<v));yes=rec(mask&~((1<<v)|adj[v]))
        ans=list(no)+[0]*max(0,len(yes)+1-len(no))
        for k,x in enumerate(yes):ans[k+1]+=x
        while len(ans)>1 and ans[-1]==0:ans.pop()
        return tuple(ans)
    return rec((1<<len(adj))-1)
def hcount(rows):
    n=len(rows)
    if not n:return 1
    @lru_cache(None)
    def tail(mask,v):
        if mask==(1<<n)-1:return 1
        options=rows[v]&~mask;ans=0
        while options:
            bit=options&-options;options-=bit
            ans+=tail(mask|bit,bit.bit_length()-1)
        return ans
    return sum(tail(1<<v,v) for v in range(n))
def cycles(rows):
    found=[]
    for start in range(len(rows)):
        def dfs(path,used):
            if len(path)>=3 and len(path)%2 and rows[path[-1]]>>start&1:found.append(tuple(path))
            for v in range(start+1,len(rows)):
                if not(used>>v&1) and rows[path[-1]]>>v&1:dfs(path+[v],used|1<<v)
        dfs([start],1<<start)
    return found

def verify(m,edges,mode=0):
    rows,C,D,A,free=make(m,edges,mode);n=len(rows)
    need(n==len(edges)+5*m+sum(sum(i in e for e in edges)%2 for i in range(m)),'size')
    supports=[set(c) for c in C+D]
    need(len({frozenset(s) for s in supports})==2*m,'distinct cycles')
    for c in C+D:need(len(c)%2==1 and all(rows[a]>>b&1 for a,b in zip(c,c[1:]+c[:1])),'cycle')
    graph=[0]*(2*m)
    for i,j in it.combinations(range(2*m),2):
        meet=bool(supports[i]&supports[j])
        expected=((i,j) in edges if j<m else i<m and j==i+m)
        need(meet==expected,'corona edges')
        if meet:graph[i]|=1<<j;graph[j]|=1<<i
    cg=poly_graph(graph)
    # Independent subset formula, evaluated by expansion with integer products.
    formula=[0]*(m+1)
    for bits in it.product((0,1),repeat=m):
        if any(bits[a] and bits[b] for a,b in edges):continue
        p=[1]
        for b in bits:
            if b:p=[0]+p
            else:p=[p[0]]+[p[j]+p[j-1] for j in range(1,len(p))]+[p[-1]]
        for i,a in enumerate(p):formula[i]+=a
    need(tuple(formula)==cg,'corona polynomial')
    gammas=[]
    for c,a,d in zip(C,A,D):
        poss=[(a,b,z) for b,z in it.permutations(c,2) if len({a,b,z})==3 and rows[a]>>b&1 and rows[b]>>z&1 and rows[z]>>a&1]
        need(bool(poss),'extra triangle')
        gammas.append(frozenset(poss[0]))
    need(len(set(gammas))==m and all(g not in supports for g in gammas),'extra triangles distinct and new')
    return n,rows

count=0;maxn=0
for m in range(6):
    pairs=list(it.combinations(range(m),2))
    for mask in range(1<<len(pairs)):
        edges=[p for i,p in enumerate(pairs) if mask>>i&1]
        n,rows=verify(m,edges);maxn=max(maxn,n);count+=1
random_completions=0
for m in range(1,8):
    for edges in ([],list(it.combinations(range(m),2))):
        for j in range(8):verify(m,edges,1);random_completions+=1
rows,C,D,A,free=make(1,[]);need(len(free)==4,'K1 free arcs')
values=[];cycle_numbers=[]
for bits in it.product((0,1),repeat=4):
    rr=rows[:]
    for (a,b),bit in zip(free,bits):
        rr[a]&=~(1<<b);rr[b]&=~(1<<a)
        if bit:rr[a]|=1<<b
        else:rr[b]|=1<<a
    cs=cycles(rr);adj=[0]*len(cs)
    for i,j in it.combinations(range(len(cs)),2):
        if set(cs[i])&set(cs[j]):adj[i]|=1<<j;adj[j]|=1<<i
    p=poly_graph(adj);hc=hcount(rr)
    perm=sum(all(rr[a]>>b&1 for a,b in zip(w,w[1:])) for w in it.permutations(range(5)))
    need(hc==perm==sum(a*2**k for k,a in enumerate(p)),'K1 OCF')
    need(p[1]>=3 and hc>=7,'K1 tax');values.append(hc);cycle_numbers.append(len(cs))
controls=[]
for m,edges in [(0,[]),(1,[]),(2,[]),(2,[(0,1)])]:
    rows,*_=make(m,edges);controls.append({'m':m,'edges':edges,'n':len(rows),'H':hcount(rows)})
need([r['H'] for r in controls]==[1,9,81,4617],'producer H values')
print(json.dumps({'status':'PASS','labelled_graphs_including_empty':count,'max_n_through_order5':maxn,'random_complete_and_empty_completions':random_completions,'K1_all16_H_values':values,'K1_odd_cycle_numbers':cycle_numbers,'hamiltonian_controls':controls},indent=2))
