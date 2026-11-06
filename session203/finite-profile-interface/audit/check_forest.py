"""Independent endpoint-pattern audit against explicit successor forests."""
from collections import Counter
from functools import cache
from itertools import product
import json, random
B=-1

def need(q,m):
    if not q:raise RuntimeError(m)
def canonical(parts):return tuple(sorted(parts))
def intro(state,u,A):
    out=Counter()
    for parts,count in state.items():
        out[canonical(parts+((u,u),))]+=count
        for i,(s,e) in enumerate(parts):
            rest=parts[:i]+parts[i+1:]
            if e==B or e not in A:out[canonical(rest+((s,u),))]+=count
            if s in A:out[canonical(rest+((u,e),))]+=count
            for j,(t,f) in enumerate(parts):
                if i!=j and (e==B or e not in A) and t in A:
                    others=tuple(part for k,part in enumerate(parts) if k not in (i,j))
                    out[canonical(others+((s,f),))]+=count
    return out

def forget(state,u,p):
    out=Counter()
    for parts,count in state.items():
        parts=canonical(tuple((B if s==u else s,B if e==u else e) for s,e in parts))
        if sum(s==B for s,e in parts)<=p:out[parts]+=count
    return out

def cap(state,c):return Counter({x:min(n,c) for x,n in state.items()})
@cache
def forests(rows):
    n=len(rows);successor=[-1]*n;answer=[]
    def visit(v,incoming):
        if v==n:
            paths=[];seen=0
            for start in range(n):
                if incoming>>start&1:continue
                path=[];at=start
                while at!=-1:
                    if seen>>at&1:return
                    seen|=1<<at;path.append(at);at=successor[at]
                paths.append(tuple(path))
            if seen==(1<<n)-1:answer.append(tuple(paths))
            return
        successor[v]=-1;visit(v+1,incoming)
        for w in range(n):
            if rows[v]>>w&1 and not incoming>>w&1:
                successor[v]=w;visit(v+1,incoming|(1<<w))
        successor[v]=-1
    visit(0,0)
    return tuple(answer)
def exact(rows,active,p):
    inverse={vertex:label for label,vertex in active.items()}
    out=Counter()
    for paths in forests(tuple(rows)):
        parts=canonical(tuple((inverse.get(path[0],B),inverse.get(path[-1],B)) for path in paths))
        if sum(s==B for s,e in parts)<=p:out[parts]+=1
    return out

def hp(rows):
    n=len(rows);F=(1<<n)-1;total=0
    def walk(last,used):
        nonlocal total
        if used==F:total+=1;return
        for j in range(n):
            if rows[last]>>j&1 and not used>>j&1:walk(j,used|(1<<j))
    for r in range(n):walk(r,1<<r)
    return total

def read(state,r):
    h=s=L=positiveB=0
    for parts,n in state.items():
        if len(parts)==1:
            h+=n
            if parts[0][0]==r:s+=n
            else:positiveB+=n
        if len(parts) in (1,2) and any(e==r for a,e in parts):L+=n;positiveB+=n
    need(h+L-s==positiveB,'nonnegative B readout')
    return h,s,L,positiveB

def check_root(rows,state,active,p,r):
    for u in list(active):
        if u!=r:state=forget(state,u,p)
    v=active[r];ref=exact(rows,{r:v},p)
    need(state==ref,'root-only final pattern')
    if p>=2:
        h,s,L,Bval=read(state,r)
        edge=list(rows);n=len(edge);edge[v]|=1<<n;edge.append(((1<<n)-1)^(1<<v))
        need(h==hp(rows) and Bval==hp(edge),'literal H and edge-extension B')
        for K in (1,2,4):
            need(read(cap(state,K+1),r)[1]==min(s,K+1),'capped final start')
    return 1

labeled=patterns=capchecks=roots=0
for n in range(1,6):
    for bits in range(1<<(n*(n-1)//2)):
        rows=[0]*n;at=0
        for i in range(n):
            for j in range(i+1,n):
                if bits>>at&1:rows[i]|=1<<j
                else:rows[j]|=1<<i
                at+=1
        state=Counter({():1});capped=state.copy()
        for u in range(n):
            A={i for i in range(u) if rows[u]>>i&1}
            state=intro(state,u,A);capped=cap(intro(capped,u,A),3)
            need(capped==cap(state,3),'introduction cap commutes');capchecks+=1
            prefix=[r&((1<<(u+1))-1) for r in rows[:u+1]]
            need(state==exact(prefix,{i:i for i in range(u+1)},2),'prefix full endpoint patterns');patterns+=1
        for p in ((1,2,3) if n<=4 else (2,)):
            for mask in range(1<<n):
                active={i:i for i in range(n) if mask>>i&1};st=state;ct=cap(state,3)
                for u in range(n):
                    if not mask>>u&1:
                        st=forget(st,u,p);ct=cap(forget(ct,u,p),3)
                        need(ct==cap(st,3),'forget cap commutes');capchecks+=1
                need(st==exact(rows,active,p),'all active subsets and bulk multiplicities');patterns+=1
        for r in range(n):roots+=check_root(rows,state,{i:i for i in range(n)},2,r)
        labeled+=1
rng=random.Random(20315);reuse_words=operations=0
for W in (2,3,4):
    for trial in range(30):
        rows=[];active={};states={p:Counter({():1}) for p in (1,2,3)}
        caps={p:Counter({():1}) for p in (1,2,3)}
        for vertex in range(6):
            remove=[u for u in active if rng.randrange(3)==0]
            if len(active)-len(remove)==W:remove.append(rng.choice(list(active)))
            for u in sorted(set(remove)):
                del active[u]
                for p in states:
                    states[p]=forget(states[p],u,p);caps[p]=cap(forget(caps[p],u,p),3)
                    need(states[p]==exact(rows,active,p),'intermediate forget pattern');patterns+=1
                    need(caps[p]==cap(states[p],3),'reuse forget cap');capchecks+=1
                operations+=1
            free=list(set(range(W))-set(active));u=rng.choice(free)
            A={label for label in active if rng.randrange(2)}
            out={active[label] for label in A};newrow=0
            for old in range(vertex):
                if old in out:newrow|=1<<old
                else:rows[old]|=1<<vertex
            rows.append(newrow)
            for p in states:
                states[p]=intro(states[p],u,A);caps[p]=cap(intro(caps[p],u,A),3)
            active[u]=vertex
            for p in states:
                need(states[p]==exact(rows,active,p),'reused-label prefix pattern');patterns+=1
                need(caps[p]==cap(states[p],3),'reuse introduction cap');capchecks+=1
            operations+=1
        for r in active:roots+=check_root(rows,states[2],active,2,r)
        reuse_words+=1
print(json.dumps(dict(status='PASS',labeled_graphs=labeled,exact_prefix_pattern_checks=patterns,
    cap_commutation_checks=capchecks,root_readouts_with_literal_edge_H=roots,
    reused_label_words=reuse_words,reuse_operations=operations,maximum_reuse_order=6,
    literal_forest_graphs_cached=forests.cache_info().currsize),sort_keys=True,indent=2))
