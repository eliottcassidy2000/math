"""Exact controls for the proved feedback-changing four-black macro."""
from pathlib import Path
from itertools import product
from hashlib import sha256
import json

HERE=Path(__file__).resolve().parent
def need(b,s):
    if not b:raise RuntimeError(s)

def arc(t,u,v):return bool(t[u]>>v&1)

def move(t,p):
    need(len(p)==len(t) and len(set(p))==len(t),'not Hamiltonian')
    q=list(t)
    for u,v in zip(p,p[1:]):
        need(arc(t,u,v),'absent path arc')
        q[u]^=1<<v; q[v]^=1<<u
    return tuple(q)

def blue(t,p):
    return all(arc(t,p[i],p[j])==arc(t,p[-1-j],p[-1-i])
               for i in range(len(t)) for j in range(i+1,len(t)))

def join(a,b):
    m=len(a)
    return tuple(x|(((1<<len(b))-1)<<m) for x in a)+tuple(x<<m for x in b)

def presentations(m):
    # Standard path plus vertex 0 a source; every remaining bit is free.
    pairs=[(i,j) for i in range(1,m) for j in range(i+2,m)]
    for bits in product((0,1),repeat=len(pairs)):
        t=[sum(1<<j for j in range(i+1,m)) for i in range(m)]
        for (i,j),b in zip(pairs,bits):
            if b:t[i]^=1<<j;t[j]^=1<<i
        yield tuple(t)

def check_macro(a):
    m=len(a); p=tuple(range(m)); u=p[0]
    x,y,z,b0,b1,b2=range(m,m+6)
    tail=join((6,4,0),(2,4,1))
    state=join(a,tail)
    routes=[p+(x,y,z,b1,b2,b0),
            (b0,b1,z,y,x)+p[::-1]+(b2,),
            (x,b1,b0,b2)+p+(y,z),
            (b2,b0,b1,x,z,y)+p[::-1]]
    sums=[]
    for i,q in enumerate(routes):
        need(not blue(state,q),'blue move')
        pair=(q[1],q[-2]) if i==2 else (q[0],q[-1])
        value=sum(state[v].bit_count() for v in pair)
        sums.append(value)
        need(value!=len(state)-1,'failed score certificate')
        state=move(state,q)
    target=join(a,join((2,4,1),(4,1,2)))
    need(state==target,'literal target mismatch')
    need(sums==[m+6,3,5,m+6],'score formula differs')
    return routes,target

def words(n):
    if n==0:yield ()
    for a in (1,3,4):
        if a<=n:
            for w in words(n-a):yield (a,)+w

def active(w):return 4 in w or w!=w[::-1]
def key(w):return min(w,w[::-1])

def word_graph(n):
    vertices={key(w) for w in words(n) if active(w)}
    parent={w:w for w in vertices}
    def root(w):
        while parent[w]!=w:
            parent[w]=parent[parent[w]];w=parent[w]
        return w
    edge_count=0
    def edge(w,z):
        nonlocal edge_count
        if not active(z):return
        w,z=key(w),key(z)
        need(w in parent and z in parent,'macro endpoint absent')
        parent[root(w)]=root(z);edge_count+=1
    for w in vertices:
        for i,a in enumerate(w):
            if a==4:
                for v in ((1,3),(3,1)):edge(w,w[:i]+v+w[i+1:])
        for i in range(len(w)-1):
            if {w[i],w[i+1]}=={1,3}:
                edge(w,w[:i]+(w[i+1],w[i])+w[i+2:])
    for t in range(1,(n-1)//3):
        r=n-3*(t+1)
        need(r>=1,'source guard absent')
        a=(1,)*r+(3,)*(t-1)
        edge(a+(1,1,1,3),a+(3,3))
    components=len({root(w) for w in vertices})
    need(components==int(bool(vertices)),'golden language disconnected')
    return dict(n=n,active_merged_words=len(vertices),tested_macro_edges=edge_count,
                components=components)

def main():
    counts=[]
    for m in range(1,8):
        count=0
        for a in presentations(m):check_macro(a);count+=1
        need(count==2**max(0,(m-2)*(m-3)//2) if m>=2 else count==1,
             'incomplete universe')
        counts.append(dict(A_order=m,presentations=count))
    # Erasing the nonempty context is invalid, and its desired target is pure.
    empty_source=join((6,4,0),(2,4,1))
    empty_target=join((2,4,1),(4,1,2))
    need(blue(empty_source,(0,1,2,4,5,3)) is False,
         'unexpected empty-source standard path')
    # Source TT3->C3 is active; target C3->C3 is pure blue.
    paths0=[(0,1,2,4,5,3),(3,4,2,1,0,5),
            (0,4,3,5,1,2),(5,3,4,0,2,1)]
    state=empty_source; first_invalid=None
    for i,p in enumerate(paths0):
        missing=[(u,v) for u,v in zip(p,p[1:]) if not arc(state,u,v)]
        if missing:first_invalid=dict(step=i+1,missing_arcs=missing);break
        state=move(state,p)
    need(first_invalid is not None,'empty-context hostile lost')
    data=dict(scope='all fixed-standard-path A with source, orders1..7',
              contexts=counts,word_graph=[word_graph(n) for n in range(1,22)],
              empty_context_first_failure=first_invalid,
              source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    HERE.joinpath('inventory.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))

if __name__=='__main__':main()
