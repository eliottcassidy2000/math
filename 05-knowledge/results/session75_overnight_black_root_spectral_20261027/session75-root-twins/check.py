from itertools import permutations,product
import json

def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def strong(rows,vertices=None):
    n=len(rows);vs=list(range(n)) if vertices is None else list(vertices)
    if len(vs)<2:return True
    S=sum(1<<v for v in vs);start=1<<vs[0]
    for rev in (False,True):
        seen=start;todo=start
        while todo:
            b=todo&-todo;todo-=b;v=b.bit_length()-1
            out=(S^rows[v]^(1<<v)) if rev else (rows[v]&S)
            nxt=out&~seen;seen|=nxt;todo|=nxt
        if seen!=S:return False
    return True
def paths(rows):
    n=len(rows);mask=(1<<n)-1
    def go(p,left):
        if not left:yield tuple(p);return
        choices=left if not p else left&rows[p[-1]]
        while choices:
            b=choices&-choices;choices-=b;v=b.bit_length()-1
            yield from go(p+[v],left^b)
    yield from go([],mask)
def blue(rows,p):
    r=dict(zip(p,reversed(p)));n=len(rows)
    return all(((rows[i]>>j)&1)!=((rows[r[i]]>>r[j])&1) for i in range(n) for j in range(i+1,n))
def reverse_path(rows,p):
    out=list(rows)
    for a,b in zip(p,p[1:]):
        require((out[a]>>b)&1,('illegal',a,b))
        out[a]^=1<<b;out[b]^=1<<a
    return tuple(out)
def core(rows):
    n=len(rows)
    return [v for v,r in enumerate(rows) if r.bit_count() in(0,n-1) and strong(rows,[i for i in range(n) if i!=v])]
def cuts(rows):
    n=len(rows);full=(1<<n)-1
    return [sum((rows[i]&(full^A)).bit_count() for i in range(n) if A>>i&1) for A in range(1,full)]
def automorphisms(rows):
    n=len(rows);blocks={}
    for i,r in enumerate(rows):blocks.setdefault(r.bit_count(),[]).append(i)
    out=[]
    for pp in product(*(permutations(vs) for vs in blocks.values())):
        f={}
        for vs,targets in zip(blocks.values(),pp):f.update(zip(vs,targets))
        if all(((rows[i]>>j)&1)==((rows[f[i]]>>f[j])&1) for i in range(n) for j in range(i+1,n)):
            out.append(tuple(f[i] for i in range(n)))
    return out
old=(966,1004,1016,1009,995,193,384,768,544,96)
new=(1008,1001,995,997,494,448,640,768,576,48)
routes=[[(2,3,4,6,7,8,9,5,0,1),(3,0,2,5,9,8,7,6,4,1)],[(0,5,6,7,8,9,4,3,2,1),(2,3,7,9,8,6,5,0,4,1)]]
require(cuts(old)==cuts(new),'entire named cut equality')
results=[]
for rows,route in zip((old,new),routes):
    n=len(rows);require(strong(rows),'initial strong')
    P=list(paths(rows));auto=automorphisms(rows)
    black=[p for p in P if not blue(rows,p)]
    exits=[p for p in black if not strong(reverse_path(rows,p))]
    all_core=[p for p in P if core(reverse_path(rows,p))]
    require(len(P)==657 and not all_core,'H/core head')
    current=rows;states=[]
    for p in route:
        require(not blue(current,p),'route black')
        current=reverse_path(current,p);states.append(list(current))
    require(core(current),'route core')
    require([r.bit_count() for r in rows]==[6,7,7,7,7,3,2,2,2,2],'named scores')
    require(min(cuts(rows))==1,'lambda')
    results.append({'rows':list(rows),'H':len(P),'automorphisms':len(auto),'blue':len(P)-len(black),'black':len(black),'black_one_step_exits':len(exits),'one_step_core_any_colour':len(all_core),'two_black_core_route':[list(p) for p in route],'states':states})
require([(r['automorphisms'],r['blue'],r['black'],r['black_one_step_exits']) for r in results]==[(1,33,624,6),(9,99,558,0)],'exact expected counts')
print(json.dumps({'status':'PASS','all_named_cuts_compared':1022,'sources':results},sort_keys=True,indent=2))
