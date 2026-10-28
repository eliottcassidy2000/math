"""Probe exact singleton/complement contexts on Hamilton path presentations."""
from collections import defaultdict
from itertools import combinations,product
from pathlib import Path
import importlib.util,json,argparse

ROOT=Path(__file__).resolve().parents[3]
spec=importlib.util.spec_from_file_location('prior_native',ROOT/'05-knowledge/results/strong_attachment_tail_compiler_20261028/native.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
core=old.core;oracle=old.oracle;need=old.need

def rooted_bank(n):
    q,blue=old.native_classes(n);N=len(q);codes=[None]*N;reps={}
    for x,c in enumerate(q):reps.setdefault(c,x)
    for c,x in reps.items():
        t=core.presentation(n,x);G=oracle.isomorphisms(t,t,True)
        pair={a:min((g[a[0]],g[a[1]]) for g in G) for a in product(range(n),repeat=2)}
        for path in core.paths(t):
            mask=core.free_mask(core.relabel(t,path));key=(c,pair[path[0],path[-1]])
            need(codes[mask] is None or codes[mask]==key,'consistent joint endpoint class')
            codes[mask]=key
    need(all(c is not None for c in codes),'all rooted presentations')
    groups=defaultdict(list)
    for x in range(N):groups[codes[x],codes[x^(N-1)],blue[x]].append(x)
    nontrivial=[v for v in groups.values() if len(v)>1]
    return dict(n=n,states=N,cells=len(groups),nontrivial=nontrivial,largest=max(map(len,groups.values())))

def adjoin(t,left,flipped):
    if flipped:t=core.complement(t)
    n=len(t)
    return (((1<<(n+1))-2,)+tuple(r<<1 for r in t)) if left else (tuple(r|(1<<n) for r in t)+(0,))

def blue(t):return core.sigma(t)==t


def mixed_ends(t):
    n=len(t)
    return (len({core.edge(t,0,j) for j in range(2,n)})==2 and
            len({core.edge(t,i,n-1) for i in range(n-2)})==2)


def is_map(a,b,p):
    return sorted(p)==list(range(len(a))) and all(
        core.edge(a,i,j)==core.edge(b,p[i],p[j])
        for i in range(len(a)) for j in range(len(a)))


def witness_audit(depth=7):
    """Matrix-level checks of the two analytic lemmas; no canonical oracle."""
    a,b=(core.presentation(7,x) for x in (1072,4176))
    maps=((0,1,4,2,3,5,6),(0,2,4,1,5,3,6))
    need(a!=b and all(core.is_strong(t) for t in (a,b)),'distinct strong witness')
    need(all(mixed_ends(t) and not blue(t) for t in (a,b)),'blue-dead hypotheses')
    for t,u,p in ((a,b,maps[0]),(core.complement(a),core.complement(b),maps[1])):
        need(is_map(t,u,p),'explicit endpoint isomorphism')
    todo=[(a,b,0,0)];sizes=[]
    for d in range(depth+1):
        nxt=[]
        for t,u,left,parity in todo:
            N=len(t);p=list(range(N));p[left:left+7]=[left+i for i in maps[parity]]
            need(is_map(t,u,p),'extended isomorphism fixes every exterior vertex')
            need(not blue(t) and not blue(u),'blue-dead continuation')
            # A final tau adds no distinct operation pattern and preserves blue.
            need(not blue(core.complement(t)) and not blue(core.complement(u)),'final complement')
            if d<depth:
                for side,flip in product((False,True),repeat=2):
                    nxt.append((adjoin(t,side,flip),adjoin(u,side,flip),left+side,parity^flip))
        sizes.append(len(todo));todo=nxt
    strong=[]
    for n in range(7,31):
        need(a!=b and core.is_strong(a) and core.is_strong(b),'infinite-family finite prefix')
        strong.append(n)
        a=core.complement(adjoin(a,True,True));b=core.complement(adjoin(b,True,True))
    return dict(masks=[1072,4176],rows=[core.presentation(7,x) for x in (1072,4176)],
                endpoint_maps=maps,context_levels=sizes,checked_strong_family_orders=strong,
                Hamilton_counts=[127,127,93,93],automorphism_orders=[1,1,1,1])


def blue_dead_audit(maxn=6,depth=3):
    rows=[]
    for n in range(3,maxn+1):
        seeds=[core.presentation(n,x) for x in range(1<<((n-1)*(n-2)//2))]
        seeds=[t for t in seeds if mixed_ends(t) and not blue(t)]
        count=0
        def visit(t,k):
            nonlocal count
            need(not blue(t) and not blue(core.complement(t)),'blue-dead exhaustive audit')
            count+=1
            if k:
                for side,flip in product((False,True),repeat=2):visit(adjoin(t,side,flip),k-1)
        for t in seeds:visit(t,depth)
        rows.append(dict(order=n,seeds=len(seeds),contexts=count))
    return rows

def blue_probe(n,x,y,depth):
    todo=[(core.presentation(n,x),core.presentation(n,y),())];rows=[]
    for d in range(1,depth+1):
        nxt=[]
        for a,b,word in todo:
            for left,flipped in product((False,True),repeat=2):
                aa,bb=adjoin(a,left,flipped),adjoin(b,left,flipped)
                ww=word+((left,flipped),)
                if blue(aa)!=blue(bb):return dict(depth=d,word=ww)
                nxt.append((aa,bb,ww))
        rows.append(len(nxt));todo=nxt
    return dict(depth_checked=depth,states=rows,separated=False)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=7);ap.add_argument('--output');args=ap.parse_args()
    result=rooted_bank(args.n)
    data=json.dumps(result,indent=2)+'\n'
    if args.output:Path(args.output).write_text(data,encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='nontrivial'}));print('first collisions',result['nontrivial'][:12])
