#!/usr/bin/env python3
"""Exact, standard-library codecs into disconnected simple graphs with e=v-1.
Run python3 check_sparse_carrier.py and python3 -O check_sparse_carrier.py.
No removable assertions. All parsers validate their image, including padding.
"""
import itertools, json, random
from collections import Counter
CHECKS=0

def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok: raise RuntimeError(message)

def vertex(g):
    g.append(set()); return len(g)-1

def edge(g,a,b):
    if a==b or b in g[a]: raise ValueError('loop or duplicate edge')
    g[a].add(b); g[b].add(a)

def path(g,root,length):
    for _ in range(length):
        v=vertex(g); edge(g,root,v); root=v

def tag(g,root,lengths):
    for a in lengths: path(g,root,a)

def components(g):
    seen=set(); out=[]
    for v in range(len(g)):
        if v in seen: continue
        todo=[v]; seen.add(v); c=[]
        while todo:
            x=todo.pop(); c.append(x)
            for y in g[x]:
                if y not in seen: seen.add(y); todo.append(y)
        out.append(set(c))
    return out

def valid_graph(g):
    if not g: raise ValueError('empty carrier is not an image')
    for v,ns in enumerate(g):
        if v in ns or any(not isinstance(w,int) or w<0 or w>=len(g) or v not in g[w] for w in ns):
            raise ValueError('not a simple undirected graph')

def tree_code(g,c):
    # AHU center code, used only for the tree-padding image, not the core.
    if sum(len(g[v]) for v in c)!=2*(len(c)-1): raise ValueError('padding is not a tree')
    alive=set(c); deg={v:len(g[v]) for v in c}
    while len(alive)>2:
        leaves=[v for v in alive if deg[v]<=1]
        if not leaves: raise ValueError('no tree center')
        alive.difference_update(leaves)
        for v in leaves:
            for w in g[v]:
                if w in alive: deg[w]-=1
    def rooted(v,p): return '('+''.join(sorted(rooted(w,v) for w in g[v] if w!=p))+')'
    return min(rooted(v,-1) for v in alive)

def padding_tree(k,j):
    if not 0<=j<k: raise ValueError('bad padding index')
    ell=(k-1).bit_length(); L=ell+5
    g=[set() for _ in range(L+1)]
    for i in range(L): edge(g,i,i+1)
    edge(g,2,vertex(g))
    for i in range(ell):
        if (j>>i)&1: edge(g,3+i,vertex(g))
    return g

def append_graph(g,h):
    off=len(g); g.extend({off+w for w in ns} for ns in h)

def pad(g,k):
    for j in range(k): append_graph(g,padding_tree(k,j))
    return g

def pairs(n): return list(itertools.combinations(range(n),2))

def arc(mask,n,u,v):
    if u==v: return False
    a,b=sorted((u,v)); bit=(mask>>pairs(n).index((a,b)))&1
    return bool(bit) if u==a else not bool(bit)

def encode_tournament(n,mask):
    if not 0<=mask<(1<<(n*(n-1)//2)): raise ValueError('bad tournament')
    if n==0: return [set()]
    g=[]; centers=[]; tails=[]; heads=[]
    for _ in range(n):
        c,t,h=vertex(g),vertex(g),vertex(g)
        centers.append(c); tails.append(t); heads.append(h)
        edge(g,c,t); edge(g,c,h)
        tag(g,c,(1,2,3)); tag(g,t,(1,2,4)); tag(g,h,(1,2,5))
    for z,(u,v) in enumerate(pairs(n)):
        if not (mask>>z)&1: u,v=v,u
        edge(g,tails[u],heads[v])
    return pad(g,(n-1)*(n-2)//2)

def encode_graph(n,mask):
    if not 0<=mask<(1<<(n*(n-1)//2)): raise ValueError('bad graph')
    if n==0: return [set()]
    g=[]; hub=vertex(g); tag(g,hub,(1,2,3,5)); roots=[]
    for _ in range(n):
        v=vertex(g); roots.append(v); edge(g,hub,v); tag(g,v,(1,2,3,4))
    for z,(u,v) in enumerate(pairs(n)):
        if (mask>>z)&1: edge(g,roots[u],roots[v])
    return pad(g,mask.bit_count())

def parse(g,kind):
    valid_graph(g)
    if len(g)==1: return 0,0
    cs=components(g); mains=[c for c in cs if any(len(g[v])>=4 for v in c)]
    if len(mains)!=1: raise ValueError('expected one tagged component')
    main=mains[0]; core={v for v in main if len(g[v])>=4}; used=set(core); roles={}
    for v in core:
        lengths=[]
        for w in g[v]-core:
            p=v; x=w; length=0
            while True:
                if x in used or x in core: raise ValueError('overlap or cycle in pendant path')
                used.add(x); length+=1; nxt=g[x]-{p}
                if not nxt: break
                if len(nxt)!=1: raise ValueError('branched tag')
                p,x=x,next(iter(nxt))
            lengths.append(length)
        roles[v]=tuple(sorted(lengths))
    if used!=main: raise ValueError('unparsed main vertices')
    if kind=='tournament':
        allowed={(1,2,3),(1,2,4),(1,2,5)}
        if any(r not in allowed for r in roles.values()): raise ValueError('bad role')
        centers=sorted(v for v in core if roles[v]==(1,2,3)); n=len(centers)
        if n==0 or len(core)!=3*n: raise ValueError('bad core sizes')
        ts=[]; hs=[]
        for c in centers:
            ns=g[c]&core
            t=[v for v in ns if roles[v]==(1,2,4)]
            h=[v for v in ns if roles[v]==(1,2,5)]
            if len(ns)!=2 or len(t)!=1 or len(h)!=1: raise ValueError('bad matching')
            ts+=t; hs+=h
        if len(set(ts+hs))!=2*n: raise ValueError('shared port')
        expected={tuple(sorted((c,p))) for c,t,h in zip(centers,ts,hs) for p in (t,h)}
        mask=0
        for z,(u,v) in enumerate(pairs(n)):
            a=hs[v] in g[ts[u]]; b=hs[u] in g[ts[v]]
            if a==b: raise ValueError('missing or bidirected tournament arc')
            if a: mask|=1<<z; expected.add(tuple(sorted((ts[u],hs[v]))))
            else: expected.add(tuple(sorted((ts[v],hs[u]))))
        actual={(u,v) for u in core for v in g[u]&core if u<v}
        if actual!=expected: raise ValueError('unexpected core edge')
        k=(n-1)*(n-2)//2
    elif kind=='graph':
        hubs=[v for v in core if roles[v]==(1,2,3,5)]
        roots=sorted(v for v in core if roles[v]==(1,2,3,4))
        if len(hubs)!=1 or len(roots)+1!=len(core) or not roots: raise ValueError('bad cone roles')
        hub=hubs[0]; n=len(roots)
        if g[hub]&core!=set(roots): raise ValueError('bad cone hub')
        mask=sum((1<<z) for z,(u,v) in enumerate(pairs(n)) if roots[v] in g[roots[u]])
        k=mask.bit_count()
    else: raise ValueError('unknown kind')
    actual=Counter(tree_code(g,c) for c in cs if c!=main)
    expected=Counter(tree_code(h,set(range(len(h)))) for h in (padding_tree(k,j) for j in range(k)))
    if actual!=expected: raise ValueError('incorrect padding multiset')
    if sum(map(len,g))//2!=len(g)-1: raise ValueError('wrong Euler budget')
    return n,mask

def permute(g,p):
    h=[set() for _ in g]
    for u,ns in enumerate(g): h[p[u]]={p[v] for v in ns}
    return h

def transform(mask,n,p,tournament):
    out=0; pp=pairs(n); idx={e:i for i,e in enumerate(pp)}
    for z,(u,v) in enumerate(pp):
        b=(mask>>z)&1; a,c=p[u],p[v]
        if a>c:
            a,c=c,a
            if tournament: b^=1
        out|=b<<idx[(a,c)]
    return out

def canonical(mask,n,tournament):
    return min(transform(mask,n,p,tournament) for p in itertools.permutations(range(n)))

def generic_iso_count(g,h,limit=1000000):
    # Independent ordinary graph-isomorphism backtracker. No codec or parser.
    if len(g)!=len(h): return 0
    n=len(g); colors=[list(map(len,g)),list(map(len,h))]
    while True:
        sigs=[]
        for a,graph in enumerate((g,h)):
            sigs.append([(colors[a][v],tuple(sorted(Counter(colors[a][w] for w in graph[v]).items()))) for v in range(n)])
        palette={s:i for i,s in enumerate(sorted(set(sigs[0]+sigs[1])))}
        nxt=[[palette[s] for s in row] for row in sigs]
        if Counter(nxt[0])!=Counter(nxt[1]): return 0
        oldclasses=len(set(colors[0]+colors[1])); newclasses=len(palette)
        colors=nxt
        if oldclasses==newclasses: break
    classes={c:{v for v in range(n) if colors[1][v]==c} for c in set(colors[1])}
    mapping={}; unused=set(range(n)); total=0
    def rec():
        nonlocal total
        if len(mapping)==n:
            total+=1
            if total>limit: raise RuntimeError('isomorphism limit')
            return
        best=None; opts=None
        for v in range(n):
            if v in mapping: continue
            candidates=[w for w in classes[colors[0][v]]&unused if all((a in g[v])==(b in h[w]) for a,b in mapping.items())]
            if not candidates: return
            if opts is None or len(candidates)<len(opts): best,opts=v,candidates
            if len(opts)==1: break
        for w in opts:
            mapping[best]=w; unused.remove(w); rec(); unused.add(w); del mapping[best]
    rec(); return total

def reversal_character(g,p):
    # Royle's character: endpoint reversals, NOT permutation sign on edges.
    return -1 if sum(p[u]>p[v] for u in range(len(g)) for v in g[u] if u<v)%2 else 1

def input_graph(n,mask):
    g=[set() for _ in range(n)]
    for z,(u,v) in enumerate(pairs(n)):
        if (mask>>z)&1: edge(g,u,v)
    return g

def expected_lift(g,n,p,tournament):
    out=list(range(len(g)))
    if tournament:
        for u in range(n):
            for j in range(24): out[24*u+j]=24*p[u]+j
    elif n:
        for u in range(n):
            for j in range(11): out[12+11*u+j]=12+11*p[u]+j
    return out

def expect_reject(g,kind,label):
    try: parse(g,kind)
    except (ValueError,RecursionError): check(True,label); return
    check(False,'accepted hostile '+label)

def run():
    rng=random.Random(9001); stats={}; bank={}
    for kind,encoder,tournament in [('tournament',encode_tournament,True),('graph',encode_graph,False)]:
        count=0; maxvertices=0
        for n in range(6):
            reps={}
            for mask in range(1<<(n*(n-1)//2)):
                g=encoder(n,mask); check(parse(g,kind)==(n,mask),'round trip')
                check(sum(map(len,g))//2==len(g)-1,'edge budget')
                expected_components=1+((n-1)*(n-2)//2 if tournament and n else mask.bit_count() if not tournament else 0)
                check(len(components(g))==expected_components,'component budget')
                p=list(range(len(g))); rng.shuffle(p); decoded=parse(permute(g,p),kind)
                check(decoded[0]==n and canonical(decoded[1],n,tournament)==canonical(mask,n,tournament),'scrambled round trip')
                count+=1; maxvertices=max(maxvertices,len(g)); reps.setdefault(canonical(mask,n,tournament),mask)
            if n<=4: bank[kind,n]=list(reps)
        stats[kind]={'labelled_sources_through_5':count,'max_carrier_vertices':maxvertices}
    # Generic uncoloured graph isomorphism and automorphism checks, no parser.
    for (kind,n),reps in bank.items():
        enc=encode_tournament if kind=='tournament' else encode_graph
        tour=kind=='tournament'
        for i,a in enumerate(reps):
            g=enc(n,a)
            aut=sum(transform(a,n,p,tour)==a for p in itertools.permutations(range(n)))
            check(generic_iso_count(g,g)==aut,'full automorphism group size')
            for b in reps[i+1:]: check(generic_iso_count(g,enc(n,b))==0,'distinct source types')
    # Exact endpoint-reversal character, including even-order graph groups.
    character_tests=0
    for tournament,enc in [(True,encode_tournament),(False,encode_graph)]:
        for n in range(5):
            for mask in range(1<<(n*(n-1)//2)):
                g=enc(n,mask)
                for p in itertools.permutations(range(n)):
                    if transform(mask,n,p,tournament)!=mask: continue
                    lift=expected_lift(g,n,p,tournament)
                    check(permute(g,lift)==g,'explicit automorphism lift')
                    expected=1 if tournament else reversal_character(input_graph(n,mask),p)
                    check(reversal_character(g,lift)==expected,'Royle endpoint character')
                    character_tests+=1
    stats['character_preserving_lifts']=character_tests
    k2=input_graph(2,1)
    check(reversal_character(k2,(1,0))==-1,'K2 edge fixed but endpoint reversed')
    # All padding families through k=32, generic automorphism check samples.
    for k in range(1,33):
        hs=[padding_tree(k,j) for j in range(k)]
        codes=[tree_code(h,set(range(len(h)))) for h in hs]
        check(len(set(codes))==k,'distinct padding types')
        for j,h in enumerate(hs):
            check(max(map(len,h))==3,'padding degree')
            if k<=8 or j in (0,k-1): check(generic_iso_count(h,h)==1,'rigid padding')
    # Hostiles: preserve e=v-1 while corrupting padding, and corrupt core.
    for kind,enc in [('tournament',encode_tournament),('graph',encode_graph)]:
        n,mask=(4,0) if kind=='tournament' else (4,7)
        g=enc(n,mask); cs=components(g); main=max(cs,key=len)
        plain=[set(g[v]) for v in range(max(main)+1)]
        k=len(cs)-1
        for _ in range(k): append_graph(plain,padding_tree(k,0))
        check(sum(map(len,plain))//2==len(plain)-1,'hostile keeps Euler budget')
        expect_reject(plain,kind,'duplicated padding')
        bad=[set(ns) for ns in g]; bad[0].add(0); expect_reject(bad,kind,'loop')
        bad=[set(ns) for ns in g]; bad[0].discard(next(iter(bad[0]))); expect_reject(bad,kind,'asymmetric adjacency')
    # n=3 tournaments: Q(C3) sole nontrivial cycle has length9; Q(TT3) length5.
    for mask,want in [(0,5),(2,9)]:
        g=encode_tournament(3,mask); main=[c for c in components(g) if any(len(g[v])>=4 for v in c)][0]
        alive=set(main); deg={v:len(g[v]) for v in main}; todo=[v for v in main if deg[v]<=1]
        while todo:
            v=todo.pop()
            if v not in alive: continue
            alive.remove(v)
            for w in g[v]:
                if w in alive:
                    deg[w]-=1
                    if deg[w]==1: todo.append(w)
        check(len(alive)==want,'cycle length detects transitivity')
    stats['checks']=CHECKS
    return stats

if __name__=='__main__': print(json.dumps(run(),sort_keys=True,indent=2))
