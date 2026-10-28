#!/usr/bin/env python3
"""Exact, dependency-free checks for rigid sparse tournament carriers."""
import itertools, json, random, math

def need(ok, message):
    if not ok: raise RuntimeError(message)

def vertex(g):
    g.append(set()); return len(g)-1

def edge(g,u,v):
    g[u].add(v);g[v].add(u)

def arm(g,u,length):
    for _ in range(length):
        v=vertex(g);edge(g,u,v);u=v

def rigid_tree(word):
    g=[set()]
    for length in (1,2,3):arm(g,0,length)
    u=0
    for bit in word:
        v=vertex(g);edge(g,u,v);u=v
        if bit:arm(g,u,1)
    arm(g,u,4)
    return g

def append_graph(g,h):
    offset=len(g)
    g.extend({offset+v for v in row} for row in h)

def tournament(n,mask):
    a=[set() for _ in range(n)]
    for bit,(u,v) in enumerate(itertools.combinations(range(n),2)):
        if (mask>>bit)&1:a[u].add(v)
        else:a[v].add(u)
    return a

def encode(a):
    n=len(a);g=[];ports=[]
    for _ in a:
        c,t,h=[vertex(g) for _ in range(3)]
        edge(g,c,t);edge(g,c,h)
        for u,lengths in ((c,(1,)),(t,(1,2)),(h,(1,2,3))):
            for length in lengths:arm(g,u,length)
        ports.append((c,t,h))
    for u in range(n):
        for v in a[u]:edge(g,ports[u][1],ports[v][2])
    k=(n-1)*(n-2)//2
    length=max(1,(k-1).bit_length()) if k else 0
    for j in range(k):append_graph(g,rigid_tree([(j>>i)&1 for i in range(length)]))
    return g

def components(g):
    seen=set();out=[]
    for start in range(len(g)):
        if start in seen:continue
        todo=[start];seen.add(start);part=[]
        while todo:
            u=todo.pop();part.append(u)
            for v in g[u]-seen:seen.add(v);todo.append(v)
        out.append(part)
    return out

def decode(g):
    parts=components(g)
    if len(parts)==1:part=parts[0]
    else:
        cyc=[p for p in parts if sum(len(g[v]) for v in p)//2>=len(p)]
        need(len(cyc)==1,'unique cyclic core');part=cyc[0]
    core={u for u in part if len(g[u])>=3};roles={}
    for u in core:
        lengths=[]
        for v in g[u]-core:
            old=u;length=1
            while len(g[v])==2:
                w=next(iter(g[v]-{old}));old,v=v,w;length+=1
            need(len(g[v])==1,'pendant termination');lengths.append(length)
        signature=tuple(sorted(lengths))
        need(signature in ((1,),(1,2),(1,2,3)),'valid port role')
        roles[u]=len(signature)
    centers=sorted(u for u in core if roles[u]==1)
    tails=[];heads=[]
    for c in centers:
        tails.append(next(v for v in g[c]&core if roles[v]==2))
        heads.append(next(v for v in g[c]&core if roles[v]==3))
    return [{j for j,h in enumerate(heads) if h in g[t]} for t in tails]

def canonical_tournament(a):
    n=len(a)
    return min(tuple(int(p[j] in a[p[i]]) for i in range(n) for j in range(i+1,n)) for p in itertools.permutations(range(n)))

def aut_tournament(a):
    n=len(a)
    return sum(all((j in a[i])==(p[j] in a[p[i]]) for i in range(n) for j in range(n)) for p in itertools.permutations(range(n)))

def aut_graph(g):
    n=len(g);colors=[len(x) for x in g]
    while True:
        signatures=[(colors[u],tuple(sorted(colors[v] for v in g[u]))) for u in range(n)]
        palette={s:i for i,s in enumerate(sorted(set(signatures)))}
        new=[palette[s] for s in signatures]
        if len(set(new))==len(set(colors)):colors=new;break
        colors=new
    classes={c:[i for i in range(n) if colors[i]==c] for c in set(colors)}
    mapping={};used=set()
    def rec():
        if len(mapping)==n:return 1
        choices=[]
        for u in range(n):
            if u in mapping:continue
            cand=[v for v in classes[colors[u]] if v not in used and all((x in g[u])==(y in g[v]) for x,y in mapping.items())]
            if not cand:return 0
            choices.append((len(cand),-len(g[u]&mapping.keys()),u,cand))
        _,_,u,cand=min(choices);total=0
        for v in cand:
            mapping[u]=v;used.add(v);total+=rec();used.remove(v);del mapping[u]
        return total
    return rec()

def permute(g,p):
    out=[set() for _ in g]
    for u in range(len(g)):out[p[u]]={p[v] for v in g[u]}
    return out

def main():
    treechecks=0
    for length in range(1,7):
        # Intrinsic root and ordered backbone decode all words, plus graph automorphisms.
        for word in itertools.product((0,1),repeat=length):
            g=rigid_tree(word)
            need(sum(map(len,g))//2==len(g)-1,'tree size')
            need(aut_graph(g)==1,'rigid tree')
            treechecks+=1
    checks=0;autchecks=0;rng=random.Random(8101)
    for n in range(1,6):
        masks=range(1<<(n*(n-1)//2)) if n<=4 else [rng.randrange(1024) for _ in range(24)]
        for mask in masks:
            a=tournament(n,mask);g=encode(a)
            need(sum(map(len,g))//2==len(g)-1,'global sparse equation')
            need(canonical_tournament(decode(g))==canonical_tournament(a),'decode')
            p=list(range(len(g)));rng.shuffle(p)
            need(canonical_tournament(decode(permute(g,p)))==canonical_tournament(a),'relabeled decode')
            need(aut_graph(g)==aut_tournament(a),'exact automorphism order');autchecks+=1
            checks+=1
    # A repeated rigid component is not harmless: it creates an involution.
    g=rigid_tree((0,));append_graph(g,rigid_tree((0,)))
    need(aut_graph(g)==2,'duplicate padding hostile')
    print(json.dumps({'status':'PASS','rigid_trees':treechecks,'tournament_roundtrips':checks,'automorphism_order_checks':autchecks,'duplicate_tree_hostile':True},sort_keys=True))
if __name__=='__main__':main()
