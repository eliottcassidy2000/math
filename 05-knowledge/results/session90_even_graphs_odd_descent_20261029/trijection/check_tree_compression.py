#!/usr/bin/env python3
"""Natural first/last-edge compression for all spanning-tree orientations.
Standard library, explicit active checks; normal and -O runs must agree.
"""
import itertools, json, random, math
from fractions import Fraction
from collections import deque
CHECKS=0

def check(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)

def pairs(n): return list(itertools.combinations(range(n),2))

def tree_from_prufer(n,word):
    if n==1: return ()
    deg=[1]*n
    for v in word: deg[v]+=1
    edges=[]
    for v in word:
        u=next(i for i in range(n) if deg[i]==1)
        edges.append(tuple(sorted((u,v)))); deg[u]-=1; deg[v]-=1
    u,v=[i for i in range(n) if deg[i]==1]
    edges.append((u,v)); return tuple(sorted(edges))

def adj(n,edges):
    g=[set() for _ in range(n)]
    for u,v in edges: g[u].add(v); g[v].add(u)
    return g

def route(g,u,v):
    parent={u:None}; todo=deque([u])
    while todo:
        a=todo.popleft()
        if a==v: break
        for b in g[a]:
            if b not in parent: parent[b]=a; todo.append(b)
    out=[v]
    while out[-1]!=u: out.append(parent[out[-1]])
    return list(reversed(out))

def outer(g,e,f):
    candidates=[(len(route(g,u,v)),u,v) for u in e for v in f]
    best=max(x[0] for x in candidates)
    far=[(u,v) for d,u,v in candidates if d==best]
    if len(far)!=1: raise ValueError('not unique outer endpoints')
    return far[0]

def arc(mask,n,u,v):
    if u==v: raise ValueError('loop')
    a,b=sorted((u,v)); z=pairs(n).index((a,b)); x=(mask>>z)&1
    return x if u==a else 1-x

def compress(n,edges,mask):
    g=adj(n,edges); out=0
    for z,(i,j) in enumerate(pairs(n-1)):
        u,v=outer(g,edges[i],edges[j])
        out|=arc(mask,n,u,v)<<z
    return out

def expand(n,edges,tree_mask,edge_mask):
    m=n-1; g=adj(n,edges); target={e:i for i,e in enumerate(pairs(n))}; out=0
    for i,e in enumerate(edges): out|=((tree_mask>>i)&1)<<target[e]
    for z,(i,j) in enumerate(pairs(m)):
        u,v=outer(g,edges[i],edges[j]); x=(edge_mask>>z)&1
        if u>v: u,v=v,u; x^=1
        out|=x<<target[u,v]
    return out

def relabel(mask,n,p):
    out=0; target={e:i for i,e in enumerate(pairs(n))}
    for z,(u,v) in enumerate(pairs(n)):
        a,b=p[u],p[v]; x=(mask>>z)&1
        if a>b: a,b=b,a; x^=1
        out|=x<<target[a,b]
    return out

def run():
    rng=random.Random(9017); tree_counts={}; completion_counts={}
    for n in range(2,7):
        total=0; completions=0
        for word in itertools.product(range(n),repeat=n-2):
            es=tree_from_prufer(n,word); g=adj(n,es); total+=1
            observed={}
            for u,v in pairs(n):
                if (u,v) in es: continue
                r=route(g,u,v); e=tuple(sorted(r[:2])); f=tuple(sorted(r[-2:])); i,j=es.index(e),es.index(f)
                check(outer(g,e,f)==(u,v),'outer inverse')
                observed[tuple(sorted((i,j)))]=(u,v)
            check(set(observed)==set(pairs(n-1)),'nonedge/edge-pair bijection')
            if n<=5:
                # One arbitrary orientation per free completion, varied tree bits.
                for edge_mask in range(1<<((n-1)*(n-2)//2)):
                    tree_mask=(edge_mask*13+7)&((1<<(n-1))-1)
                    t=expand(n,es,tree_mask,edge_mask)
                    check(compress(n,es,t)==edge_mask,'round trip')
                    check(all(arc(t,n,u,v)==((tree_mask>>i)&1) for i,(u,v) in enumerate(es)),'fixed tree arcs')
                    completions+=1
            # All relabelings through n=5, five fixed-seed relabelings at n=6.
            perms=list(itertools.permutations(range(n))) if n<=5 else [tuple(rng.sample(range(n),n)) for _ in range(5)]
            mask=rng.randrange(1<<(n*(n-1)//2)); compressed=compress(n,es,mask)
            for p in perms:
                renamed=tuple(sorted(tuple(sorted((p[u],p[v]))) for u,v in es))
                pe=tuple(renamed.index(tuple(sorted((p[u],p[v])))) for u,v in es)
                check(compress(n,renamed,relabel(mask,n,p))==relabel(compressed,n-1,pe),'full naturality')
        tree_counts[n]=total; completion_counts[n]=completions
    # Exact marked stabilizer and Burnside counts for every oriented tree at n<=4.
    orbit_counts={}
    for n in range(2,5):
        accum=0
        for word in itertools.product(range(n),repeat=n-2):
            es=tree_from_prufer(n,word)
            for tm in range(1<<(n-1)):
                h=[]
                for p in itertools.permutations(range(n)):
                    if set(tuple(sorted((p[u],p[v]))) for u,v in es)!=set(es): continue
                    if not all(arc(expand(n,es,tm,0),n,p[u],p[v])==((tm>>i)&1) for i,(u,v) in enumerate(es)): continue
                    pe=tuple(es.index(tuple(sorted((p[u],p[v])))) for u,v in es)
                    h.append((p,pe))
                states=range(1<<((n-1)*(n-2)//2)); reps=set(); fixed=0
                for x in states:
                    reps.add(min(relabel(x,n-1,pe) for _,pe in h))
                    t=expand(n,es,tm,x)
                    a=sum(relabel(t,n,p)==t for p,_ in h)
                    b=sum(relabel(x,n-1,pe)==x for _,pe in h)
                    check(a==b,'marked stabilizer')
                    fixed+=b
                check(fixed==len(h)*len(reps),'marked Burnside')
                accum+=len(reps)
        orbit_counts[n]=accum
    # Forgetting the mark is not unit weight on unmarked tournament classes.
    marked_counts={}; masses={}
    for n in range(2,5):
        perms=list(itertools.permutations(range(n)))
        reps=sorted({min(relabel(t,n,p) for p in perms) for t in range(1<<(n*(n-1)//2))})
        trees=[tree_from_prufer(n,w) for w in itertools.product(range(n),repeat=n-2)]
        total=0; mass=Fraction(0)
        for t in reps:
            aut=[p for p in perms if relabel(t,n,p)==t]
            marked={}
            for es in trees:
                orbit=[tuple(sorted(tuple(sorted((p[u],p[v]))) for u,v in es)) for p in aut]
                key=min(orbit)
                if key not in marked: marked[key]=sum(z==es for z in orbit)
            total+=len(marked); mass+=sum((Fraction(1,a) for a in marked.values()),Fraction(0))
        check(mass==Fraction(n**(n-2)*2**(n*(n-1)//2),math.factorial(n)),'marked groupoid mass')
        marked_counts[n]=total; masses[n]=str(mass)
    check(marked_counts[3]==4,'marked versus unmarked hostile')
    return {'marked_tournament_counts':marked_counts,'marked_groupoid_masses':masses,'checks':CHECKS,'labelled_trees':tree_counts,'completion_roundtrips':completion_counts,'sum_over_labelled_oriented_tree_frames_of_marked_orbits':orbit_counts}

if __name__=='__main__': print(json.dumps(run(),sort_keys=True,indent=2))
