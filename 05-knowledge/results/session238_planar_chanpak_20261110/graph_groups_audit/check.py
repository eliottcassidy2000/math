#!/usr/bin/env python3
"""Independent construction and checks; imports no producer code. No assert gates."""
import itertools as it
import json
import sys
from collections import Counter


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pair(x, y):
    return frozenset((x, y))


def group(moduli):
    elements = list(it.product(*(range(k) for k in moduli)))
    def add(a, b):
        return tuple((x+y) % k for x,y,k in zip(a,b,moduli))
    return elements, add


def graph_orbits(vertices, elements, add):
    def act(g, v):
        return (v[0], add(g, v[1]))
    remaining = {pair(a,b) for a,b in it.combinations(vertices,2)}
    orbits = []
    while remaining:
        a,b = sorted(next(iter(remaining)))
        orbit = {pair(act(g,a),act(g,b)) for g in elements}
        remaining -= orbit
        orbits.append(orbit)
    return orbits


def build(moduli, copies, edges):
    elements, add = group(moduli)
    vertices = [(i,g) for i in range(copies) for g in elements]
    zero = elements[0]
    def ga(g,v):
        return (v[0],add(g,v[1]))
    def act(g,x):
        if x[0]=='s':
            return ('s',tuple(sorted(ga(g,v) for v in x[1])))
        return (x[0],ga(g,x[1]))
    edge_names = {f:('s',tuple(sorted(f))) for f in edges}
    cycles = {}
    private = {}
    for i in range(copies):
        v = (i,zero)
        incidence = sorted(edge_names[f] for f in edges if v in f)
        names = list('pqrab') + (['c'] if len(incidence)%2 else [])
        base = incidence + [(s,v) for s in names]
        for g in elements:
            w = ga(g,v)
            cycles[w] = [act(g,x) for x in base]
            private[w] = [(s,w) for s in 'pqr']
    universe = sorted(set(x for c in cycles.values() for x in c))
    arcs = set()
    for c in list(cycles.values())+list(private.values()):
        for a,b in zip(c,c[1:]+c[:1]):
            need((b,a) not in arcs,'contradictory prescription')
            arcs.add((a,b))
    prescribed = set(arcs)
    pending = {pair(a,b) for a,b in it.combinations(universe,2)} - {pair(a,b) for a,b in arcs}
    free_orbits = []
    while pending:
        a,b = sorted(min(pending,key=lambda f:tuple(sorted(f))))
        orbit = {(act(g,a),act(g,b)) for g in elements}
        need(len(orbit)==len(elements),'pair stabilizer')
        need(not any((b,a) in orbit for a,b in orbit),'pair reversal')
        need({pair(a,b) for a,b in orbit}<=pending,'partial pair orbit')
        pending -= {pair(a,b) for a,b in orbit}
        free_orbits.append(orbit)
        arcs |= orbit
    return vertices,universe,cycles,private,prescribed,free_orbits,arcs,elements,act


def audit(moduli,copies,edges):
    vs,us,cs,ds,ps,ors,arcs,els,act=build(moduli,copies,edges)
    if '--inject-failure' in sys.argv:
        a,b=min(ps)
        arcs.remove((a,b))
        arcs.add((b,a))
    need(ps <= arcs,'prescribed directed cycle arc lost')
    m,n,q=len(vs),len(us),len(els)
    o=sum(sum(v in f for f in edges)%2 for v in vs)
    need(n==len(edges)+5*m+o<=len(edges)+6*m,'vertex budget')
    need(len(ps)==2*len(edges)+6*m+o,'forced pairs')
    need(len(ors)*q==n*(n-1)//2-len(ps),'completion exponent')
    need(len(arcs)==n*(n-1)//2,'not a tournament')
    for g in els:
        need({(act(g,a),act(g,b)) for a,b in arcs}==arcs,'action not automorphic')
        if g!=els[0]:
            need(all(act(g,x)!=x for x in us),'nonfree action')
    selected=[set(cs[v]) for v in vs]+[set(ds[v]) for v in vs]
    for i,j in it.combinations(range(2*m),2):
        expected=(pair(vs[i],vs[j]) in edges) if j<m else (i<m and j==i+m)
        need(bool(selected[i]&selected[j])==expected,'wrong induced corona')
    for v in vs:
        a=('a',v)
        need(any((a,x) in arcs and (x,y) in arcs and (y,a) in arcs for x,y in it.permutations(cs[v],2)),'missing private tax triangle')
    return {'q':q,'m':m,'e':len(edges),'N':n,'k':len(ors)}


def hamilton(arcs,n):
    dp={(1<<v,v):1 for v in range(n)}
    for mask in range(1,1<<n):
        for v in range(n):
            count=dp.get((mask,v),0)
            if count:
                for w in range(n):
                    if not mask>>w&1 and (v,w) in arcs:
                        key=(mask|1<<w,w)
                        dp[key]=dp.get(key,0)+count
    return sum(dp.get(((1<<n)-1,v),0) for v in range(n))


def cycle_polynomial(arcs,n):
    # Cyclic rotations fixed by least vertex; same-support distinct cycles retained.
    cycles=[]
    for k in range(3,n+1,2):
        for subset in it.combinations(range(n),k):
            for tail in it.permutations(subset[1:]):
                c=(subset[0],)+tail
                if all((a,b) in arcs for a,b in zip(c,c[1:]+c[:1])):
                    cycles.append(sum(1<<v for v in c))
    states={(0,0):1}
    for support in cycles:
        extra=Counter()
        for (used,k),count in states.items():
            if not used&support:
                extra[(used|support,k+1)]+=count
        for key,count in extra.items():
            states[key]=states.get(key,0)+count
    coefficients=Counter()
    for (_,k),count in states.items():
        coefficients[k]+=count
    return coefficients


def main():
    rows=[]
    for q,copies in [(1,1),(1,2),(1,3),(3,1),(3,2),(5,1),(7,1)]:
        els,add=group((q,))
        vs=[(i,g) for i in range(copies) for g in els]
        gos=graph_orbits(vs,els,add)
        for word in range(1<<len(gos)):
            edges=set().union(*(o for i,o in enumerate(gos) if word>>i&1))
            rows.append(audit((q,),copies,edges))
    els,add=group((3,3))
    vs=[(0,g) for g in els]
    gos=graph_orbits(vs,els,add)
    for word in range(1<<len(gos)):
        edges=set().union(*(o for i,o in enumerate(gos) if word>>i&1))
        rows.append(audit((3,3),1,edges))
    vs,us,cs,ds,ps,ors,arcs,els,act=build((1,),1,set())
    hhist=Counter()
    for word in range(1<<len(ors)):
        arcs=set(ps)
        for i,orbit in enumerate(ors):
            arcs |= orbit if word>>i&1 else {(b,a) for a,b in orbit}
        index={v:i for i,v in enumerate(us)}
        indexed={(index[a],index[b]) for a,b in arcs}
        poly=cycle_polynomial(indexed,len(us))
        h=hamilton(indexed,len(us))
        need(h==sum(c*2**k for k,c in poly.items()),'OCF mismatch')
        need(poly[0]==1 and poly[1]>=3 and h>=7,'positive tax failure')
        hhist[h]+=1
    # Even-order control: regular C2 swaps its only pair, so neither orientation survives.
    need(all({(1-a,1-b) for a,b in orientation}!=orientation for orientation in [{(0,1)},{(1,0)}]),'even obstruction')
    print(json.dumps({'status':'PASS','graph_cases':len(rows),'cases':rows,'all_single_vertex_completions':sum(hhist.values()),'single_vertex_H_histogram':dict(hhist),'even_C2_control':'PASS'},sort_keys=True,indent=2))


if __name__=='__main__':
    main()
