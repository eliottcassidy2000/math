#!/usr/bin/env python3
"""Complete invariant-universe census for all C3/C3xC3 actions on six vertices.
Independent of repository atlas and nauty. Normal and -O are identical.
"""
import itertools, json
from collections import Counter
N=6
PAIRS=list(itertools.combinations(range(N),2)); IDX={e:i for i,e in enumerate(PAIRS)}
ID=tuple(range(N))
def require(c,s):
    if not c: raise RuntimeError(s)
def compose(p,q): return tuple(p[q[i]] for i in range(N))
def group(gens):
    out={ID}; todo=[ID]
    while todo:
        a=todo.pop()
        for b in gens:
            c=compose(a,b)
            if c not in out: out.add(c); todo.append(c)
    return out
def prep(p):
    ims=[]; rev=0; destrev=0
    for k,(a,b) in enumerate(PAIRS):
        x,y=p[a],p[b]; j=IDX[tuple(sorted((x,y)))]; ims.append(1<<j)
        if x>y: rev|=1<<k; destrev|=1<<j
    return p,ims,rev,destrev
def move(x,r,t=False):
    y=0
    while x:
        b=x&-x; y|=r[1][b.bit_length()-1]; x-=b
    return y^(r[3] if t else 0)
PERMS=[prep(p) for p in itertools.permutations(range(N))]
def sign(x,r): return -1 if (x&r[2]).bit_count()%2 else 1
def canonical(x,t): return min(move(x,r,t) for r in PERMS)
def even(x,rs): return all(sign(x,r)==1 for r in rs if move(x,r)==x)
def invariant_masks(H,t):
    # Pair orbits; orient each from one directed representative consistently.
    seen=set(); orbitmasks=[]; tbase=0
    for a,b in PAIRS:
        if (a,b) in seen: continue
        directed={(h[a],h[b]) for h in H}
        require(not any((y,x) in directed for x,y in directed),'odd action no edge reversal')
        orb={tuple(sorted(e)) for e in directed}; seen|=orb
        mask=sum(1<<IDX[e] for e in orb); orbitmasks.append(mask)
        tbase|=sum(1<<IDX[tuple(sorted((x,y)))] for x,y in directed if x>y)
    vals=[]
    for sel in range(1<<len(orbitmasks)):
        m=tbase if t else 0
        for i,orb in enumerate(orbitmasks):
            if sel>>i&1: m^=orb
        vals.append(m)
    return vals

def case(name,gens):
    H=group(gens)
    norm=[]
    for r in PERMS:
        p=r[0]; inv=tuple(p.index(i) for i in range(N))
        if {compose(compose(p,h),inv) for h in H}==H: norm.append(r)
    ts=invariant_masks(H,True); gs=invariant_masks(H,False)
    tr={canonical(t,True) for t in ts}
    er={canonical(x,False) for x in gs if even(x,PERMS)}
    nt={min(move(t,r,True) for r in norm) for t in ts}
    ne={min(move(x,r) for r in norm) for x in gs if even(x,norm)}
    escaped=[x for x in sorted(ne) if not even(x,PERMS)]
    tm=Counter(canonical(t,True) for t in nt)
    em=Counter(canonical(x,False) for x in ne if even(x,PERMS))
    records=[]
    for t in sorted(tr):
        aut=[r[0] for r in PERMS if move(t,r,True)==t]
        seen=set(); orb=[]
        for v in range(N):
            if v in seen: continue
            s={p[v] for p in aut}; seen|=s; orb.append(len(s))
        records.append({'mask':t,'aut_order':len(aut),'aut_vertex_orbits':sorted(orb),'normalizer_multiplicity':tm[t]})
    require(len(nt)==len(ne),'relative Burnside')
    require(len(tr)==len(er),'p-group theorem exact result')
    return {'action':name,'order':len(H),'normalizer_order':len(norm),'invariant_tournaments':len(ts),'invariant_graphs':len(gs),
    'tournament_types':len(tr),'even_graph_types':len(er),'normalizer_tournament_orbits':len(nt),'normalizer_even_graph_orbits':len(ne),
    'globally_odd_relative_even_masks':escaped,'tournaments':records,'even_graph_masks':sorted(er),'even_graph_normalizer_multiplicity':dict(sorted(em.items()))}

def main():
    a=(1,2,0,3,4,5); b=(0,1,2,4,5,3)
    out=[case('single_3_cycle',[a]),case('double_3_cycle',[compose(a,b)]),case('independent_3_cycles',[a,b])]
    require([d['tournament_types'] for d in out]==[11,3,1],'expected predicate counts')
    T=[{x['mask'] for x in d['tournaments']} for d in out]; E=[set(d['even_graph_masks']) for d in out]
    require(T[0]&T[1]==T[2] and E[0]&E[1]==E[2],'joint predicate intersection')
    print(json.dumps({'status':'PASS','universe':'all invariant labelled objects; full S6 equivalence and automorphism signs','cases':out},indent=2))
if __name__=='__main__': main()
