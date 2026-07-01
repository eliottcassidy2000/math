#!/usr/bin/env python3
"""
flip_rank_realizing_subcube_klein.py  --  klein-2026-07-01-S71

THE FLIP-RANK of tournament iso classes in the tiling/cube model.

Model: a labeled tournament on n vertices = an orientation of the C(n,2) edges = a bit per edge (i<j):
bit 1 = i->j, bit 0 = j->i. The NAIVE Hamiltonian-path tiling model FIXES the n-1 base-path arcs and
FLIPS the m=C(n-1,2) tiles. But the owner's n=4 observation: all 4 iso classes are reachable by flipping
only 2 arcs, if the 4 FIXED arcs are in a special configuration.

DEFINITION. A REALIZING SUBCUBE of dimension k = a choice of k "free" edges + a fixed orientation of the
other C(n,2)-k edges, whose 2^k completions include a representative of EVERY isomorphism class. The
FLIP-RANK rho(n) = the minimum such k. Lower bound: rho(n) >= ceil(log2 |G_n|) (A000568). This studies
rho(n), and the SHAPE of the fixed/free arcs that achieves it, as n grows.

Computes (n=4,5 exhaustive; n=6 randomized): rho(n); the achieving configs; the free-arc graph shape
(matching/path/star/triangle) and the fixed sub-tournament; redundancy 2^rho - |G_n|.
"""
import itertools, math, random

def edges(n): return [(i,j) for i in range(n) for j in range(i+1,n)]

def perm_maps(n, E):
    idx={e:k for k,e in enumerate(E)}
    maps=[]
    for pi in itertools.permutations(range(n)):
        m=[]  # for source edge k: (target k', flip)
        for (i,j) in E:
            a,b=pi[i],pi[j]
            if a<b: m.append((idx[(a,b)],0))
            else:   m.append((idx[(b,a)],1))
        maps.append(m)
    return maps

def canon(t, maps, ne):
    best=None
    for m in maps:
        v=0
        for k in range(ne):
            bit=(t>>k)&1
            tk,fl=m[k]
            if fl: bit^=1
            v|=bit<<tk
        if best is None or v<best: best=v
    return best

def classes_of(n):
    E=edges(n); ne=len(E); maps=perm_maps(n,E)
    cls={}; cid={}; nc=0
    for t in range(1<<ne):
        c=canon(t,maps,ne)
        if c not in cid: cid[c]=nc; nc+=1
        cls[t]=cid[c]
    return E, ne, cls, nc

def completions(free, fixed_bits):
    k=len(free)
    for a in range(1<<k):
        t=fixed_bits
        for b in range(k):
            if (a>>b)&1: t|=1<<free[b]
        yield t

def realizes(free, fixed_bits, cls, nc):
    seen=set()
    for t in completions(free, fixed_bits):
        seen.add(cls[t])
        if len(seen)==nc: return True
    return len(seen)==nc

def free_shape(free, E, n):
    fe=[E[k] for k in free]
    verts=set(v for e in fe for v in e)
    # component structure
    import collections
    adj=collections.defaultdict(list)
    for a,b in fe: adj[a].append(b); adj[b].append(a)
    seen=set(); comps=0; sizes=[]
    for v in verts:
        if v in seen: continue
        comps+=1; st=[v]; seen.add(v); sz=0
        while st:
            u=st.pop(); sz+=1
            for w in adj[u]:
                if w not in seen: seen.add(w); st.append(w)
        sizes.append(sz)
    degs=sorted(len(adj[v]) for v in verts)
    kind="matching" if max(degs)==1 else ("star" if (comps==1 and degs.count(1)==len(verts)-1) else ("triangle" if len(fe)==3 and len(verts)==3 and comps==1 else ("path" if max(degs)<=2 and comps==1 else "other")))
    return len(verts), comps, tuple(degs), kind

def flip_rank(n, exhaustive=True, tries=200000):
    E, ne, cls, nc = classes_of(n)
    lb=math.ceil(math.log2(nc))
    print(f"n={n}: |G_n|={nc} edges={ne} lower-bound ceil(log2)={lb}  naive C(n-1,2)={math.comb(n-1,2)}")
    for k in range(lb, ne+1):
        found=[]
        alledges=list(range(ne))
        if exhaustive and math.comb(ne,k)*(1<<(ne-k)) <= 3_000_000:
            for free in itertools.combinations(alledges,k):
                rest=[e for e in alledges if e not in free]
                for fa in range(1<<len(rest)):
                    fb=0
                    for b,e in enumerate(rest):
                        if (fa>>b)&1: fb|=1<<e
                    if realizes(list(free), fb, cls, nc):
                        found.append((free,fb))
            mode="exhaustive"
        else:
            seen_hit=set()
            for _ in range(tries):
                free=tuple(sorted(random.sample(alledges,k)))
                rest=[e for e in alledges if e not in free]
                fb=0
                for e in rest:
                    if random.random()<0.5: fb|=1<<e
                if realizes(list(free), fb, cls, nc):
                    found.append((free,fb));
                    if len(found)>=50: break
            mode=f"random({tries})"
        if found:
            print(f"  rho({n}) = {k}  ({'=LB, TIGHT' if k==lb else f'LB+{k-lb}'}); redundancy 2^{k}-|G|={(1<<k)-nc}; #achieving configs found ({mode}): {len(found)}")
            # characterize free-arc shapes among achievers
            shapes={}
            for free,fb in found[:5000]:
                s=free_shape(free,E,n); shapes[s[3]]=shapes.get(s[3],0)+1
            print(f"     free-arc shape kinds: {shapes}")
            # sample one config
            free,fb=found[0]; fe=[E[k] for k in free]
            print(f"     example: free arcs {fe} (shape {free_shape(free,E,n)}); a fixed orientation exists")
            return k, found, E, cls, nc
    return None

if __name__=="__main__":
    for n in [4,5]:
        flip_rank(n, exhaustive=True)
        print()
    # n=6: classes exhaustive, flip-rank randomized
    print("n=6 (randomized flip-rank search):")
    flip_rank(6, exhaustive=False, tries=60000)
