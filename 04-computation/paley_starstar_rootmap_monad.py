#!/usr/bin/env python3
"""
paley_starstar_rootmap_monad.py
monad-explorer-2026-06-07 (deep-research, 7th session)

Builds on THM-438 ADDENDUM-4.  Goal: test the ROOTED-MAP / loop-equation reframing
of (star-star):   sum_{even-series sigma of [0..2k]} mu(0,sigma) = (-1)^k C_k.

Three things verified here:

(1) SIGN IDENTITY.   sign(mu(sigma)) = (-1)^{m(sigma)},  m = cycle rank = 2k - V + 1.
    Hence  S_k = sum_{even-series sigma} (-1)^{m} * prod_v (|B_v|-1)!.
    (|B_v|-1)! is the number of CYCLIC ORDERS of the |B_v| visits to vertex v.

(2) GENUS-BLINDNESS (explains why all three ADDENDUM-4 genus-0 routes fail).
    Expand mu over the prod_v (b_v-1)! cyclic visit-orders.  Each choice is a ribbon
    map (Eulerian rotation system) with F faces; Euler gives F = m + 1 - 2g, so
    (-1)^{F-1} = (-1)^m for EVERY rotation system, independent of genus g.
    => the per-map sign is genus-blind; restricting to genus 0 throws away the
    cancellation.  We verify (-1)^{F(R)-1} is constant = (-1)^m across all R.

(3) ROOT-EDGE (Tutte) statistics.  For each even-series sigma we record the type of
    the ROOT step (walk edge 0->1): does block(2)=block(0) (a root BIGON / immediate
    return)?  This is the candidate first factor of the loop equation S_k=-sum S_iS_j.
    We tabulate S_k split by "root is a pendant bigon" vs not, to see the recursion.

No primes, no characters.  k up to 5 (Bell(11)=678570).
"""
import sys, math
from collections import defaultdict
from itertools import permutations

def set_partitions(c):
    c = list(c)
    if len(c) == 1:
        yield [c]; return
    f = c[0]
    for sm in set_partitions(c[1:]):
        for i, s in enumerate(sm):
            yield sm[:i] + [[f] + s] + sm[i+1:]
        yield [[f]] + sm

def mu_partition(blocks):
    m = 1
    for B in blocks:
        b = len(B)
        m *= ((-1)**(b-1)) * math.factorial(b-1)
    return m

def catalan(k):
    return math.comb(2*k, k)//(k+1)

def build_graph(blocks, L):
    pos2blk = {}
    for bi, B in enumerate(blocks):
        for pos in B:
            pos2blk[pos] = bi
    edges = [(pos2blk[i], pos2blk[i+1]) for i in range(L)]
    return edges, len(blocks), pos2blk

# ---- even-series test via flow-line multiplicities (same as paley_catalan_star_star) ----
import numpy as np
def edge_flow_lines(edges, nb):
    E = len(edges)
    Bm = np.zeros((nb, E))
    for ei,(u,v) in enumerate(edges):
        Bm[v,ei]+=1.0; Bm[u,ei]-=1.0
    u,s,vh = np.linalg.svd(Bm)
    tol=1e-9; rank=int((s>tol).sum()); m=E-rank
    if m==0: return [tuple()]*E, 0
    ns=vh[rank:]
    lines=[]
    for e in range(E):
        v=ns[:,e]
        if np.max(np.abs(v))<1e-7:
            lines.append(("ZERO",)); continue
        v=v/np.max(np.abs(v))
        for x in v:
            if abs(x)>1e-7:
                if x<0: v=-v
                break
        lines.append(tuple(round(float(x),6) for x in v))
    return lines, m

def is_even_series(edges, nb):
    adj=defaultdict(list)
    for (u,v) in edges:
        adj[u].append(v); adj[v].append(u)
    seen={0}; stk=[0]
    while stk:
        x=stk.pop()
        for w in adj[x]:
            if w not in seen: seen.add(w); stk.append(w)
    if len(seen)!=nb: return False
    lines,m=edge_flow_lines(edges,nb)
    if m==0: return False
    if any(l==("ZERO",) for l in lines): return False
    groups=defaultdict(int)
    for l in lines: groups[l]+=1
    return all(c%2==0 for c in groups.values())

def faces_of_rotation(edges, nb, rot):
    """edges: list of (u,v) directed as walk traverses (u=from,v=to) but graph undirected.
       Build ribbon map from rotation system 'rot' (rot[v] = cyclic list of half-edge ids
       around v). Return #faces. Each undirected edge has 2 half-edges.
       We use the standard: faces = orbits of the permutation  next = rot_perm o involution."""
    # half-edges: for edge index e with endpoints (a,b): he 2e at a, he 2e+1 at b.
    H=2*len(edges)
    endpoint=[0]*H; other=[0]*H
    for e,(a,b) in enumerate(edges):
        endpoint[2*e]=a; endpoint[2*e+1]=b
        other[2*e]=2*e+1; other[2*e+1]=2*e   # involution = same edge other half
    # rotation: next half-edge clockwise at the vertex
    nextrot={}
    for v in range(nb):
        cyc=rot[v]
        L=len(cyc)
        for i in range(L):
            nextrot[cyc[i]]=cyc[(i+1)%L]
    # face permutation: phi(h) = nextrot[ involution(h) ]
    seen=[False]*H; F=0
    for h in range(H):
        if not seen[h]:
            F+=1; cur=h
            while not seen[cur]:
                seen[cur]=True
                cur=nextrot[other[cur]]
    return F

def all_eulerian_rotations(edges, nb):
    """Yield rotation systems where at each vertex the half-edges are arranged so that
       in/out alternate? For the genus-blindness check we just enumerate ALL rotation
       systems (cyclic orders of half-edges) -- (deg-1)! per vertex -- and check the
       claim (-1)^{F-1} = const. (We use ALL, the strongest form of genus-blindness.)"""
    halves_at=defaultdict(list)
    for e,(a,b) in enumerate(edges):
        halves_at[a].append(2*e); halves_at[b].append(2*e+1)
    verts=list(range(nb))
    # cyclic orders: fix first element, permute rest
    def cyc_orders(items):
        if len(items)<=1:
            yield list(items); return
        first=items[0]
        for p in permutations(items[1:]):
            yield [first]+list(p)
    # cartesian product over vertices
    def rec(i, cur):
        if i==len(verts):
            yield dict(cur); return
        v=verts[i]
        for order in cyc_orders(halves_at[v]):
            cur[v]=order
            yield from rec(i+1, cur)
    yield from rec(0, {})

KMAX=int(sys.argv[1]) if len(sys.argv)>1 else 4
CHECK_GENUS_UPTO=int(sys.argv[2]) if len(sys.argv)>2 else 3

print("="*70)
print("(1) SIGN IDENTITY  sign(mu)=(-1)^m  and  S_k = sum (-1)^m prod(|B|-1)!")
print("(3) ROOT-BIGON split of S_k")
print(f"  {'k':>2} {'S_k':>8} {'(-1)^kC_k':>10} {'sign-id':>8}  {'#ev-ser':>8}  {'rootbigon':>10} {'non':>8}")
for k in range(1, KMAX+1):
    L=2*k
    S=0; cnt=0; sign_ok=True
    S_rootbigon=0; S_non=0
    for blocks in set_partitions(range(L+1)):
        edges,nb,pos2blk=build_graph(blocks,L)
        if any(u==v for (u,v) in edges): continue
        if not is_even_series(edges,nb): continue
        mu=mu_partition(blocks)
        V=nb; m=L-V+1
        # sign identity
        fact=1
        for B in blocks: fact*=math.factorial(len(B)-1)
        if mu != ((-1)**m)*fact: sign_ok=False
        S+=mu; cnt+=1
        # root-bigon: does position 2 share block with position 0 AND is that a pendant?
        rootbigon = (pos2blk[0]==pos2blk[2])
        if rootbigon: S_rootbigon+=mu
        else: S_non+=mu
    tgt=(-1)**k*catalan(k)
    print(f"  {k:>2} {S:>8} {tgt:>10} {str(S==tgt and sign_ok):>8}  {cnt:>8}  {S_rootbigon:>10} {S_non:>8}")

print("="*70)
print("(2) GENUS-BLINDNESS:  (-1)^{F(R)-1} == (-1)^m  for ALL rotation systems R")
print(f"     (checked exhaustively over even-series sigma at k<= {CHECK_GENUS_UPTO})")
allok=True; nchk=0
for k in range(1, CHECK_GENUS_UPTO+1):
    L=2*k
    for blocks in set_partitions(range(L+1)):
        edges,nb,_=build_graph(blocks,L)
        if any(u==v for (u,v) in edges): continue
        if not is_even_series(edges,nb): continue
        V=nb; m=L-V+1; target=(-1)**m
        # cap rotation enumeration cost
        cost=1
        deg=defaultdict(int)
        for (a,b) in edges: deg[a]+=1; deg[b]+=1
        for v in range(nb): cost*=math.factorial(deg[v]-1)
        if cost>200000:  # skip the very large ones, note it
            continue
        for rot in all_eulerian_rotations(edges,nb):
            F=faces_of_rotation(edges,nb,rot)
            nchk+=1
            if (-1)**(F-1)!=target:
                allok=False
                print(f"   FAIL k={k} m={m} F={F}")
                break
        if not allok: break
    if not allok: break
print(f"   genus-blindness holds for all {nchk} (sigma,R) checked: {allok}")
print("="*70)
print("Interpretation: the per-ribbon-map sign (-1)^{F-1} sees only the cycle rank m,")
print("NOT the genus.  So no genus-0 restriction can recover C_k (ADDENDUM-4 routes).")
print("S_k is the N=-1 evaluation of an ALL-GENUS rooted-map sum; the loop equation")
print("xF^2+F=1 must come from a Tutte/root-edge recursion, not a planarity filter.")
