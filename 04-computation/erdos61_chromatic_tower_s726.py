#!/usr/bin/env python3
"""
S726 — Erdos problem 61 (OPEN, infinite graph theory) via the chromatic tower from C_5.

61: if G1,G2 both have chromatic number aleph_1, must they have a COMMON subgraph G with chi(G)=4 (or
even aleph_0)? Erdos' remark: probably every aleph_1-chromatic graph contains EVERY chromatic-4 graph of
large girth (=> any two share one => 61 follows from this UNAVOIDABILITY).

Leveraging the repo (S725 + S720/721 + THM-209):
  - the unavoidable TARGETS (chromatic-4 triangle-free graphs) are a RECURSIVE TOWER from C_5 via the
    MYCIELSKI operation M: chi(M(G))=chi(G)+1, triangle-free preserved. C_5 = the Paley P_5 SEED of S725.
    M(C_5)=Grotzsch (chi=4, the canonical chromatic-4 triangle-free graph) = the finite seed of 61's
    target. The Mycielski tower is a graph TRANSFER (the S720/721 ladder); chi is the additive (+1) ladder.
  - "aleph_1-chromatic forces the whole tower" is the S725 renormalization: the uncountable (limit) object
    forces the finite/countable obstruction tower (cf. limit ordinals restore partition relations).
  - independence (THM-209, H=I(Omega,2)) tracks up the tower; common subgraph <-> shared tower level.

This session computes the Mycielski chromatic tower from C_5 (chi, girth, independence, the transfer
structure) -- the finite skeleton of 61's unavoidable targets -- and frames the reduction + heuristic.

No numpy/sympy. Exact.
"""
from itertools import combinations

def C5():
    adj={i:set() for i in range(5)}
    for i in range(5):
        adj[i].add((i+1)%5); adj[i].add((i-1)%5)
    return adj
def mycielski(adj):
    n=len(adj); new={}
    # originals 0..n-1, shadows n..2n-1, apex 2n
    for v in range(n): new[v]=set(adj[v])
    for v in range(n): new[n+v]=set()
    new[2*n]=set()
    for v in range(n):
        for u in adj[v]:
            new[v].add(n+u); new[n+u].add(v)   # shadow u' ~ v for each edge uv (and originals keep G edges)
        new[n+v].add(2*n); new[2*n].add(n+v)   # apex ~ all shadows
    return new
def chromatic_number(adj):
    n=len(adj); V=list(adj)
    def kcol(k):
        color={}
        order=sorted(V,key=lambda v:-len(adj[v]))
        def bt(i):
            if i==len(order): return True
            v=order[i]
            used={color[u] for u in adj[v] if u in color}
            for c in range(k):
                if c not in used:
                    color[v]=c
                    if bt(i+1): return True
                    del color[v]
            return False
        return bt(0)
    k=1
    while not kcol(k): k+=1
    return k
def girth(adj):
    n=len(adj); best=float('inf')
    from collections import deque
    for s in adj:
        dist={s:0}; par={s:None}; dq=deque([s])
        while dq:
            v=dq.popleft()
            for u in adj[v]:
                if u not in dist:
                    dist[u]=dist[v]+1; par[u]=v; dq.append(u)
                elif par[v]!=u:
                    best=min(best, dist[v]+dist[u]+1)
        # (BFS girth per source; min over sources)
    return best if best<float('inf') else None
def count_independent(adj):
    """return (alpha, total #independent sets, #independent sets of size 2 = I-coeff)."""
    V=list(adj)
    # branch-and-count
    best=[0]; total=[0]; bysize={}
    Vs=sorted(V,key=lambda v:-len(adj[v]))
    def rec(idx, chosen, forbidden):
        if idx==len(Vs):
            s=len(chosen); total[0]+=1; bysize[s]=bysize.get(s,0)+1; best[0]=max(best[0],s); return
        v=Vs[idx]
        # exclude v
        rec(idx+1, chosen, forbidden)
        # include v if allowed
        if v not in forbidden:
            rec(idx+1, chosen+[v], forbidden|adj[v]|{v})
    rec(0,[],set())
    return best[0], total[0], bysize.get(2,0), bysize

if __name__=="__main__":
    print("="*86)
    print("S726 — Erdos 61: the chromatic-4 targets are a Mycielski tower from C_5 (= Paley P_5, S725)")
    print("="*86)
    print("\n61 (OPEN): G1,G2 with chi=aleph_1 => common subgraph G with chi(G)=4 (or aleph_0)?")
    print("Erdos remark: probably every aleph_1-chromatic graph contains EVERY chromatic-4 LARGE-GIRTH")
    print("graph => any two share one => 61 reduces to UNAVOIDABILITY of chromatic-4-large-girth graphs.\n")

    g=C5(); names=["C_5 (=Paley P_5, S725 seed)","M(C_5)=Grotzsch","M^2(C_5)"]
    for lvl in range(3):
        n=len(g); chi=chromatic_number(g); gir=girth(g)
        alpha,tot,i2,bysize=count_independent(g)
        print(f"  level {lvl}: {names[lvl]:28} |V|={n:2d}  chi={chi}  girth={gir}  alpha={alpha}  "
              f"#ind-sets={tot}  I(G,2)-coeff_2={i2}")
        if lvl<2: g=mycielski(g)
    print("\n  => the Mycielski tower from C_5 climbs chi = 3,4,5,... (additive ladder), triangle-free")
    print("     preserved (girth >=4). LEVEL 1 = Grotzsch = the canonical chromatic-4 (girth-4) graph =")
    print("     the finite SEED of 61's unavoidable target. The tower is a graph TRANSFER (S720/721).")

    print("\nTHE REDUCTION + HEURISTIC (creative, leveraging S725):")
    print("  - 61's common-subgraph question <= UNAVOIDABILITY: every aleph_1-chromatic graph contains")
    print("    every chromatic-4-large-girth graph (then two such graphs share one).")
    print("  - the targets form a RECURSIVE TOWER from C_5 (Mycielski for girth 4; high-girth needs")
    print("    Ramanujan/probabilistic towers -- the hard 'large girth' case).")
    print("  - S725 RENORMALIZATION reading: aleph_1 (uncountable = a 'limit/fixed point') forces the")
    print("    whole finite obstruction tower, as limit ordinals restore partition relations (Chang).")
    print("    The finite levels are 'hot' (avoidable individually); the uncountable limit is the 'cold")
    print("    fixed point' that cannot avoid the tower -- the renormalization that should force the target.")
    print("  - COMMON-subgraph as INTERSECTION: G1 cap G2 contains the tower iff BOTH contain it; the")
    print("    unavoidability makes the intersection nonempty at chromatic-4. (THM-209 independence tracks")
    print("    the shared level.)")
    print("  - GIRTH is the obstruction (as in S722/S723: large girth = locally tree/ultrametric, the S725")
    print("    laminar structure). chi-4 + large girth = a high-temperature target; aleph_1 = the limit")
    print("    forcing it. The program: prove unavoidability of the Mycielski/high-girth tower in")
    print("    aleph_1-chromatic graphs via the renormalization/elementary-submodel-tower method.")
    print("="*86)
