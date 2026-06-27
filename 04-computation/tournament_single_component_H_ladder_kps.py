"""
tournament_single_component_H_ladder_kps.py  (kind-pasteur-2026-06-27-S31ah)

Probe the SINGLE-COMPONENT H spectrum (connected Omega). By OCF, H = prod over
Omega-components of I(C_i,2); the "indecomposable" achievable values are the
single-component ones. Their GAPS are the irreducible impossibility certificates.

Reframing (S31ah): K_m realizable as Omega <=> 1+2m is single-component-achievable
<=> m not in {3,10} (the H=7, H=21 gaps). Here we map the single-component H spectrum
to find ALL such gaps (does anything beyond {7,21} appear? does the clique ladder cap?).
"""
import sys, random, itertools
from tournament_certificate_engine_kps import (
    all_tournaments, random_tournament, conflict_graph, independence_poly, H_value)

def omega_components(m, E):
    """connected components of Omega (graph on m vertices, adjacency E)."""
    seen=[False]*m; comps=[]
    for s in range(m):
        if seen[s]: continue
        stack=[s]; comp=[]; seen[s]=True
        while stack:
            v=stack.pop(); comp.append(v)
            for w in range(m):
                if E[v][w] and not seen[w]: seen[w]=True; stack.append(w)
        comps.append(comp)
    return comps

def is_clique(verts, E):
    return all(E[a][b] for a,b in itertools.combinations(verts,2))

def single_comp_H_and_clique(adj):
    """if Omega is connected, return (H, m, is_clique); else None."""
    m,E=conflict_graph(adj)
    if m==0: return None
    comps=omega_components(m,E)
    if len(comps)!=1: return None
    a=independence_poly(m,E)
    H=sum(a[k]*2**k for k in range(len(a)))
    return (H, m, is_clique(list(range(m)), E))

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(210)
    print("Single-component (connected Omega) H spectrum:")
    sc_H=set(); clique_sizes=set(); H_by_clique={}
    # exhaustive n<=6
    for n in range(3,7):
        for adj in all_tournaments(n):
            r=single_comp_H_and_clique(adj)
            if r:
                H,m,clq=r; sc_H.add(H)
                if clq: clique_sizes.add(m); H_by_clique[m]=H
    # sample n=7,8,9
    for n in (7,):
        for _ in range(30000):
            r=single_comp_H_and_clique(random_tournament(n,random))
            if r:
                H,m,clq=r; sc_H.add(H)
                if clq: clique_sizes.add(m); H_by_clique[m]=H
    mx=max(sc_H)
    gaps=[v for v in range(1,mx+1,2) if v not in sc_H]
    print(f"  single-component H achieved (sample, up to {mx}): {sorted(v for v in sc_H if v<=60)} ...")
    print(f"  ODD GAPS in single-component H below {mx}: {[g for g in gaps if g<=60]}")
    print(f"  => {7 in gaps and 21 in gaps and 'BOTH 7 and 21 are single-comp gaps (=K3,K10 non-realizable)'}")
    print(f"\n  CLIQUE-Omega (K_m) realizable sizes m: {sorted(clique_sizes)}")
    print(f"  K_m gives H=1+2m: {[(m,1+2*m) for m in sorted(clique_sizes)]}")
    missing_m=[m for m in range(1,13) if m not in clique_sizes]
    print(f"  K_m NOT found as Omega for m in {missing_m}  (expect 3,10 -> H=7,21; others may need bigger n/more samples)")
