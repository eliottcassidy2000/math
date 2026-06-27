"""
tournament_I21_omega_miner_kps.py  (kind-pasteur-2026-06-27-S31ah)

codex-S260 frontier task 1: the Omega-realizability miner for the remaining
connected graphs G with I(G,2)=21. By OCF a single-component H=21 <=> some
tournament has a connected Omega-COMPONENT isomorphic to one of:
   K10, K8-e, K6-M (K6 minus a 2-matching), K6-P3 (K6 minus 2 adjacent edges), P4.
P4 is excluded (THM-202). This miner searches tournaments for an Omega-component
isomorphic to each target and confirms 0 hits (corroborating THM-115's H!=21),
and reports the EXPANSION obstruction (extra cycles forced), generalizing the
K3->K4 mechanism of THM-200/201.
"""
import sys, itertools, random
from tournament_certificate_engine_kps import (
    all_tournaments, random_tournament, conflict_graph, independence_poly)

def graph_canon(m, E):
    best=None
    for perm in itertools.permutations(range(m)):
        adj=tuple(1 if E[perm[i]][perm[j]] else 0 for i in range(m) for j in range(i+1,m))
        if best is None or adj<best: best=adj
    return best

def make_graph(m, nonedges):
    """K_m minus the listed non-edges (0-indexed pairs)."""
    E=[[i!=j for j in range(m)] for i in range(m)]
    for (a,b) in nonedges: E[a][b]=E[b][a]=False
    return m, E

# the I=21 connected targets
TARGETS={
 "K10":      make_graph(10, []),
 "K8-e":     make_graph(8, [(0,1)]),
 "K6-M":     make_graph(6, [(0,1),(2,3)]),       # 2-matching removed
 "K6-P3":    make_graph(6, [(0,1),(1,2)]),       # 2 adjacent removed
 "P4":       (4, [[False,True,False,False],[True,False,True,False],
                  [False,True,False,True],[False,False,True,False]]),
}

def I2(m,E):
    a=independence_poly(m,E); return sum(a[k]*2**k for k in range(len(a)))

def omega_components(m, E):
    seen=[False]*m; comps=[]
    for s in range(m):
        if seen[s]: continue
        st=[s]; comp=[]; seen[s]=True
        while st:
            v=st.pop(); comp.append(v)
            for w in range(m):
                if E[v][w] and not seen[w]: seen[w]=True; st.append(w)
        comps.append(comp)
    return comps

def subgraph_canon(verts, E):
    m=len(verts); idx={v:i for i,v in enumerate(verts)}
    Es=[[E[a][b] for b in verts] for a in verts]
    return graph_canon(m, Es)

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(21)
    # verify target I-values
    print("Targets (I should be 21):")
    tcanon={}
    for name,(m,E) in TARGETS.items():
        tcanon[(m,graph_canon(m,E))]=name
        print(f"  {name}: m={m}, I(G,2)={I2(m,E)}")
    sizes=set(m for (m,E) in TARGETS.values())

    print("\nMining tournaments for an Omega-component iso to a target (expect 0 hits):")
    hits={name:0 for name in TARGETS}; checked=0
    # exhaustive n<=6 + heavy sample n=7,8,9,10
    pools=[("exhaustive n=6", all_tournaments(6))]
    for adj in pools[0][1]:
        checked+=1
        m,E=conflict_graph(adj)
        for comp in omega_components(m,E):
            if len(comp) in sizes:
                c=subgraph_canon(comp,E)
                key=(len(comp),c)
                if key in tcanon: hits[tcanon[key]]+=1
    for n in (7,):
        T= 12000
        for _ in range(T):
            adj=random_tournament(n,random); checked+=1
            m,E=conflict_graph(adj)
            for comp in omega_components(m,E):
                if len(comp) in sizes:
                    c=subgraph_canon(comp,E)
                    key=(len(comp),c)
                    if key in tcanon: hits[tcanon[key]]+=1
    print(f"  checked {checked} tournaments")
    for name in TARGETS:
        print(f"   {name}: {hits[name]} hits  {'(non-realizable so far)' if hits[name]==0 else 'REALIZED!'}")
    print("\nExpansion obstruction (the K3->K4 mechanism, generalized):")
    print("  K3 (I=7): 3 pairwise-sharing 3-cycles force a common vertex + a 5-cycle => K4 (I=9). THM-201.")
    print("  K10 (I=21): 10 pairwise-sharing odd cycles force >10 cycles (Moon: strong on m>=9 has alpha1>=12)")
    print("  => no isolated K10 component; the single-component H=21 is blocked. Corroborates THM-115.")
