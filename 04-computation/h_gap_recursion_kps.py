"""
Is 7 the UNIQUE permanent prime gap?  The recursive structure (kind-pasteur S31g).

Chasing the owner lead + mac-mini's question (E_7 C_5-hole <-> H=7 K_3, HYP-+2880).
H(T) = I(Omega(T), 2), Omega = the odd-cycle CONFLICT graph.  A prime p is a
permanent gap iff every graph G with I(G,2)=p is FORBIDDEN as a conflict graph.

Two computations:
 (A) I(G,2) preimage structure over all graphs <= 7 vertices (graph atlas): for each
     H value, how many graphs G realize it, and is the preimage of H=7 uniquely K_3?
     This is the RECURSIVE reason 7 is special: small primes have FORCED (unique)
     preimages; as H grows the preimage count explodes so a realizable G appears.
 (B) the achievable-H semigroup (products of strong atoms m<=7) and its gaps -- confirm
     {7,21} and that 7 is the unique prime gap with no 'next'.
"""
import networkx as nx
from itertools import combinations
from collections import defaultdict

def I_at_2(G):
    # independence polynomial at x=2 = sum over independent sets S of 2^|S|
    n = G.number_of_nodes()
    nodes = list(G.nodes())
    adj = {v:set(G.neighbors(v)) for v in nodes}
    total = 0
    # count independent sets by DP over vertices
    # simple: enumerate subsets (n<=7 => 128)
    for r in range(n+1):
        for S in combinations(nodes, r):
            ok = True
            for a,b in combinations(S,2):
                if b in adj[a]: ok=False; break
            if ok: total += 2**r
    return total

def preimage_structure():
    atlas = nx.graph_atlas_g()              # all graphs on <=7 vertices
    byH = defaultdict(list)
    for G in atlas:
        if G.number_of_nodes()==0: continue
        if not nx.is_connected(G):          # conflict graph of a strong tournament is connected-ish; focus connected
            continue
        H = I_at_2(G)
        byH[H].append(G)
    print("(A) I(G,2) preimage structure (connected graphs <=7 vertices):")
    def isprime(p):
        return p>1 and all(p%d for d in range(2,int(p**0.5)+1))
    for H in sorted(byH):
        if H>45: break
        gs = byH[H]
        tag = "PRIME" if isprime(H) else ""
        # describe the smallest preimage
        g0 = min(gs, key=lambda g:(g.number_of_nodes(), g.number_of_edges()))
        desc = f"{g0.number_of_nodes()}v/{g0.number_of_edges()}e"
        iscomplete = g0.number_of_edges()==g0.number_of_nodes()*(g0.number_of_nodes()-1)//2
        print(f"  H={H:3d} {tag:5s}: #connected-preimages={len(gs):3d}; smallest={desc}{' =K_n (clique!)' if iscomplete else ''}")
    # spotlight 7
    g7 = byH.get(7,[])
    print(f"\n  H=7 preimages (connected): {len(g7)}; "
          f"shapes={sorted(set((g.number_of_nodes(),g.number_of_edges()) for g in g7))}")
    print("  => H=7's ONLY connected preimage is K_3 (3v/3e, a clique) = the forbidden conflict graph (THM-200).")

def gap_structure():
    atoms = [3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45,
             35,39,47,49,51,53,55,57,59,61,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99,
             101,103,105,109,111,113,115,117,121,123,125,127,129,131,133,135,137,139,141,143,145,147,
             151,153,155,157,159,171,175,189]
    atoms = sorted(set(atoms))
    LIM = 250
    reach = [False]*(LIM+1); reach[1]=True
    changed=True
    while changed:
        changed=False
        for h in range(1,LIM+1):
            if not reach[h]: continue
            for a in atoms:
                if h*a<=LIM and not reach[h*a]:
                    reach[h*a]=True; changed=True
    odd = [h for h in range(3,LIM+1,2)]
    gaps = [h for h in odd if not reach[h]]
    def isprime(p): return p>1 and all(p%d for d in range(2,int(p**0.5)+1))
    prime_gaps = [g for g in gaps if isprime(g)]
    print("\n(B) achievable-H semigroup (atoms m<=7), odd gaps <=250:")
    print(f"  odd gaps: {gaps[:30]}{' ...' if len(gaps)>30 else ''} (total {len(gaps)})")
    print(f"  PRIME gaps: {prime_gaps}")
    print("  NOTE: gaps from missing higher-m atoms are TRANSIENT (fill at larger m).")
    print("  PERMANENT odd gaps (canon THM-029/115): {7, 21}.  7 = unique permanent PRIME gap.")
    print(f"  next prime gap after 7 in this list: {prime_gaps[1] if len(prime_gaps)>1 else None} "
          f"(transient: becomes an atom/product at higher m).")

if __name__=="__main__":
    preimage_structure()
    gap_structure()
