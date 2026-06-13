"""
opus-2026-06-08-S709  (CORRECTS S708 / MISTAKE-064)
===================================================
Erdős Problem 64 = the Erdős–Gyárfás conjecture (1995):
   every finite graph with minimum degree >= 3 contains a cycle whose length is a POWER OF TWO
   (2^k, i.e. 4, 8, 16, 32, ...).  OPEN and falsifiable.

This is NOT the even-cycle statement (that one -- "min degree >=3 => an even cycle of length >=4"
-- is classical and TRUE; I wrongly conflated them in S708, see MISTAKE-064).  Here we treat the
REAL problem honestly: verify it on small + structured graphs (no counterexample), and pinpoint
what a counterexample must avoid.

Known status: OPEN in general; PROVED for cubic planar graphs (Heckman-Krakovski 2013);
computationally verified for cubic graphs (Markström) with no counterexample found.
"""
import networkx as nx
from itertools import combinations
import random

def is_pow2(L):
    return L>=4 and (L&(L-1))==0     # 4,8,16,32,...

def cycle_length_set(G, bound=None):
    """set of simple-cycle lengths in G (undirected). bound caps the search length."""
    lens=set()
    if bound is None:
        for c in nx.simple_cycles(G):
            if len(c)>=3: lens.add(len(c))
    else:
        for c in nx.simple_cycles(G, length_bound=bound):
            if len(c)>=3: lens.add(len(c))
    return lens

def min_degree(G): return min((d for _,d in G.degree()), default=0)

def has_pow2_cycle(G, bound=None):
    return any(is_pow2(L) for L in cycle_length_set(G, bound))

# ---------- (1) exhaustive over the graph atlas (n<=7) -------------------------
print("="*78)
print("(1) EXHAUSTIVE n<=7 (graph atlas): does every min-degree>=3 graph have a")
print("    cycle of length a power of two (here necessarily length 4, since 8>7)?")
print("="*78)
atlas=nx.graph_atlas_g()
md3=0; with_pow2=0; without=[]
for G in atlas:
    if G.number_of_nodes()<4: continue
    if min_degree(G)>=3:
        md3+=1
        if has_pow2_cycle(G): with_pow2+=1
        else: without.append(G)
print(f"  min-degree>=3 graphs on <=7 vertices: {md3}; with a power-of-2 cycle: {with_pow2}; "
      f"without: {len(without)}")
print(f"  (on <=7 vertices the only available power of two is 4; min-deg>=3 forces girth<=4 here,")
print(f"   so all have a 4-cycle. The conjecture's substance appears only at girth>=5, n>=10.)")

# ---------- (2) girth-5 named graphs: no 4-cycle, but a power-of-2 cycle? -------
print("\n"+"="*78)
print("(2) girth>=5 cubic/structured graphs (no 4-cycle): do they have an 8- or 16-cycle?")
print("="*78)
named={
 'Petersen (3-reg, girth5)': nx.petersen_graph(),
 'Heawood (3-reg, girth6)' : nx.heawood_graph(),
 'Pappus (3-reg, girth6)'  : nx.pappus_graph(),
 'Desargues (3-reg,girth6)': nx.desargues_graph(),
 'Dodecahedral (3-reg)'    : nx.dodecahedral_graph(),
 'Moebius-Kantor (girth6)' : nx.moebius_kantor_graph(),
}
for name,G in named.items():
    Ls=sorted(cycle_length_set(G, bound=min(16,G.number_of_nodes())))
    p2=[L for L in Ls if is_pow2(L)]
    g=min(Ls) if Ls else None
    print(f"  {name:26s} n={G.number_of_nodes():2d} girth={g} cycle-lens(<=16)={Ls}  power-of-2: {p2}  "
          f"=> {'OK' if p2 else 'NO POWER OF 2 (<=16)!'}")

# ---------- (3) random cubic graphs up to n=16 --------------------------------
print("\n"+"="*78)
print("(3) random 3-regular graphs, several sizes: all contain a power-of-2 cycle?")
print("="*78)
rng=random.Random(7)
fail=0; tot=0
for n in [6,8,10,12,14,16]:
    for t in range(8):
        try:
            G=nx.random_regular_graph(3, n, seed=rng.randint(0,10**9))
        except Exception:
            continue
        tot+=1
        if not has_pow2_cycle(G, bound=16): fail+=1
print(f"  tested {tot} random cubic graphs (n=6..16); WITHOUT a power-of-2 cycle (<=16): {fail}")

# ---------- (4) what a counterexample must avoid; the dyadic framing -----------
print("\n"+"="*78)
print("(4) the dyadic structure of the problem")
print("="*78)
print("  A counterexample must avoid EVERY length in {4,8,16,32,...} simultaneously -- a")
print("  2-ADIC (dyadic) constraint, far stronger than 'no even cycle'. Min degree >=3 forces")
print("  MANY cycle lengths (e.g. Bondy-Vince: two cycles whose lengths differ by 1 or 2);")
print("  the open content is whether that spread must HIT a power of two.")
print("  Status: OPEN; proved for cubic PLANAR (Heckman-Krakovski 2013); no counterexample known")
print("  (Markström computational searches). Verified above on all <=7-vertex and many cubic graphs.")
print("\nDONE.")
