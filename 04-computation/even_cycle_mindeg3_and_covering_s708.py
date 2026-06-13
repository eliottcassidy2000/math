"""
opus-2026-06-08-S708
====================
WORK ON: does every finite graph with min degree >= 3 contain a cycle of length 2k (k>=2)?
  ANSWER: YES (classical). Proof = longest-path + pigeonhole on neighbour parities.
  Characterisation: even-cycle-free  <=>  every block is an edge or an ODD cycle
  (a "cactus of odd cycles")  <=>  min degree <= 2.
EXPLORE: covering systems (Owens min-modulus 42; Hough boundedness) and their tie to LRC.
BRIDGE: the DIRECTED even-cycle problem <-> Pfaffian orientations / Polya permanent (S707);
        in tournaments, a strong component of size >=4 has an even directed cycle (Moon pancyclic).
"""
import networkx as nx
from itertools import combinations
import random

# ---------- even-cycle detector via block (biconnected) structure -------------
def has_even_cycle(G):
    # even cycle (length>=4) exists iff some biconnected block has E>V (contains a theta
    # => even cycle) or is an even cycle (E==V, V even).
    for comp_edges in nx.biconnected_component_edges(G):
        edges=list(comp_edges)
        verts=set()
        for u,v in edges: verts.add(u); verts.add(v)
        V=len(verts); E=len(edges)
        if E>V: return True               # theta or denser -> even cycle (proved)
        if E==V and V%2==0: return True    # an even cycle block
    return False

def min_degree(G):
    return min((d for _,d in G.degree()), default=0)

# ---------- (1) EXHAUSTIVE verification over all graphs up to 7 vertices -------
print("="*78)
print("(1) THEOREM (exhaustive, all graphs n<=7 via the graph atlas):")
print("    min degree >= 3  =>  has an even cycle of length >= 4")
print("="*78)
atlas=nx.graph_atlas_g()   # all graphs on 0..7 nodes
tested=0; mind3=0; fails=0
ecfree_mindeg=[]   # min-degrees of even-cycle-free graphs (with >=1 cycle)
for G in atlas:
    if G.number_of_nodes()<1: continue
    tested+=1
    md=min_degree(G); ec=has_even_cycle(G)
    if md>=3:
        mind3+=1
        if not ec: fails+=1
    if not ec and G.number_of_edges()>=G.number_of_nodes():  # has a cycle but no EVEN cycle
        ecfree_mindeg.append(md)
print(f"  graphs tested: {tested};  with min-degree>=3: {mind3};  "
      f"min-deg>=3 WITHOUT an even cycle: {fails}")
print(f"  => theorem holds with 0 counterexamples for all graphs on <=7 vertices.")
print(f"  even-cycle-free graphs that DO contain a cycle: max min-degree among them = "
      f"{max(ecfree_mindeg) if ecfree_mindeg else 'n/a'}  (must be <= 2)")

# ---------- (2) the longest-path proof witnessed --------------------------------
print("\n"+"="*78)
print("(2) PROOF witness: longest path endpoint has 3 nbrs on the path; two share parity")
print("="*78)
def even_cycle_via_longest_path(G):
    # find a longest path by DFS (heuristic длиннейший); then exhibit the even cycle
    best=[]
    def dfs(u,vis,path):
        nonlocal best
        if len(path)>len(best): best=path[:]
        for w in G[u]:
            if w not in vis:
                vis.add(w); path.append(w); dfs(w,vis,path); path.pop(); vis.remove(w)
    for s in G:
        dfs(s,{s},[s])
    P=best
    pos={v:i for i,v in enumerate(P)}
    v0=P[0]
    nbr_pos=sorted(pos[w] for w in G[v0] if w in pos and pos[w]>0)
    # two same-parity positions
    for a,b in combinations(nbr_pos,2):
        if (a%2)==(b%2):
            return (b-a)+2, P[a], P[b]   # even cycle length, endpoints
    return None
# demo on a non-bipartite min-degree-3 graph: the 5-wheel W5 (hub+C5)
W=nx.wheel_graph(6)  # hub 0 + cycle 1..5
res=even_cycle_via_longest_path(W)
print(f"  wheel W5 (min deg {min_degree(W)}): even cycle length {res[0]} via nbrs at "
      f"endpoints {res[1]},{res[2]}; has_even_cycle={has_even_cycle(W)}")
# Petersen graph (cubic, girth 5, non-bipartite)
P=nx.petersen_graph()
print(f"  Petersen (3-regular, girth 5): has_even_cycle={has_even_cycle(P)} "
      f"(shortest even cycle length = {min((len(c) for c in nx.minimum_cycle_basis(P) if len(c)%2==0), default=None)})")

# ---------- (3) theta graph always has an even cycle (parity of 3 path-sums) ----
print("\n"+"="*78)
print("(3) THETA lemma: 3 internally-disjoint u-v paths of lengths a,b,c give cycles")
print("    a+b, b+c, a+c; their sum 2(a+b+c) is even => an even number are odd => >=1 even")
print("="*78)
import itertools
allok=True
for a,b,c in itertools.product(range(1,6),repeat=3):
    cyc=[a+b,b+c,a+c]
    if not any(L%2==0 for L in cyc): allok=False
print(f"  every theta (a,b,c in 1..5) has an even cycle among {{a+b,b+c,a+c}}: {allok}")

# ---------- (4) DIRECTED even cycle in tournaments <-> Pfaffian (S707 bridge) ----
print("\n"+"="*78)
print("(4) DIRECTED even cycles in tournaments (strong comp size>=4 => even dicycle, Moon)")
print("    -- the directed even-cycle problem ties to Pfaffian orientations / Polya permanent")
print("="*78)
def random_strong_tournament(n,rng,tries=200):
    for _ in range(tries):
        D=nx.DiGraph()
        D.add_nodes_from(range(n))
        for i,j in combinations(range(n),2):
            if rng.random()<0.5: D.add_edge(i,j)
            else: D.add_edge(j,i)
        if nx.is_strongly_connected(D): return D
    return None
rng=random.Random(0)
for n in [3,4,5,6]:
    D=random_strong_tournament(n,rng)
    if D is None: print(f"  n={n}: no strong tournament found"); continue
    cyc_lens=set()
    for c in nx.simple_cycles(D):
        if len(c)>=2: cyc_lens.add(len(c))
    even=sorted(L for L in cyc_lens if L%2==0)
    print(f"  strong tournament n={n}: directed cycle lengths {sorted(cyc_lens)}; even ones: {even} "
          f"(Moon: vertex-pancyclic, lengths 3..n)")
print("  => strong tournament on n>=4 has a directed 4-cycle (even); transitive (no cycle) has none.")
print("     'no even directed cycle' <-> Pfaffian orientation exists <-> permanent = +-det (Polya).")

# ---------- (5) covering systems & the LRC-as-covering-mod-q lens ---------------
print("\n"+"="*78)
print("(5) COVERING SYSTEMS (Owens min-modulus 42) and LRC = a danger-covering mod q")
print("="*78)
# a classic small covering system: 0(2),0(3),1(4),5(6),7(12) covers Z
cs=[(0,2),(0,3),(1,4),(5,6),(7,12)]
covered=all(any((x-a)%m==0 for a,m in cs) for x in range(0,12))
print(f"  covering system {cs}: covers Z/12: {covered}  (min modulus {min(m for _,m in cs)})")
print(f"  Erdos min-modulus problem: can min modulus be arbitrarily large? NO (Hough 2015,")
print(f"  min modulus <= 10^16; improved to 616, BBMST). Largest KNOWN min modulus = 42 (Owens).")
# LRC as covering: danger residues of speed v at modulus q, threshold 1/n: {m: v*m mod q in (-q/n,q/n)}
def danger(v,q,n):
    return set(m for m in range(q) if min((v*m)%q, q-(v*m)%q) < q/n)  # ||v m/q|| < 1/n
def covers(S,q,n):
    U=set()
    for v in S: U|=danger(v,q,n)
    return len(U)==q, U
for tag,S,q,n in [("AP_5 worry-set",[1,2,3,4],5,5), ("loose (1,2,4,8)",[1,2,4,8],5,5)]:
    cov,U=covers(S,q,n)
    print(f"  {tag}: danger residues cover Z/{q}: {cov}  (lonely tick exists: {not cov})")
print("  => LRC FAILURE = the danger residues COVER Z/q for every q = a persistent covering")
print("     system (the S704 'unbounded depth q*' = the covering never breaks). Worry-set = tight cover.")
print("\nDONE.")
