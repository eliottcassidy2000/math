# Six private interiors: classify 10-vertex 7-critical quotients (4 port groups + 6 interiors) surviving
# planar unit-distance necessary conditions.  Model = FIVE_INTERIORS.md (paired_budget_selectors_matrix_20261028).
import subprocess, itertools, networkx as nx, sys
from functools import lru_cache
g6=subprocess.run(['geng','-q','10','-D3'],capture_output=True,text=True).stdout.split()
def chrom_le(G,k):
    # exact k-colourability by backtracking (10 vertices)
    nodes=sorted(G.nodes, key=lambda v:-G.degree(v)); col={}
    def bt(i):
        if i==len(nodes): return True
        v=nodes[i]; used={col[u] for u in G[v] if u in col}
        for c in range(k):
            if c not in used:
                col[v]=c
                if bt(i+1): return True
                del col[v]
        return False
    return bt(0)
crit=[]
for s in g6:
    Gam=nx.from_graph6_bytes(s.encode()); G=nx.complement(Gam)
    if chrom_le(G,6): continue                     # need chi(G) >= 7
    if all(chrom_le(G.subgraph([u for u in G if u!=v]).copy(),6) for v in G):   # vertex-critical
        crit.append(G)
print("10-vertex vertex-7-critical graphs:",len(crit))
for G in crit:
    Gam=nx.complement(G)
    print("   complement components:",sorted([len(c) for c in nx.connected_components(Gam)]),
          "complement edges",Gam.number_of_edges(), "degree seq of complement", sorted(d for _,d in Gam.degree()))
# --- port/interior split and geometric filters ---
def has_K4(H): return any(all(H.has_edge(a,b) for a,b in itertools.combinations(Q,2)) for Q in itertools.combinations(H.nodes,4))
def has_K23(H):
    for a,b in itertools.combinations(H.nodes,2):
        if len(set(H[a])&set(H[b]))>=3: return True
    return False
def circle_ok(H,T):
    # points on one unit circle: induced unit graph is a linear forest, or exactly a 6-cycle (regular hexagon)
    if len(T)<=1: return True
    S=H.subgraph(T)
    if len(T)==6 and S.number_of_edges()==6 and all(d==2 for _,d in S.degree()) and nx.is_connected(S): return True
    return nx.is_forest(S) and all(d<=2 for _,d in S.degree())
survivors=[]
for gi,G in enumerate(crit):
    for P in itertools.combinations(G.nodes,4):
        if not all(G.has_edge(a,b) for a,b in itertools.combinations(P,2)): continue
        I=[v for v in G if v not in P]; H=G.subgraph(I).copy()
        if has_K4(H) or has_K23(H): continue
        S={g:[x for x in I if G.has_edge(g,x)] for g in P}
        # each group: choose two point-contact sets T1,T2 with union S[g]
        opts={}
        for g in P:
            o=[]
            Sg=S[g]
            for assign in itertools.product((1,2,3),repeat=len(Sg)):   # 1: point1 only, 2: point2 only, 3: both
                T1=[x for x,a in zip(Sg,assign) if a in (1,3)]; T2=[x for x,a in zip(Sg,assign) if a in (2,3)]
                if tuple(sorted(T2))<tuple(sorted(T1)): continue
                if circle_ok(H,T1) and circle_ok(H,T2) and len(set(T1)&set(T2))<=2: o.append((frozenset(T1),frozenset(T2)))
            opts[g]=o
        if any(not opts[g] for g in P): continue
        ok=False
        for choice in itertools.product(*[opts[g] for g in P]):
            pts=[T for pair in choice for T in pair]            # eight port points' contact sets
            bad=False
            for a,b in itertools.combinations(range(8),2):
                if len(pts[a]&pts[b])>2: bad=True; break
            if bad: continue
            for p in pts:                                       # interior x vs port point p: common interior nbrs <= 2
                for x in I:
                    if len(set(H[x])&p)>2: bad=True; break
                if bad: break
            if bad: continue
            for x,y in itertools.combinations(I,2):             # two interiors: common nbrs (interiors + port points) <= 2
                c=len(set(H[x])&set(H[y]))+sum(1 for p in pts if x in p and y in p)
                if c>2: bad=True; break
            if bad: continue
            ok=True; break
        if ok: survivors.append((gi,P))
print("(critical graph, port 4-clique) pairs surviving the planar filters:",len(survivors))
seen=set()
for gi,P in survivors:
    G=crit[gi]; I=[v for v in G if v not in P]; H=G.subgraph(I)
    key=(gi,)
    if key in seen: continue
    seen.add(key)
    print(f"  graph {gi}: interiors induce {H.number_of_edges()} unit edges, degrees {sorted(d for _,d in H.degree())}, contacts per group {[len([x for x in I if G.has_edge(g,x)]) for g in P]}")
