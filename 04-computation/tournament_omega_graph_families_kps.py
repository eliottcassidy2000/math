"""
tournament_omega_graph_families_kps.py  (kind-pasteur-2026-06-27-S31ah)

Generalize THM-201/202 (K3, P4 forbidden as Omega) to GRAPH FAMILIES: which
structural classes can be a tournament's conflict graph Omega(T)?
Realized small Omega are clique-like (HYP-3101); test which families appear:
  clique, cycle C_m, star K_{1,m}, path P_m, tree (acyclic), bipartite,
  disconnected. Forbidden families = new certificate classes.

Strategy: near-transitive tournaments (transitive base + k random flips) have
SMALL Omega; enumerate/sample them at n<=8, classify each connected Omega.
"""
import sys, random, itertools
from tournament_certificate_engine_kps import conflict_graph, independence_poly

def classify_graph(m, E):
    """structural labels for a small connected graph (m vertices, adj E)."""
    if m==0: return "empty"
    if m==1: return "K1"
    deg=[sum(1 for w in range(m) if E[v][w]) for v in range(m)]
    edges=sum(deg)//2
    full=m*(m-1)//2
    # connected?
    seen=[False]*m; st=[0]; seen[0]=True; cc=1
    while st:
        v=st.pop()
        for w in range(m):
            if E[v][w] and not seen[w]: seen[w]=True; st.append(w);
    connected=all(seen)
    labels=[]
    if not connected: return "DISCONNECTED"
    if edges==full: labels.append(f"CLIQUE K{m}")
    if edges==m and all(d==2 for d in deg): labels.append(f"CYCLE C{m}")
    if edges==m-1:
        labels.append("TREE")
        if sorted(deg)==[1]*(m-1)+[m-1]: labels.append(f"STAR K1,{m-1}")
        if sorted(deg)==[1,1]+[2]*(m-2): labels.append(f"PATH P{m}")
    # bipartite check (2-coloring)
    color=[-1]*m; bip=True; color[0]=0; q=[0]
    while q and bip:
        v=q.pop()
        for w in range(m):
            if E[v][w]:
                if color[w]==-1: color[w]=1-color[v]; q.append(w)
                elif color[w]==color[v]: bip=False
    if bip and edges>0: labels.append("BIPARTITE")
    return ", ".join(labels) if labels else f"{m}v-{edges}e (other)"

def near_transitive(n, k, rng):
    """transitive tournament + k random non-base flips."""
    adj=[[1 if i<j else 0 for j in range(n)] for i in range(n)]
    pairs=list(itertools.combinations(range(n),2))
    for (i,j) in rng.sample(pairs, min(k,len(pairs))):
        adj[i][j]=1-adj[i][j]; adj[j][i]=1-adj[j][i]
    return adj

def omega_connected_class(adj):
    m,E=conflict_graph(adj)
    if m==0: return None
    # only classify if connected (single component)
    seen=[False]*m; st=[0]; seen[0]=True
    while st:
        v=st.pop()
        for w in range(m):
            if E[v][w] and not seen[w]: seen[w]=True; st.append(w)
    if not all(seen): return None
    a=independence_poly(m,E); H=sum(a[k]*2**k for k in range(len(a)))
    return (m, classify_graph(m,E), H)

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(202)
    realized={}  # (m, class) -> example H
    family_count={}
    for n in range(4,9):
        for k in range(1,6):
            for _ in range(4000):
                adj=near_transitive(n,k,random)
                r=omega_connected_class(adj)
                if r and r[0]<=6:
                    m,cls,H=r
                    realized[(m,cls)]=H
                    for tok in cls.split(", "):
                        base=tok.split()[0]
                        family_count[base]=family_count.get(base,0)+1
    print("Connected Omega structural classes REALIZED (m<=6), with example H:")
    for (m,cls),H in sorted(realized.items()):
        print(f"  m={m}: {cls:30s} H={H}")
    print(f"\nFamily realized-counts: {dict(sorted(family_count.items(), key=lambda x:-x[1]))}")
    print("\nForbidden families (NEVER realized as connected Omega in this scan):")
    for fam in ["CYCLE","STAR","PATH","TREE","BIPARTITE"]:
        seen=any(fam in k[1] for k in realized)
        print(f"  {fam}: {'realized' if seen else 'NOT realized (candidate forbidden family)'}")
