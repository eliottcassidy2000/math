"""The three exceptional unit-circle graphs: G_4 (Clebsch), G_8 (Hamming H(3,4)), G_9 (ternary 5-frame).
Checks: explicit isomorphisms via subfield lines, chromatic numbers by SAT, and the frame colouring of G_9."""
import itertools, networkx as nx
from pysat.solvers import Glucose3
import importlib.util, sys
spec=importlib.util.spec_from_file_location("gqa","gq_automorphisms.py")
# reuse field/gq constructors without running the census: copy minimal code
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
def chromatic(adj):
    n=len(adj)
    for k in range(2,8):
        s=Glucose3(); var=lambda v,c: v*k+c+1
        for v in range(n): s.add_clause([var(v,c) for c in range(k)])
        for v in range(n):
            for w_ in adj[v]:
                if w_>v:
                    for c in range(k): s.add_clause([-var(v,c),-var(w_,c)])
        ok=s.solve(); s.delete()
        if ok: return k
for q in (4,8,9):
    adj,p,e,mu,elems,idx,add=gq(q)
    G=nx.Graph([(i,j) for i in adj for j in adj[i] if j>i])
    print("q=%d: |V|=%d, degree=%d, chi=%d"%(q,len(adj),len(mu),chromatic(adj)))
    if q==4:
        # Clebsch = folded 5-cube: vertices F_2^4, generators e1..e4 and e1+e2+e3+e4
        H=nx.Graph(); vs=list(itertools.product(range(2),repeat=4)); gens=[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(1,1,1,1)]
        for v in vs:
            for g in gens: H.add_edge(v,tuple((a+b)%2 for a,b in zip(v,g)))
        print("   G_4 isomorphic to Cay(F_2^4, {e1,e2,e3,e4,e1+e2+e3+e4}) (Clebsch):", nx.is_isomorphic(G,H))
        # mu_5 sums to zero
        tot=elems[0]
        for m_ in mu: tot=add(tot,m_)
        print("   sum of mu_5 in F_16 =", tot)
    if q==8:
        H=nx.cartesian_product(nx.cartesian_product(nx.complete_graph(4),nx.complete_graph(4)),nx.complete_graph(4))
        print("   G_8 isomorphic to K4 x K4 x K4 = H(3,4):", nx.is_isomorphic(G,H))
    if q==9:
        # ternary 5-frame: F_3^4 with generators +-w_j, w_0..w_3 = basis, w_4 = -(w_0+..+w_3)
        vs=list(itertools.product(range(3),repeat=4)); basis=[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]; w4=(2,2,2,2)
        gens=basis+[w4]; gens=gens+[tuple((-g_)%3 for g_ in g) for g in gens]
        H=nx.Graph()
        for v in vs:
            for g in gens: H.add_edge(v,tuple((a+b)%3 for a,b in zip(v,g)))
        print("   G_9 isomorphic to Cay(F_3^4, {+-w_0..+-w_4}, sum w_j = 0):", nx.is_isomorphic(G,H))
        # additive 3-colouring: lambda(w_0..w_3) = (1,1,1,1) => lambda(w_4) = -4 = -1: nowhere zero
        lam=lambda v: (v[0]+v[1]+v[2]+v[3])%3
        print("   lambda nowhere zero on generators:", all(lam(g)!=0 for g in gens), "=> proper 3-colouring of H by lambda:", all(lam(v)!=lam(w_) for v,w_ in H.edges()))
