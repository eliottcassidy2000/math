"""kappa(q) = chi(Cay(F_{q^2}, mu_{q+1})): decide 5-colourability of G_13 and G_17 with CaDiCaL (pysat), symmetry breaking
colour(0)=0, colour(1)=1 (the affine group is vertex-transitive and mu_{q+1} is transitive on the neighbours of 0: THM-4624)."""
import sys, time
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
from pysat.solvers import Solver
q=int(sys.argv[1]); k=int(sys.argv[2]); solver_name=sys.argv[3] if len(sys.argv)>3 else 'cadical153'
adj,p,e,mu,elems,idx,add=gq(q)
n=len(adj); one=idx[tuple([1]+[0]*(2*e-1))]
var=lambda v,col: v*k+col+1
s=Solver(name=solver_name)
for v in range(n):
    s.add_clause([var(v,col) for col in range(k)])
    for col in range(k):
        for col2 in range(col+1,k): s.add_clause([-var(v,col),-var(v,col2)])
for v in range(n):
    for u_ in adj[v]:
        if u_>v:
            for col in range(k): s.add_clause([-var(v,col),-var(u_,col)])
s.add_clause([var(0,0)]); s.add_clause([var(one,1)])
t=time.time(); res=s.solve(); dt=time.time()-t
print("q=%d k=%d solver=%s: %s (%.1fs)"%(q,k,solver_name,"SAT (k-colourable)" if res else "UNSAT (not k-colourable)",dt))
if res:
    m=s.get_model(); col={v:[c_ for c_ in range(k) if m[var(v,c_)-1]>0][0] for v in range(n)}
    assert all(col[v]!=col[u_] for v in range(n) for u_ in adj[v]); print("colouring verified")
