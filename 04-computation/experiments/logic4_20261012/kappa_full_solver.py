"""Full 5-colourability instance of G_q with a chosen pysat solver, symmetry breaking colour(0)=0, colour(1)=1.  Usage: q k solver"""
import sys, time
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
from pysat.solvers import Solver
q=int(sys.argv[1]); k=int(sys.argv[2]); name=sys.argv[3]
adj,p,e,mu,elems,idx,add=gq(q); n=len(adj); one=idx[tuple([1]+[0]*(2*e-1))]
var=lambda v,col: v*k+col+1
s=Solver(name=name)
for v in range(n):
    s.add_clause([var(v,col) for col in range(k)])
    for col in range(k):
        for col2 in range(col+1,k): s.add_clause([-var(v,col),-var(v,col2)])
for v in range(n):
    for u_ in adj[v]:
        if u_>v:
            for col in range(k): s.add_clause([-var(v,col),-var(u_,col)])
s.add_clause([var(0,0)]); s.add_clause([var(one,1)])
t=time.time(); r=s.solve(); dt=time.time()-t
print("q=%d k=%d solver=%s: %s (%.1fs)"%(q,k,name,"SAT" if r else "UNSAT",dt))
if r:
    m=s.get_model(); col={v:[c for c in range(k) if m[var(v,c)-1]>0][0] for v in range(n)}
    assert all(col[v]!=col[u_] for v in range(n) for u_ in adj[v]); print("colouring:", [col[v] for v in range(n)])
