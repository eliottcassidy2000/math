"""Look for a small non-5-colourable induced subgraph of G_13: balls around 0 and unions of norm circles.
Each test uses kissat404 with a conflict budget."""
import sys, time, itertools
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
from pysat.solvers import Solver
q=13; k=5
adj,p,e,mu,elems,idx,add=gq(q); n=len(adj)
elems2,idx2,mul,add2,pw=field(p,2*e)
norm=[pw(z,q+1)[0] for z in elems]
def test(V, budget=2_000_000, name='kissat404'):
    V=sorted(V); pos={v:i for i,v in enumerate(V)}; var=lambda v,col: pos[v]*k+col+1
    s=Solver(name=name)
    for v in V:
        s.add_clause([var(v,col) for col in range(k)])
        for col in range(k):
            for col2 in range(col+1,k): s.add_clause([-var(v,col),-var(v,col2)])
    for v in V:
        for u_ in adj[v]:
            if u_ in pos and u_>v:
                for col in range(k): s.add_clause([-var(v,col),-var(u_,col)])
    if 0 in pos: s.add_clause([var(0,0)])
    s.conf_budget(budget); t=time.time(); r=s.solve_limited(); dt=time.time()-t; s.delete()
    return ("SAT" if r else "UNSAT" if r is False else "UNKNOWN"), dt
# balls
B={0}; frontier={0}
for r in range(1,4):
    frontier={u_ for v in frontier for u_ in adj[v]}-B; B|=frontier
    print("ball radius %d: %d vertices"%(r,len(B))); sys.stdout.flush()
    if len(B)<n:
        res,dt=test(B); print("   5-colourable? %s (%.1fs)"%(res,dt)); sys.stdout.flush()
# circle unions: vertex 0 plus circles with norms in T
circles={m:[i for i in range(n) if norm[i]==m] for m in range(1,q)}
for size in range(6,12):
    for T in itertools.combinations(range(1,q),size):
        V={0}|{v for m in T for v in circles[m]}
        res,dt=test(V, budget=300_000)
        if res!="SAT":
            print("circles T=%s (%d vertices): %s (%.1fs)"%(T,len(V),res,dt)); sys.stdout.flush()
    print("size %d done"%size); sys.stdout.flush()
