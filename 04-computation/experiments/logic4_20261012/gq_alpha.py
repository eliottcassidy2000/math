"""Independence number alpha(G_q) of the unit-circle graph G_q = Cay(F_{q^2}, mu_{q+1}) by CP-SAT, with vertex 0 forced into the set
(vertex-transitive).  If alpha(G_q) < q^2/k then G_q is not k-colourable.  Usage: python3 gq_alpha.py q [workers] [time_limit_s]"""
import sys, time
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
from ortools.sat.python import cp_model
q=int(sys.argv[1]); workers=int(sys.argv[2]) if len(sys.argv)>2 else 4; tl=float(sys.argv[3]) if len(sys.argv)>3 else 3600.0
adj,p,e,mu,elems,idx,add=gq(q); n=len(adj)
m=cp_model.CpModel(); x=[m.NewBoolVar("x%d"%i) for i in range(n)]
for i in range(n):
    for j in adj[i]:
        if j>i: m.AddBoolOr([x[i].Not(), x[j].Not()])
m.Add(x[0]==1)
m.Maximize(sum(x))
s=cp_model.CpSolver(); s.parameters.num_search_workers=workers; s.parameters.max_time_in_seconds=tl
t=time.time(); st=s.Solve(m)
print("q=%d n=%d: status=%s alpha=%s bound=%s time=%.1fs  (n/5=%.2f, n/6=%.2f)"%(q,n,s.StatusName(st),int(s.ObjectiveValue()),int(s.BestObjectiveBound()),time.time()-t,n/5,n/6))
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    S=[i for i in range(n) if s.Value(x[i])]
    assert all(j not in adj[i] for i in S for j in S)
    print("independent set:", S)
