"""Decision version: does G_q have an independent set of size >= t containing 0?  (CP-SAT feasibility.)
Usage: python3 gq_indep_decision.py q t [workers] [time_limit]"""
import sys, time
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
from ortools.sat.python import cp_model
q=int(sys.argv[1]); t=int(sys.argv[2]); workers=int(sys.argv[3]) if len(sys.argv)>3 else 1; tl=float(sys.argv[4]) if len(sys.argv)>4 else 7200.0
adj,p,e,mu,elems,idx,add=gq(q); n=len(adj)
m=cp_model.CpModel(); x=[m.NewBoolVar("x%d"%i) for i in range(n)]
for i in range(n):
    for j in adj[i]:
        if j>i: m.AddBoolOr([x[i].Not(), x[j].Not()])
m.Add(x[0]==1); m.Add(sum(x)>=t)
s=cp_model.CpSolver(); s.parameters.num_search_workers=workers; s.parameters.max_time_in_seconds=tl
t0=time.time(); st=s.Solve(m)
print("q=%d n=%d size>=%d: %s (%.1fs)"%(q,n,t,s.StatusName(st),time.time()-t0))
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    S=[i for i in range(n) if s.Value(x[i])]; assert all(j not in adj[i] for i in S for j in S); print("set of size %d:"%len(S), S)
