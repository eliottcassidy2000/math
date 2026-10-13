"""alpha(G_q) by CP-SAT with progress logging (bounds) and symmetry level 2.  Usage: q workers time_limit"""
import sys, time
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
from ortools.sat.python import cp_model
q=int(sys.argv[1]); workers=int(sys.argv[2]); tl=float(sys.argv[3])
adj,p,e,mu,elems,idx,add=gq(q); n=len(adj)
m=cp_model.CpModel(); x=[m.NewBoolVar("x%d"%i) for i in range(n)]
for i in range(n):
    for j in adj[i]:
        if j>i: m.AddBoolOr([x[i].Not(), x[j].Not()])
m.Add(x[0]==1)
# redundant: every 13-cycle line in a square-norm direction contains at most floor(q/2) points of an independent set
elems2,idx2,mul,add2,pw=field(p,2*e)
for u in mu:
    line=[idx[mul(elems[i],u)] for i in range(n) if pw(elems[i],q)==elems[i]]
    seen=set()
    for i in range(n):
        if i in seen: continue
        cs=[idx[add(elems[i],elems[j])] for j in line]; seen.update(cs)
        m.Add(sum(x[v] for v in cs)<=q//2)
m.Maximize(sum(x))
s=cp_model.CpSolver(); s.parameters.num_search_workers=workers; s.parameters.max_time_in_seconds=tl
s.parameters.log_search_progress=True; s.parameters.symmetry_level=2
t=time.time(); st=s.Solve(m)
print("RESULT q=%d: status=%s alpha=%s bound=%s time=%.1fs"%(q,s.StatusName(st),int(s.ObjectiveValue()),int(s.BestObjectiveBound()),time.time()-t))
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    S=[i for i in range(n) if s.Value(x[i])]; assert all(j not in adj[i] for i in S for j in S); print("independent set:", S)
