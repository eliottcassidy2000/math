"""k-colourability of G_q by CP-SAT with symmetry detection, logging.  Usage: q k workers time_limit"""
import sys, time
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
from ortools.sat.python import cp_model
q=int(sys.argv[1]); k=int(sys.argv[2]); workers=int(sys.argv[3]); tl=float(sys.argv[4])
adj,p,e,mu,elems,idx,add=gq(q); n=len(adj)
m=cp_model.CpModel(); c=[m.NewIntVar(0,k-1,"c%d"%i) for i in range(n)]
for i in range(n):
    for j in adj[i]:
        if j>i: m.Add(c[i]!=c[j])
m.Add(c[0]==0); one=idx[tuple([1]+[0]*(2*e-1))]; m.Add(c[one]==1)
s=cp_model.CpSolver(); s.parameters.num_search_workers=workers; s.parameters.max_time_in_seconds=tl
s.parameters.log_search_progress=True; s.parameters.symmetry_level=2
t=time.time(); st=s.Solve(m)
print("RESULT q=%d k=%d: %s (%.1fs)"%(q,k,s.StatusName(st),time.time()-t))
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    col=[s.Value(c[i]) for i in range(n)]; assert all(col[i]!=col[j] for i in range(n) for j in adj[i]); print("colouring:", col)
