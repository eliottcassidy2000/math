"""Phase 1 of the largest-class decomposition for G_13 and k = 5: enumerate the MAXIMAL independent sets I with 0 in I and
|I| >= 34 (= ceil(169/5)); reduce by the stabiliser of 0 (order 28) to canonical representatives; write them to a file.
(A 5-colouring can be translated so that a largest colour class contains 0, and that class can be enlarged to a maximal
independent set by moving vertices into it, staying proper.)"""
import sys, time
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
from ortools.sat.python import cp_model
q=13; adj,p,e,mu,elems,idx,add=gq(q); n=len(adj); thr=34
elems2,idx2,mul,add2,pw=field(p,2*e)
def frob(t,k):
    r=t
    for _ in range(k): r=pw(r,p)
    return r
stab=[[idx[mul(a,frob(t,k))] for t in elems] for a in mu for k in range(2*e)]
def canon(S):
    return min(tuple(sorted(g[v] for v in S)) for g in stab)
m=cp_model.CpModel(); x=[m.NewBoolVar("x%d"%i) for i in range(n)]
for i in range(n):
    for j in adj[i]:
        if j>i: m.AddBoolOr([x[i].Not(), x[j].Not()])
for v in range(n):  # maximality: v in I or some neighbour in I
    m.AddBoolOr([x[v]]+[x[u] for u in adj[v]])
m.Add(x[0]==1); m.Add(sum(x)>=thr)
class Coll(cp_model.CpSolverSolutionCallback):
    def __init__(s): super().__init__(); s.sets={}; s.count=0; s.t=time.time()
    def on_solution_callback(s):
        S=[i for i in range(n) if s.Value(x[i])]; s.count+=1
        c=canon(S); s.sets[c]=len(S)
        if s.count%2000==0:
            print("  %d raw solutions, %d canonical classes, sizes %s, %.0fs"%(s.count,len(s.sets),sorted(set(s.sets.values())),time.time()-s.t)); sys.stdout.flush()
sol=cp_model.CpSolver(); sol.parameters.enumerate_all_solutions=True; sol.parameters.num_search_workers=1
c=Coll(); t=time.time(); st=sol.Solve(m,c)
print("DONE: %d raw maximal independent sets of size >= %d containing 0 (%s), %d canonical classes, sizes histogram %s, %.0fs"%(
    c.count,thr,sol.StatusName(st),len(c.sets),{s_:sum(1 for v in c.sets.values() if v==s_) for s_ in sorted(set(c.sets.values()))},time.time()-t))
with open("kappa13_mis_classes.txt","w") as f:
    for S in sorted(c.sets): f.write(" ".join(map(str,S))+"\n")
