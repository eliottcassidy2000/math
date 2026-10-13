"""kappa(q) <= 5?  Largest-class decomposition: a 5-colouring of G_q (vertex-transitive) has a colour class of size
>= ceil(q^2/5); translate it to contain 0.  So G_q is 5-colourable iff some independent set I containing 0 with
|I| >= ceil(q^2/5) has G_q - I 4-colourable.  Phase 1 enumerates such I (CP-SAT, symmetry-reduced by the stabiliser
mu_{q+1} x| Gal of 0: we only keep canonical representatives), phase 2 tests 4-colourability of G - I by CaDiCaL.
Usage: python3 kappa_largest_class.py q [max_sets] [workers]"""
import sys, time, itertools
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
from ortools.sat.python import cp_model
from pysat.solvers import Solver
q=int(sys.argv[1]); max_sets=int(sys.argv[2]) if len(sys.argv)>2 else 10**6; workers=int(sys.argv[3]) if len(sys.argv)>3 else 4
adj,p,e,mu,elems,idx,add=gq(q); n=len(adj); thr=-(-n//5)
# stabiliser of 0: z -> a * sigma(z), a in mu, sigma Frobenius powers
def frob(t,k):  # t^(p^k)
    r=t
    for _ in range(k): r=pw(r,p)
    return r
# field ops from gq's closure are not exported; rebuild
elems2,idx2,mul,add2,pw=field(p,2*e)
assert elems2==elems
stab=[]
for a in mu:
    for k in range(2*e):
        stab.append([idx[mul(a,frob(t,k))] for t in elems])
print("q=%d n=%d threshold=%d |Stab_0|=%d"%(q,n,thr,len(stab)))
def canon(S):
    best=None
    for g in stab:
        img=tuple(sorted(g[v] for v in S))
        if best is None or img<best: best=img
    return best
# phase 1: enumerate independent sets containing 0 of size >= thr (all, then reduce by stabiliser)
m=cp_model.CpModel(); x=[m.NewBoolVar("x%d"%i) for i in range(n)]
for i in range(n):
    for j in adj[i]:
        if j>i: m.AddBoolOr([x[i].Not(), x[j].Not()])
m.Add(x[0]==1); m.Add(sum(x)>=thr)
class Coll(cp_model.CpSolverSolutionCallback):
    def __init__(s): super().__init__(); s.sets=set(); s.count=0
    def on_solution_callback(s):
        S=[i for i in range(n) if s.Value(x[i])]; s.count+=1
        s.sets.add(canon(S))
        if s.count>=max_sets: s.StopSearch()
s=cp_model.CpSolver(); s.parameters.enumerate_all_solutions=True; s.parameters.num_search_workers=1
c=Coll(); t=time.time(); st=s.Solve(m,c)
print("phase 1: %d independent sets >= %d containing 0 enumerated (%s), %d canonical classes, %.1fs"%(c.count,thr,s.StatusName(st),len(c.sets),time.time()-t)); sys.stdout.flush()
# phase 2: for each canonical I, is G - I 4-colourable?
found=False; t=time.time()
for k,I in enumerate(sorted(c.sets)):
    Iset=set(I); rest=[v for v in range(n) if v not in Iset]; pos={v:i for i,v in enumerate(rest)}
    var=lambda v,col: pos[v]*4+col+1
    sol=Solver(name='cadical153')
    for v in rest:
        sol.add_clause([var(v,col) for col in range(4)])
    for v in rest:
        for u_ in adj[v]:
            if u_ in pos and u_>v:
                for col in range(4): sol.add_clause([-var(v,col),-var(u_,col)])
    r=sol.solve(); sol.delete()
    if r:
        print("5-COLOURABLE: largest class", I); found=True; break
    if k%50==0: print("  tested %d classes (%.1fs)"%(k+1,time.time()-t)); sys.stdout.flush()
print("RESULT q=%d: %s (%d classes tested, %.1fs)"%(q, "5-colourable" if found else "NOT 5-colourable (every candidate largest class leaves a non-4-colourable remainder)", k+1 if c.sets else 0, time.time()-t))
