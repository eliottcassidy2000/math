"""Timing probe: cubes = colourings of the 7 consecutive neighbours zeta^0..zeta^6 of 0 (colour(0)=0), canonical under S_4."""
import sys, time, random
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
from pysat.solvers import Solver
q=13; k=5
adj,p,e,mu,elems,idx,add=gq(q); n=len(adj)
elems2,idx2,mul,add2,pw=field(p,2*e)
# a generator zeta of mu_14: element of order 14
one=tuple([1]+[0]*(2*e-1))
zeta=next(u for u in mu if all(pw(u,d)!=one for d in (1,2,7)))
path=[idx[pw(zeta,i)] for i in range(7)]
assert all(path[j] not in adj[path[i]] for i in range(7) for j in range(7))
var=lambda v,col: v*k+col+1
base=[]
for v in range(n):
    base.append([var(v,col) for col in range(k)])
    for col in range(k):
        for col2 in range(col+1,k): base.append([-var(v,col),-var(v,col2)])
for v in range(n):
    for u_ in adj[v]:
        if u_>v:
            for col in range(k): base.append([-var(v,col),-var(u_,col)])
base.append([var(0,0)])
solver=sys.argv[1] if len(sys.argv)>1 else 'cadical153'
random.seed(7)
for trial in range(int(sys.argv[2]) if len(sys.argv)>2 else 3):
    # random restricted-growth colouring of the 7 path vertices with colours 1..4
    cols=[1]
    for i in range(1,7): cols.append(random.randint(1,min(4,max(cols)+1)))
    s=Solver(name=solver, bootstrap_with=base)
    for v,c_ in zip(path,cols): s.add_clause([var(v,c_)])
    t=time.time(); r=s.solve(); dt=time.time()-t; s.delete()
    print("solver %s cube %s -> %s in %.1fs"%(solver,''.join(map(str,cols)), "SAT" if r else "UNSAT", dt)); sys.stdout.flush()
