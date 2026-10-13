"""Timing test for cube-and-conquer on 5-colourability of G_13: fix colour(0)=0 and random 4-colourings of the 14
neighbours of 0 (an independent set), solve the remainder with CaDiCaL."""
import sys, time, random
exec(open("gq_automorphisms.py").read().split("qs=[q for q")[0])
from pysat.solvers import Solver
q=13; k=5
adj,p,e,mu,elems,idx,add=gq(q); n=len(adj)
nb=sorted(adj[0]); assert len(nb)==14 and all(v not in adj[u] for u in nb for v in nb)
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
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
times=[]
for trial in range(int(sys.argv[2]) if len(sys.argv)>2 else 10):
    cols=[random.randint(1,4) for _ in nb]
    s=Solver(name='cadical153', bootstrap_with=base)
    for v,c_ in zip(nb,cols): s.add_clause([var(v,c_)])
    t=time.time(); r=s.solve(); dt=time.time()-t; s.delete()
    times.append(dt); print("cube %s -> %s in %.2fs"%(''.join(map(str,cols)), "SAT" if r else "UNSAT", dt)); sys.stdout.flush()
print("mean %.2fs max %.2fs"%(sum(times)/len(times), max(times)))
