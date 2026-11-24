#!/usr/bin/env python3
"""Two exhaustive seed classes, plus an optional bounded SciPy/HiGHS search."""
import argparse
from collections import Counter
from pathlib import Path
from runpy import run_path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--class-index", type=int, choices=[0,1])
p.add_argument("--seconds", type=float, default=180)
p.add_argument("--witness", default="gewirtz-cover.npy")
a = p.parse_args()
m = run_path(str(Path(__file__).with_name("gewirtz-model.py")))
gens = m["G"]
G = m["group"](gens)
linear = set(m["group"](gens[:-1]))
I = tuple(range(56))
adj = [{j for j,y in enumerate(m["H"]) if not x&y} for x in m["H"]]
D = [g for g in G if m["derange"](g)]

def moved(g):
    return sum(j in adj[i] for i,j in enumerate(g))

assert len(G)==80640 and len(linear)==40320
assert Counter((g not in linear,moved(g)) for g in D)=={
    (True,8):16800, (False,14):11520, (False,8):7560, (True,20):5040}
pool = {g for g in D if g in linear and moved(g)==8}
X = [(h,m["inverse"](h)) for h in gens]
classes = []
while pool:
    g = min(pool)
    orbit = [g]
    pool.remove(g)
    for v in orbit:
        for h,hi in X:
            w = tuple(h[v[i]] for i in hi)
            if w in pool:
                pool.remove(w)
                orbit.append(w)
    classes.append((len(orbit),g))
classes.sort()
assert [size for size,g in classes]==[2520,5040]
data = []
for k,(size,seed) in enumerate(classes):
    candidates = [g for g in D if all(x!=y for x,y in zip(g,seed))]
    print(f"seed {k}: class {size}, {len(candidates)} candidates after identity and fixed seed.",flush=True)
    data.append((seed,candidates))
if a.class_index is None:
    raise SystemExit(0)

import numpy as np
from scipy.sparse import coo_matrix
from scipy.optimize import milp, Bounds, LinearConstraint
seed,candidates = data[a.class_index]
C=np.array(candidates,dtype=np.uint8)
N=56
position=np.arange(N)*N
holes=np.ones(N*N,dtype=bool)
holes[position+I]=holes[position+seed]=False
assert holes.sum()==3024
r=(C+position).ravel()
j=np.repeat(np.arange(len(C)),N)
M=coo_matrix((np.ones(len(j)),(r,j)),shape=(N*N,len(C))).tocsr()[holes].tocsc()
sign=np.array([g not in linear for g in candidates],dtype=int)
# Since both fixed maps are collineations, the original correlation bound is unchanged.
res=milp(c=np.zeros(len(C)), integrality=np.ones(len(C)),
         bounds=Bounds(0,1),
         constraints=[LinearConstraint(M,1,1),LinearConstraint(sign[None,:],7,49)],
         options={"time_limit":a.seconds,"disp":True})
print(res.message,flush=True)
if res.x is not None:
    S=np.vstack([I,seed,C[res.x>.5]])
    assert S.shape==(56,56) and np.all(np.sort(S,axis=0)==np.arange(56)[:,None])
    np.save(a.witness,S)
    print(f"Verified 56-map witness saved: {a.witness}",flush=True)
