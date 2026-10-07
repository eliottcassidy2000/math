#!/usr/bin/env python3
"""CEGAR for the fully translation-invariant tree-grid game Q_inv(n,t) (THM-470's Finv):
E = set of lex-positive gap vectors in Z^n (|coords|<=t-1); triangle-free on [t]^n <=> no realizable
d1,d2,d1+d2 in E; must hit every binary subgrid of [t]^n.  Complete verifier: enumerates every
placement of the two height-(n-1) row subgrids in [t]^(n-1) and every root gap h.
Usage: finv_cegar.py n t [timelimit_s] [batch]
"""
import sys, itertools, time
import numpy as np
from pysat.solvers import Solver

n=int(sys.argv[1]); t=int(sys.argv[2])
TL=float(sys.argv[3]) if len(sys.argv)>3 else 3600
BATCH=int(sys.argv[4]) if len(sys.argv)>4 else 64
t0=time.time()
R=range(-(t-1),t)
def lexpos(d):
    for x in d:
        if x>0: return True
        if x<0: return False
    return False
gaps=[d for d in itertools.product(R,repeat=n) if lexpos(d)]
gid={d:i+1 for i,d in enumerate(gaps)}
print(f"n={n} t={t}: {len(gaps)} gap classes",flush=True)
s=Solver(name='cadical195')
ntri=0
for d1 in gaps:
    for d2 in gaps:
        ok=True
        for i in range(n):
            lo=min(0,d1[i],d1[i]+d2[i]); hi=max(0,d1[i],d1[i]+d2[i])
            if hi-lo>t-1: ok=False;break
        if not ok: continue
        d3=tuple(d1[i]+d2[i] for i in range(n))
        s.add_clause([-gid[d1],-gid[d2],-gid[d3]]); ntri+=1
print(f"triangle clauses: {ntri}",flush=True)

# enumerate height-(n-1) binary subgrids of [t]^(n-1) as tuples of leaves (lex order)
def subgrids(h):
    if h==0:
        yield ((),)
        return
    for b0,b1 in itertools.combinations(range(t),2):
        for X0 in subgrids(h-1):
            for X1 in subgrids(h-1):
                yield tuple((b0,)+x for x in X0)+tuple((b1,)+x for x in X1)
rows=list(subgrids(n-1))
print(f"row subgrids (height {n-1}): {len(rows)}",flush=True)
m=n-1
pts=list(itertools.product(range(t),repeat=m)); pid={p:i for i,p in enumerate(pts)}
P=len(pts)
rowmask=np.array([sum(1<<pid[p] for p in X) for X in rows],dtype=object)
# in-row gap ids for each row subgrid
def ingaps(X):
    out=[]
    for i in range(len(X)):
        for j in range(i+1,len(X)):
            d=(0,)+tuple(X[j][k]-X[i][k] for k in range(m))
            out.append(gid[d])
    return out
row_in=[ingaps(X) for X in rows]
# cross gap id table: cg[h][px][py] = gid of (h, py-px)
cg=np.zeros((t,P,P),dtype=np.int64)
for h in range(1,t):
    for px,p in enumerate(pts):
        for py,q in enumerate(pts):
            cg[h,px,py]=gid[(h,)+tuple(q[k]-p[k] for k in range(m))]
rowpts=[[pid[p] for p in X] for X in rows]
rowbits=np.array([sum(1<<pid[p] for p in X) for X in rows],dtype=np.uint64) if P<=64 else None
assert rowbits is not None
it=0; nclauses=0
while True:
    if time.time()-t0>TL:
        print(f"TIMEOUT after {time.time()-t0:.0f}s, {it} iterations, {nclauses} CEGAR clauses",flush=True); break
    ok=s.solve()
    if not ok:
        print(f"UNSAT: Q_inv({n},{t}) is FALSE  ({it} iterations, {nclauses} CEGAR clauses, {time.time()-t0:.0f}s)",flush=True); break
    model=s.get_model()
    inE=np.zeros(len(gaps)+1,dtype=bool)
    for v in model:
        if v>0 and v<=len(gaps): inE[v]=True
    # independent row subgrids
    indep=[k for k in range(len(rows)) if not any(inE[g] for g in row_in[k])]
    found=[]
    if indep:
        ib=rowbits[indep]
        for h in range(1,t):
            Eh=inE[cg[h]]           # P x P bool: (h, y-x) in E
            compat=np.zeros(P,dtype=np.uint64)
            for px in range(P):
                bits=0
                row=Eh[px]
                for py in range(P):
                    if not row[py]: bits|=(1<<py)
                compat[px]=np.uint64(bits)
            for a,k in enumerate(indep):
                mx=np.uint64((1<<P)-1)
                for px in rowpts[k]: mx&=compat[px]
                hits=np.nonzero((ib & ~mx)==0)[0]
                if len(hits):
                    for b in hits[: max(1,BATCH//8)]:
                        found.append((h,k,indep[b]))
                    if len(found)>=BATCH: break
            if len(found)>=BATCH: break
    if not found:
        print(f"SAT: Q_inv({n},{t}) TRUE; witness |E|={int(inE.sum())} ({it} iterations, {time.time()-t0:.0f}s)",flush=True)
        Elist=[gaps[i-1] for i in range(1,len(gaps)+1) if inE[i]]
        with open(f"finv_witness_n{n}_t{t}.txt","w") as f:
            for d in Elist: f.write(" ".join(map(str,d))+"\n")
        break
    seen=set()
    def norm(d):
        return d if lexpos(d) else tuple(-x for x in d)
    refl=[S for r in range(1,n+1) for S in itertools.combinations(range(n),r)]
    for (h,k1,k2) in found:
        cl=set(row_in[k1])|set(row_in[k2])
        for px in rowpts[k1]:
            for py in rowpts[k2]:
                cl.add(int(cg[h,px,py]))
        key=frozenset(cl)
        if key in seen: continue
        seen.add(key)
        s.add_clause(sorted(cl)); nclauses+=1
        vecs=[gaps[c-1] for c in cl]
        for S in refl:
            img=frozenset(gid[norm(tuple(-v[i] if i in S else v[i] for i in range(n)))] for v in vecs)
            if img in seen: continue
            seen.add(img); s.add_clause(sorted(img)); nclauses+=1
    it+=1
    if it%10==0:
        print(f"  iter {it}: clauses {nclauses}, |E|={int(inE.sum())}, indep rows {len(indep)}, {time.time()-t0:.0f}s",flush=True)
