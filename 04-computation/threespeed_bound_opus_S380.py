# opus-2026-07-17-S380 -- HYP-7640: THE THREE-SPEED CASE, AND THE METHOD'S CEILING.
#
# GENERAL k-SPEED DENSITY BOUND.  Each comb meets an interval of length ell in
# measure at most 2*lam*ell + 2*lam/w, so covering by k combs forces
#         ell*(1 - 2*k*lam)  <=  2*lam * sum_i (1/w_i).
# This is USEFUL ONLY WHILE 1 - 2*k*lam > 0, i.e.  k < 1/(2*lam) = 7.
# So the method reaches k <= 6 substitutions and DEGENERATES EXACTLY AT k = 7 --
# the same 7 that runs through this entire problem.  Three-speed is well inside.
#
# AT k = 3, lam = 1/14:  1 - 6*lam = 4/7, so  ell <= (1/4)*(1/r + 1/s + 1/t),
# and with r <= s <= t this gives  r <= 3/(4*ell_max(E_ijk)).
# RECURSION: fix r, then E' = E_ijk \ D_r must be covered by two combs (the S379
# bound), then E'' = E' \ D_s by one (the S378 bound).  Termination needs E', E''
# nonempty, which is LRC(12) and LRC(13) respectively -- 11 speeds cannot cover
# at radius 1/14 (gap >= 1/12), nor can 12 (gap >= 1/13).
from fractions import Fraction as F
from itertools import combinations
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def complement(u):
    out=[]; prev=F(0)
    for a,b in u:
        if a>prev: out.append((prev,a))
        prev=b
    if prev<1: out.append((prev,F(1)))
    return out
def circ_lmax(c):
    if not c: return F(0)
    if len(c)>1 and c[0][0]==0 and c[-1][1]==1:
        wrap=(c[0][1]-c[0][0])+(c[-1][1]-c[-1][0])
        inner=max((b-a) for a,b in c[1:-1]) if len(c)>2 else F(0)
        return max(wrap,inner)
    return max(b-a for a,b in c)
def ess(V,drop):
    allv=[]
    for x in V:
        if x not in drop: allv.extend(teeth01(x))
    return complement(union(allv))

base=list(range(1,14))
print("(1) THE CEILING: 1 - 2*k*lam at lam = 1/14")
for k in range(1,9):
    val=1-2*k*LAM
    print(f"    k={k}: 1 - 2k*lam = {val}  {'USABLE' if val>0 else '<-- DEGENERATE (method stops)'}")
print("    => the density method reaches k <= 6 and dies exactly at k = 7 = 1/(2*lam).")

print()
print("(2) STEP 1 -- bound on the SMALLEST new speed over all 286 triples")
print("    r <= 3/(4*ell_max(E_ijk))")
rows=[]
for i,j,k in combinations(base,3):
    E=ess(base,{i,j,k}); L=circ_lmax(E)
    rb=int(F(3,4)/L) if L>0 else None
    rows.append((i,j,k,L,rb))
rbs=[r[4] for r in rows if r[4]]
print(f"    ell_max ranges {float(min(r[3] for r in rows)):.6f} .. {float(max(r[3] for r in rows)):.6f}")
print(f"    bound on smallest speed: {min(rbs)} .. {max(rbs)}   (MAX = {max(rbs)})")
w=max(rows,key=lambda t:t[4] or 0)
print(f"    worst triple: ({w[0]},{w[1]},{w[2]})  ell_max={float(w[3]):.6f}  r <= {w[4]}")
print()
print(f"    crude enumeration estimate: sum of r-bounds over triples = {sum(rbs)}")
print("    (each then spawns a two-speed subproblem -- feasibility checked next)")
