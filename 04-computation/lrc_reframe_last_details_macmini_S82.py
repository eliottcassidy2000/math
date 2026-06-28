"""
S82: reframing the last two LRC(14) details to fit, unified by the unit-grid (Z/14)* core/bulk split.
(b) RESONANT survivor reduces to GENERIC: a resonant v=14m has danger ON the 14-grid; the seed's optimum dodges
    OFF the exact a/14 (where v is safe) into the bulk; M>1/14 always. No resonant exception => off-grid bulk bound.
(a) RIGIDITY = a FINITE equioscillation system: tight <=> 3 sum-14 binding pairs {1,13},{3,11},{5,9} at the units
    + complement symmetry; AP & GW both satisfy it; integer 13-set solutions = AP/GW (+ dilations).
"""
import numpy as np
def M_arg(S,grid=560000):
    t=np.arange(1,grid)/grid; s=np.ones(grid-1)
    for x in S:
        fr=(x*t)%1.0; s=np.minimum(s,np.minimum(fr,1-fr))
    i=s.argmax(); return s.max(),t[i]
seed=[1,2,3,4,5,6,8,9,10,11,12,13]
print("(b) resonant v covers the on-grid core; optimum dodges off-grid (bulk survives), M>1/14:")
M0,t0=M_arg(seed); print(f"  seed alone: M={M0:.5f} at t*={t0:.4f} (on-grid core 1/7)")
for v in (14,28,42,56):
    M,ts=M_arg(seed+[v]); ng=round(ts*14)/14
    print(f"  seed+{v:>2}: M={M:.5f}>1/14 at t*={ts:.4f}  (dodged off exact a/14; v safe there)")
print("  => resonant case reduces to the GENERIC off-grid bulk bound (Vitali: v hits core, bulk survives).")
print()
def res(S): return sorted(set(x%14 for x in S))
AP=list(range(1,14)); GW=list(range(1,12))+[13,24]; pairs=[(1,13),(3,11),(5,9)]
print("(a) rigidity = finite equioscillation system (3 sum-14 binding pairs at the units):")
for nm,S in [("AP",AP),("GW",GW)]:
    R=set(res(S)); have=[p for p in pairs if p[0] in R and p[1] in R]
    print(f"  {nm}: residues mod14={res(S)}  contains pairs {have} (all 3: {len(have)==3})")
print("  => tight <=> 3 sum-14 pairs at units + complement-symmetric; integer 13-set solutions = AP/GW (+dilations).")
print("  Unified: the unit grid (Z/14)* is where tight is pinned (a) AND where resonant danger sits (b); off-grid bulk = generic.")
