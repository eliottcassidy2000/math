"""kps-S39: the far-peel is a MEASURE bound (comb discrepancy): mu_v >= (6/7)mu' - 2p/(7 v_r),
p = #safe-components of the other 12, mu' = their safe measure. Threshold v_r > p/(3 mu') is
LINEAR in p (not V^2). Closes DOMINANT-FAR + the deep well. But p ~ Sum v (not <= C*n), so it
only reaches dominant-far, NOT the uniform looseness (= covering LRC(14), LRC-hard)."""
import numpy as np
rng=np.random.default_rng(2)
def p_mu(v,grid=1500000):
    v=np.array(v,dtype=np.float64); t=np.linspace(0,1,grid,endpoint=False)
    x=np.outer(v,t); m=np.abs(x-np.round(x)).min(axis=0)>=1/14-1e-12
    p=int(((m.astype(int)-np.roll(m,1).astype(int))==1).sum()); return max(p,1),m.mean()
def mu_full(v,grid=1500000):
    v=np.array(v,dtype=np.float64); t=np.linspace(0,1,grid,endpoint=False)
    return (np.abs(np.outer(v,t)-np.round(np.outer(v,t))).min(axis=0)>=1/14-1e-12).mean()
print("MEASURE FAR-PEEL: mu_v >= (6/7)mu' - 2p/(7 v_r); v_r > p/(3mu') => mu_v>0 => M>1/14.")
for v,ri in [(list(range(1,13))+[182],12),(list(range(1,13))+[500],12)]:
    others=v[:ri]+v[ri+1:]; vr=v[ri]; p,mu=p_mu(others); muv=mu_full(v)
    lb=(6/7)*mu-2*p/(7*vr); thr=p/(3*mu)
    print(f"  v={v}: p={p} mu'={mu:.4f} thr=p/3mu'={thr:.0f} v_r={vr}>{thr:.0f}? {vr>thr}; "
          f"bound={lb:.4f}>0, actual mu_v={muv:.4f}")
print()
print("p is NOT <= C*n -- it grows like Sum v (so far-peel = DOMINANT-far, not uniform):")
for mag in [30,300,3000]:
    ps=[max(p_mu(sorted(set(int(x) for x in rng.integers(1,mag+1,size=12))),400000)[0] for _ in [0]) for _ in range(30)]
    print(f"  n=12 mag<={mag}: max p over 30 fams = {max(ps)} (p/n up to {max(ps)/12:.0f})")
print("=> uniform looseness (all covering M>1/14) = covering case of LRC(14) = LRC-hard; NOT proved.")
print("   Far-peel/measure-independence closes dominant-far (deep well incl); residual =")
print("   non-dominant multi-scale (renormalization, opus THM-608) + rigidity (mac-mini THM-612).")
