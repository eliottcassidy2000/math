#!/usr/bin/env python3
"""mac-mini-S80: pursue the DELSARTE dual-existence LP concretely.
Certificate: a degree-D positive density nu (trig poly, nu>=0, INT nu=1) with INT W dnu < 1
=> W=0 on a positive nu-set => safe point => M(S) >= level. W = danger count at 'level'.
LP: min over c_h (h=1..D) of INT W*nu, s.t. nu(t_j)>=0 on a grid. Compute the Delsarte bound
vs 1 at increasing D for the deep well (knife-edge at its own M) and other covering families.
TEST: (a) at level 1/14, do all covering families certify (<1) at finite D? (b) degree scaling
~ 1/(M-level)? (c) at level 14/183, is the deep well exactly at the knife-edge (bound ->1)?"""
import numpy as np
from scipy.optimize import linprog
def W_grid(S,level,N):
    ts=(np.arange(N)+0.5)/N
    W=np.zeros(N)
    for v in S:
        r=(v*ts)%1.0; d=np.minimum(r,1-r)
        W+=(d<level).astype(float)
    return ts,W
def delsarte_bound(S,level,D,N=1200):
    ts,W=W_grid(S,level,N)
    # variables: a_h=Re(c_h), b_h=Im(c_h), h=1..D (2D vars). nu(t)=1+2*sum(a_h cos - b_h sin).
    C=np.cos(2*np.pi*np.outer(ts,np.arange(1,D+1)))  # N x D
    Sn=np.sin(2*np.pi*np.outer(ts,np.arange(1,D+1)))
    # nu_j = 1 + 2*(C@a - Sn@b); constraint nu_j>=0 => -2C a +2Sn b <= 1
    A_ub=np.hstack([-2*C, 2*Sn]); b_ub=np.ones(N)
    # objective INT W nu = mean(W*nu) = mean(W) + 2*mean(W*C)@a - 2*mean(W*Sn)@b
    obj=np.hstack([2*(W[:,None]*C).mean(0), -2*(W[:,None]*Sn).mean(0)])
    const=W.mean()
    res=linprog(obj, A_ub=A_ub, b_ub=b_ub, bounds=[(None,None)]*(2*D), method='highs')
    if not res.success: return None
    return const+res.fun
def is_cov(S,n=14): return all(any(v%q==0 for v in S) for q in range(2,n+1))

deep=[*range(1,13),182]; tight=[*range(1,12),13,84]; loose=[2,3,5,7,11,13,17,19,23,29,31,37,42]
print("DELSARTE bound = min_{deg<=D positive nu} INT W dnu. <1 => degree-D existence certificate.\n")
print("(a) at level 1/14 (LRC threshold): all covering M>=14/183>1/14 => loose, should certify:")
for nm,S in [("deep well {1..12,182}",deep),("{1..11,13,84}",tight),("{2..14}",list(range(2,15)))]:
    S=sorted(set(S)); row=[]
    for D in [20,40,80,160]:
        b=delsarte_bound(S,1/14,D)
        row.append(f"D={D}:{b:.3f}" if b is not None else f"D={D}:fail")
    cert=next((D for D in [20,40,80,160] if (delsarte_bound(S,1/14,D) or 9)<1),None)
    print(f"   {nm:22s} cov={is_cov(S)}: "+"  ".join(row)+f"  (certifies<1 by D~{cert})")
print("\n(b) at level 14/183 (covering-MIN): deep well M=14/183 EXACTLY => knife-edge, bound->1;")
print("    other covering (M>14/183) should still certify:")
for nm,S in [("deep well (M=14/183)",deep),("{1..11,13,84}(M=7/89)",tight)]:
    S=sorted(set(S)); row=[]
    for D in [40,80,160,320]:
        b=delsarte_bound(S,14/183,D)
        row.append(f"D={D}:{b:.4f}" if b is not None else "fail")
    print(f"   {nm:22s}: "+"  ".join(row))
print("\nREAD: if deep-well bound at 14/183 stays >=~1 at all D (knife-edge) but <1 at level 1/14")
print("at finite D, and cert degree ~1/(M-level), the Delsarte LP is EQUIVALENT to the finite-Vmax")
print("glue (degree=Vmax); certifies loose constructively, deep well needs rigidity. Honest test.")
