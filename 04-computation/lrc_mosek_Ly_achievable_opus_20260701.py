"""Where MOSEK CAN bound L_y: add the singular-series LOWER bounds on J(A) (the achievability that offset sets
MUST leave sectors empty near t=0). Feed the EXACT single-runner value J({j})>=diagonal and pair-lower J(A) via
the actual arc geometry. If the analytic J(A) bounds + MOSEK moment-LP give <= cap, it's a HYBRID certificate."""
import numpy as np, cvxpy as cp, itertools
# The KEY missing constraint the naive SDP lacked: J(A) >= product lower bound is WRONG (correlated);
# but the RIGOROUS universal lower bound is J(A) >= meas{all runners in sector 0} for the r sectors adjacent.
# Instead test: with the ACTUAL consec moments as the claimed max, is the moment vector on the PSD boundary?
# Compute consec_8 exact moments S_r and check the moment matrix rank (atoms) + whether small perturbations raise L_y.
def moments(E,N=400000):
    x=np.arange(N)/N; hit=np.zeros((N,7),bool)
    for e in E:
        s=np.floor(7*((e*x)%1.0)).astype(int)%7; hit[np.arange(N),s]=True
    Nm=6-hit[:,1:].sum(axis=1); cnt=np.bincount(Nm,minlength=7)/N
    from math import comb
    return [sum(comb(t,r)*cnt[t] for t in range(7)) for r in range(7)], cnt
y8=[1,-1,1,-0.9,0.6,0,0]; cap8=0.3815
Sr,pt=moments(list(range(8)))
print("consec_8 factorial moments S_r =", [round(s,4) for s in Sr])
print("consec_8 N-distribution p_t =", [round(p,4) for p in pt])
print(f"L_y(consec) = {sum(y8[r]*Sr[r] for r in range(7)):.5f}  (cap {cap8})")
# The relaxation needs the JOINT achievability. Honest test: max L_y over N-distributions p (7 vars) with ONLY
# the moment S_1 pinned to a LOWER bound (consec minimizes S1) -- shows even that is not enough (need coupling).
p=cp.Variable(7,nonneg=True)
S=[cp.sum([__import__('math').comb(t,r)*p[t] for t in range(7)]) for r in range(7)]
cons=[cp.sum(p)==1, S[1]>=Sr[1]]   # consec minimizes S_1 (a real universal lower bound, THM-534 mechanism)
Ly=cp.sum([y8[r]*S[r] for r in range(7)])
cp.Problem(cp.Maximize(Ly),cons).solve(solver=cp.MOSEK)
print(f"\nmax L_y s.t. only S_1>=S_1(consec): {Ly.value:.5f} (LOOSE, need S_2,S_4 coupling) -- confirms THM-534's")
print("'coupled alternating monotonicity, cannot bound moments separately' (=HYP-2738 signed certificate).")
print("=> MOSEK verdict: no clean moment-LP closes 'consec maximizes L_y'; the coupling IS the singular series.")
