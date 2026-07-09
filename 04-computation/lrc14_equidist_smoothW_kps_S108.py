# Run the equidistribution on the SMOOTH surrogate W and show it beats the SHARP indicator (kps-S108).
# For a fixed cluster E, vary Vmax: measure the convergence rate of
#   SMOOTH:  E_grid[W](V) = (1/V) sum_j W(je/V)  ->  E_x[W] (continuum).  R_grid = E_grid - E_x ~ C/V^p.
#   SHARP:   rho_K(V) = #{good j}/V  ->  rho* = meas(G*).  disc = rho_K - rho* ~ C'/V^q.
# S107 prediction: smooth p ~ 2 (W is C^0 PL, desingularizes the pinches), sharp q ~ 1 (grid-invisible
# pinches = jump discontinuities of the indicator). => W converges an ORDER faster; the equidistribution
# is CLEAN on W. Good period exists <=> E_grid[W] > 0.
import numpy as np
TH=1.0/7
def W_and_good(E,V):
    Ea=np.array(E,float); j=np.arange(0,V)
    ph=np.mod(np.outer(j,Ea),V)/V; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    W=np.maximum(g-TH,0).sum(axis=1)
    good=(g.max(axis=1)>TH)
    return W.mean(), good.mean()   # E_grid[W], rho_K
def cont(E,N=400000):
    Ea=np.array(E,float); x=(np.arange(N)+0.5)/N
    ph=np.mod(np.outer(x,Ea),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    W=np.maximum(g-TH,0).sum(axis=1); good=(g.max(axis=1)>TH)
    return W.mean(), good.mean()   # E_x[W], rho*
def fit_p(Vs, discs):
    m=[(np.log(V),np.log(abs(d))) for V,d in zip(Vs,discs) if abs(d)>1e-9]
    if len(m)<2: return float('nan')
    lx,ly=zip(*m); A=np.vstack([lx,np.ones(len(lx))]).T
    slope,_=np.linalg.lstsq(A,np.array(ly),rcond=None)[0]; return -slope
sets={
 "dissoc k=13 (spread 35)":[0,1,4,9,11,16,20,23,25,28,30,33,35],
 "7-struct k=13 (spread 82)":[0,7,14,21,26,29,37,44,51,58,67,75,82],
}
for name,E in sets.items():
    Exw, rhostar = cont(E)
    print(f"\n{name}: E_x[W]={Exw:.5f} (continuum), rho*={rhostar:.4f}")
    print(f"  {'V':>6}{'E_grid[W]':>11}{'R_grid':>11}{'rho_K':>9}{'rho_K-rho*':>12}{'good?':>6}")
    Vs=[200,400,800,1600,3200,6400]; Rs=[]; Ds=[]
    for V in Vs:
        Egw,rk=W_and_good(E,V); R=Egw-Exw; D=rk-rhostar; Rs.append(R); Ds.append(D)
        print(f"  {V:>6}{Egw:>11.5f}{R:>+11.5f}{rk:>9.4f}{D:>+12.5f}{'YES' if Egw>0 else 'NO':>6}")
    p_smooth=fit_p(Vs,Rs); p_sharp=fit_p(Vs,Ds)
    print(f"  => SMOOTH W discrepancy R_grid ~ 1/V^{p_smooth:.2f}   |   SHARP rho_K discrepancy ~ 1/V^{p_sharp:.2f}")
print("\n=> if smooth p~2 >> sharp q~1: the smooth surrogate W equidistributes an ORDER faster (the")
print("   grid-invisible pinches are jump-discontinuities for the SHARP indicator = O(1/V), but only")
print("   corners for the C^0 W = O(1/V^2). Running the equidistribution on W is the desingularized route.")
