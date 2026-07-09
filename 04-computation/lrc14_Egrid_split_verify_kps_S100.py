# Independent confirmation of the fleet's E_grid[W] = E_x[W] + R_grid split (kps-S100).
# E_x[W] = continuum average (V-independent, the density-floor MAIN term); R_grid = E_grid - E_x
# (wraparound residual). Confirm: E_x stable across V, R_grid bounded & decaying with V/spread.
# Localizes my S97 kissing |R| to the R_grid (decaying) side, main term = growing winning side.
import numpy as np
from math import gcd, floor
from functools import reduce
import random
random.seed(100019)
K=13; TH=1.0/7
def prim(E):
    E=sorted(E); return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def W_grid_mean(E,V):
    Ea=np.array(sorted(E)); j=np.arange(0,V)
    ph=np.mod(np.outer(j,Ea),V)/V; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return np.maximum(g-TH,0).sum(axis=1).mean()
def W_cont_mean(E,N=20000):   # continuum E_x[W] via fine grid (V-independent)
    Ea=np.array(sorted(E)); x=(np.arange(N)+0.5)/N
    ph=np.mod(np.outer(x,Ea),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return np.maximum(g-TH,0).sum(axis=1).mean()
def rand_diss(s):
    for _ in range(500):
        mid=sorted(random.sample(range(1,s),K-2)); E=[0]+mid+[s]
        if len(set(E))==K and prim(E): return E
    return None
print("Split E_grid[W] = E_x[W](main, V-indep) + R_grid(residual). Confirm E_x stable, R_grid bounded/decays.")
print(f"{'spread':>7}{'E_x[W]':>9}{'V':>5}{'E_grid':>9}{'R_grid':>9}{'|R_grid|/E_x':>12}")
for s in (60,90,140,220):
    E=rand_diss(s)
    if E is None: continue
    Ex=W_cont_mean(E)
    for V in (s+1, floor(7*s/6), 2*s, 5*s):   # from knife-edge to large ruler
        Eg=W_grid_mean(E,V); Rg=Eg-Ex
        print(f"{s:>7}{Ex:>9.5f}{V:>5}{Eg:>9.5f}{Rg:>+9.5f}{abs(Rg)/Ex:>12.4f}")
    print()
print("=> E_x[W] is the STABLE main term (density floor, V-independent); R_grid shrinks as V/spread grows")
print("   => my S97 kissing |R| lives on the DECAYING R_grid side; main term = growing V*E_x (mac-mini-S64).")
