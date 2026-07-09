#!/usr/bin/env python3
"""
klein-2026-07-08-S194: the exact 𝒲̂ formula gives BOTH resonance sums a-priori.

Using 𝒲̂(n) = (-1)^r (6/7)^{k-1-r} prod b0(n_i) (1[sigma=0]-c(sigma)):
 (A) DENSITY-FLOOR decorrelation:  E[W] - (6/7)^k = Sum_{n!=0, n.e=0} 𝒲̂(n).
 (B) THM-664 GRID resonance:       E_grid[W] - (6/7)^k = Sum_{n!=0, Vmax|n.e} 𝒲̂(n).
Both are the SAME 𝒲̂ summed over a resonance lattice. Verify each truncated sum
(|n|<=H) CONVERGES to the directly-computed LHS, and exhibit the geometric
per-coordinate decay factor (7/6)/pi ~ 0.37 that makes them a-priori.
"""
import numpy as np, cmath, itertools
INV7=1.0/7.0
def W_uncovered(ph):
    p=np.sort(np.mod(np.asarray(ph,float),1.0)); g=np.diff(p); g=np.append(g,1.0-p[-1]+p[0])
    return np.maximum(g-INV7,0.0).sum()
def b0(m): return 1.0/7.0 if m==0 else (cmath.exp(2j*np.pi*m/7.0)-1.0)/(2j*np.pi*m)
def c_of(s): return 1.0/7.0 if s==0 else (1.0-cmath.exp(-2j*np.pi*s/7.0))/(2j*np.pi*s)
def What(n,k):
    nz=[int(v) for v in n if v!=0]; r=len(nz); sig=int(sum(n))
    pr=1.0+0j
    for v in nz: pr*=b0(v)
    return ((-1)**r)*(6.0/7.0)**((k-1)-r)*pr*((1.0 if sig==0 else 0.0)-c_of(sig))

def E_W_direct(e,N=400000):
    e=np.asarray(e,float); x=(np.arange(N)+0.5)/N
    Ph=np.mod(np.outer(x,e),1.0); Ph.sort(axis=1)
    g=np.diff(Ph,axis=1); g=np.concatenate([g,(1.0-Ph[:,-1]+Ph[:,0])[:,None]],axis=1)
    return np.maximum(g.max(axis=1)*0+ (g-INV7),0).sum(axis=1).mean() if False else \
           np.maximum(g-INV7,0.0).sum(axis=1).mean()
def E_grid_W(e,V):
    e=np.asarray(e,float); j=np.arange(V); x=j/V
    Ph=np.mod(np.outer(x,e),1.0); Ph.sort(axis=1)
    g=np.diff(Ph,axis=1); g=np.concatenate([g,(1.0-Ph[:,-1]+Ph[:,0])[:,None]],axis=1)
    return np.maximum(g-INV7,0.0).sum(axis=1).mean()

# cluster: e=(0,e1,...,e_{k-1}); use small k for full-lattice enumeration
for (e, Vmax, tag) in [((0,1,3), 7, "k=3 tiny"),
                       ((0,1,3,7), 13, "k=4"),
                       ((0,2,5,11), 17, "k=4 spread")]:
    k=len(e); ep=np.array(e[1:])   # e_1..e_{k-1} (the phase coords)
    print(f"\n===== {tag}: e={e}, Vmax={Vmax}, (6/7)^k={(6/7)**k:.6f} =====")
    EW=E_W_direct(e); Eg=E_grid_W(e,Vmax)
    print(f"  E[W] (direct)      = {EW:.6f}   => decorrelation E[W]-(6/7)^k = {EW-(6/7)**k:+.6f}")
    print(f"  E_grid[W] (direct) = {Eg:.6f}   => grid resid E_grid-(6/7)^k  = {Eg-(6/7)**k:+.6f}")
    for H in (4,8,14,22):
        sA=0.0+0j; sB=0.0+0j
        for n in itertools.product(range(-H,H+1),repeat=k-1):
            if all(v==0 for v in n): continue
            ne=int(np.dot(n,ep))
            wv=What(n,k)
            if ne==0: sA+=wv
            if ne % Vmax == 0: sB+=wv
        print(f"    H={H:2d}: (A) sum_{{n.e=0}} 𝒲̂ = {sA.real:+.6f}   (B) sum_{{Vmax|n.e}} 𝒲̂ = {sB.real:+.6f}")
    print(f"  target (A)={EW-(6/7)**k:+.6f}   target (B)={Eg-(6/7)**k:+.6f}")

# decay demonstration: per-coordinate suppression factor
print("\n--- a-priori decay: adding a nonzero coord multiplies |𝒲̂| bound by <= (7/6)/pi ---")
print(f"    (7/6)/pi = {(7/6)/np.pi:.4f} < 1  (geometric in # nonzero coords)")
print("    |𝒲̂(n)| <= (6/7)^{k-1-r} prod_{n_i!=0} 1/(pi|n_i|) * min(6/7, 1/(pi|sigma|)); =0 if 7|n_i or 7|sigma")
