"""
(A) SECOND-VARIATION: is the AP a strict MIN of M over REAL speed perturbations, or is the integer
    structure essential (AP a real saddle)?  Perturb {1..13} by real eps in each coordinate.
(B) CUSP: divisor-loaded {1..11,13,84m} as m->inf -- where does M, L go? (one speed -> cusp)
"""
import numpy as np
def M_grid(S,Q=200000,n=14):
    t=np.arange(1,Q)/Q
    f=np.ones(Q-1)
    for v in S:
        r=np.abs(((v*t+0.5)%1)-0.5)   # ||v t||
        f=np.minimum(f,r)
    return f.max()
def L_grid(S,Q=200000,n=14):
    t=np.arange(Q)/Q; safe=np.ones(Q,bool)
    for v in S:
        r=np.abs(((v*t+0.5)%1)-0.5)
        safe &= (r> 1/n)
    return safe.mean()
ap=list(range(1,14))
print(f"AP M={M_grid(ap):.6f} (=1/14={1/14:.6f}), L={L_grid(ap):.6f}")
print("\n(A) real perturbations of the AP -- does M dip BELOW 1/14?  (eps on one coordinate)")
print(f"{'perturb':>22} {'M':>10} {'M-1/14':>10}")
import itertools
for i,eps in [(13,0.1),(13,-0.1),(13,0.5),(1,0.1),(1,-0.1),(7,0.3),(7,-0.3),(13,1.0)]:
    S=ap.copy(); S[i-1]=i+eps
    M=M_grid(S); print(f"  v{i}->{i+eps:<6}        {M:>10.6f} {M-1/14:>+10.6f}")
# random small real perturbations
np.random.seed(0); below=0; mn=1.0
for _ in range(40):
    S=[v+np.random.uniform(-0.2,0.2) for v in ap]
    M=M_grid(S,80000); mn=min(mn,M); below+= (M<1/14-1e-6)
print(f"  40 random real perturbations (|eps|<0.2): {below} gave M<1/14; min M = {mn:.6f}")
print("\n(B) CUSP: {1..11,13,84m} as m grows (one speed -> infinity):")
print(f"{'m':>4} {'last':>7} {'M':>10} {'L':>10}")
core=[1,2,3,4,5,6,7,8,9,10,11,13]
for m in [1,2,3,5,10,20,50]:
    S=core+[84*m]; M=M_grid(S, max(200000, 84*m*40)); L=L_grid(S, max(200000,84*m*40))
    print(f"{m:>4} {84*m:>7} {M:>10.6f} {L:>10.6f}")
Mc=M_grid(core,200000); Lc=L_grid(core,200000)
print(f"  core {{1..11,13}} alone: M={Mc:.6f} (=1/12={1/12:.6f}), L={Lc:.6f}")
print("  cusp limit: huge speed decouples (thin danger) => M->1/12, L-> (6/7)*L(core)?")
print(f"  (6/7)*L(core) = {6/7*Lc:.6f}")
