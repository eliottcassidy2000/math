"""SYNTHESIS: L_y(E) two regimes. STRUCTURED (consec) = the MAX (THM-534). SPREAD (large offsets) -> the
OCCUPANCY/equidistribution limit (i.i.d. runners) = FAR BELOW consec => this IS mac-mini's far-element
equidistribution (HYP-3786). So 'consec maximizes L_y' splits: bounded (finite check, THM-534) + far (equidist)."""
import numpy as np, random
from math import comb
def Ly(E,y,N=300000):
    x=np.arange(N)/N; hit=np.zeros((N,7),bool)
    for e in E:
        s=np.floor(7*((e*x)%1.0)).astype(int)%7; hit[np.arange(N),s]=True
    Nm=6-hit[:,1:].sum(axis=1)
    Sr=[np.mean([comb(int(t),r) for t in Nm]) for r in range(len(y))]  # slow; use bincount
    cnt=np.bincount(Nm,minlength=7)/N
    return sum(y[r]*sum(comb(t,r)*cnt[t] for t in range(7)) for r in range(len(y))), cnt
# k=8 magic function: L_y = 1 - S1 + S2 - 0.9 S3 + 0.6 S4
y8=[1,-1,1,-0.9,0.6,0,0]; cap8=0.3815
def occ_Ly(k,y):  # k-1 i.i.d. runners over 7 sectors: S_r = C(6,r)*((7-r)/7)^(k-1)
    Sr=[comb(6,r)*((7-r)/7)**(k-1) for r in range(7)]; return sum(y[r]*Sr[r] for r in range(7))
consec=list(range(8)); Lc,_=Ly(consec,y8)
Locc=occ_Ly(8,y8)
print(f"k=8: L_y(consec)={Lc:.4f} (THM-534 MAX)   cap_8={cap8:.4f}   occupancy/equidist limit={Locc:.4f}")
print(f"\nL_y for increasingly SPREAD offset sets (0 + geometric-ish large offsets) -> occupancy limit:")
rng=random.Random(1)
for scale in [1,3,10,40,150]:
    Ls=[]
    for _ in range(8):
        E=[0]+sorted(rng.sample(range(1,7*scale+2),7)); L,_=Ly(E,y8,120000); Ls.append(L)
    print(f"   offsets ~ up to {7*scale}: mean L_y = {np.mean(Ls):.4f}  (consec {Lc:.4f}, occ {Locc:.4f})")
print(f"\n=> STRUCTURED consec = MAX (0.358 < cap 0.382); SPREAD -> occupancy {Locc:.3f} << consec. Two regimes:")
print(f"   bounded (THM-534 finite check) + far/spread (mac-mini HYP-3786 equidistribution on the lonely set).")
print(f"   The 'consec maximizes L_y' crux = these two regimes; the boolean-moment SDP is LOOSE without the")
print(f"   offset-set achievability = the singular series = where the far-element equidistribution bites.")
