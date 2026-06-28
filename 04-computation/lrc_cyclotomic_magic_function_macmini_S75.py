"""
S75: BUILD the cyclotomic Delsarte/Beurling-Selberg magic function for the LRC(14) cover bound.
Target (THM-576): cap_P = meas(lonely(P)), lonely(P)={x: ||p x||>=1/14 for all p in P}; we want a CERTIFIED
lower bound on cap_P via a trig-polynomial MINORANT F(x) <= 1_lonely(x), with F̂_0=∫F the bound.
The MAGIC properties to verify:
  (1) ∫F -> cap_P as degree N grows (the certificate is sharp),
  (2) the equioscillation (active F=1_lonely) points are the DE MOIVRE / cyclotomic angles (sector boundaries j/7),
  (3) f̂_n >= 0 (Bochner/PSD) -- the Cohn-Elkies "magic" positivity.
Build via LP (scipy.linprog), F(x)=a0+sum a_n cos(2pi n x) (even; lonely is x->1-x symmetric).
"""
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction as F_

def lonely_indicator(P, x):
    out=np.ones_like(x)
    for p in P:
        fr=(p*x)%1.0
        out*= ((fr>=1/14-1e-12)&(fr<=13/14+1e-12)).astype(float)
    return out

def meas_lonely(P, M=200004):
    x=(np.arange(M)+0.5)/M
    return lonely_indicator(P,x).mean()

def build_minorant(P, N, M=4004):
    """optimal even trig-poly minorant F<=1_lonely, maximize a0=∫F. returns a (coeffs), ∫F, active pts."""
    x=(np.arange(M)+0.5)/M
    L=lonely_indicator(P,x)
    # F(x_j)=sum_{n=0}^N a_n cos(2 pi n x_j) <= L(x_j);  maximize a0 = ∫F
    C=np.cos(2*np.pi*np.outer(x,np.arange(N+1)))         # M x (N+1)
    # linprog minimizes c^T a ; we maximize a0 => c=[-1,0,...,0]
    c=np.zeros(N+1); c[0]=-1.0
    res=linprog(c, A_ub=C, b_ub=L, bounds=[(None,None)]*(N+1), method="highs")
    a=res.x
    Fvals=C@a
    active=x[np.abs(Fvals-L)<1e-4]                        # equioscillation / touch points
    return a, a[0], active, x, Fvals, L

print("="*92)
print(" BUILDING the cyclotomic magic function: optimal minorant F<=1_lonely(P), ∫F -> cap")
print("="*92)
demoivre_x=sorted(set([round(j/7,4) for j in range(8)]))   # sector boundaries j/7 (the cyclotomic pts in x)
print(f" de Moivre / cyclotomic boundary points in x (j/7): {demoivre_x}")
for P,label in [((1,13),"{1,13} j=2 cap=66/91"),
                ((1,12,13),"{1,12,13} j=3"),
                ((1,11,12,13),"{1,11,12,13} j=4"),
                ((1,5,7,8,9),"{1,5,7,8,9} j=5 k=8 BINDING")]:
    cap=meas_lonely(P)
    print(f"\n P={label}:  meas(lonely)={cap:.5f}")
    for N in (7,14,28,56):
        a,intF,active,x,Fv,L=build_minorant(P,N)
        # f-hat >=0 ? (a_n are the cosine Fourier coeffs; Bochner: a_n>=0)
        fhat_ok=bool(np.all(a[1:]>=-1e-6))
        nneg=int(np.sum(a[1:]<-1e-6))
        # active points mod 1/7 (are they at sector boundaries?)
        act_mod=sorted(set(round((ax*7)%1,2) for ax in active))[:8]
        print(f"   N={N:>3}: ∫F={intF:.5f}  gap={cap-intF:+.5f}  f̂>=0:{fhat_ok}({nneg} neg)  active*7 mod1≈{act_mod}")
print("="*92)
print(" If ∫F->cap, f̂>=0, and active points sit at sector boundaries (j/7, the cyclotomic) => MAGIC FUNCTION.")
