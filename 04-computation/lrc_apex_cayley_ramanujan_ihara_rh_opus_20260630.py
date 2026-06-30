"""
MEASURE column = SPECTRAL: apex gap g(O)=min_k|sum omega^{kx}|^2 = a Cayley spectral gap; THM-590 0.198 =
2+lambda_min(C_7). Apex Cayley graphs (C_p, Paley_p, K_p) are RAMANUJAN => their Ihara zeta satisfies RH
(non-trivial poles on |u|=1/sqrt(k-1)). EXISTENCE/M column = covering radius (even-n 1/n, equally-spaced AP).
"""
import cmath, math
import numpy as np
p=7; w=cmath.exp(2j*cmath.pi/p)
def apex_gap(O): return min(abs(sum(w**((k*x)%p) for x in O))**2 for k in range(1,p))
# apex Cayley graphs on Z_7 and Ramanujan / Ihara-RH
def cayley_adj(conn,p):
    A=np.zeros((p,p))
    for i in range(p):
        for c in conn: A[i][(i+c)%p]=1
    return A
def ihara_poles_on_circle(A):
    p=A.shape[0]; k=int(A.sum(axis=1)[0]); 
    eig=np.linalg.eigvalsh(A)
    rh_ok=True; poles=[]
    for lam in eig:
        if abs(lam-k)<1e-9: continue   # trivial
        disc=lam*lam-4*(k-1)
        if disc<=1e-9:  # complex pole => |u|=1/sqrt(k-1)
            u=abs(complex(lam,math.sqrt(-disc)))/(2*(k-1))
            poles.append(round(u,5))
        else:
            rh_ok=False  # real pole off the circle (non-Ramanujan)
    return k,sorted(set(round(float(e),4) for e in eig)),rh_ok,1/math.sqrt(k-1) if k>1 else None,poles
print("APEX CAYLEY GRAPHS on Z_7 (the MEASURE column) -- Ramanujan + Ihara-RH:")
for name,conn in [("C_7 (doublet +-1)",[1,6]),("Paley_7 (QR {1,2,4})",[1,2,4]),("K_7 (full)",[1,2,3,4,5,6])]:
    A=cayley_adj(conn,p)
    k,eig,rh,bound,poles=ihara_poles_on_circle(A)
    print(f"  {name:>22}: {k}-regular, eigs={eig}, Ramanujan(|nontriv|<=2sqrt(k-1)={2*math.sqrt(k-1):.3f}): {max(abs(e) for e in eig if abs(e-k)>1e-6)<=2*math.sqrt(k-1)+1e-6}")
    print(f"        Ihara-RH (poles on |u|=1/sqrt(k-1)={bound:.4f}): {rh}; poles={poles}")
print()
print("THE APEX GAP as 2+lambda_min (THM-590, mac-mini HYP-3594):")
for name,O in [("doublet {0,1}",[0,1]),("{1,3,5}",[1,3,5]),("Z_7 (cusp)",[0,1,2,3,4,5,6])]:
    g=apex_gap(O)
    print(f"  {name:>16}: apex gap = {g:.5f}")
lamC7=min(2*math.cos(2*math.pi*k/7) for k in range(1,7))
print(f"  C_7 cycle lambda_min = {lamC7:.5f}; 2+lambda_min = {2+lamC7:.5f} = apex doublet gap 4cos^2(3pi/7)=0.198")
print(f"  => the apex gap is the C_7 SLACK from the Ramanujan/spectral bound -2: (-2)-lambda_min would be 0 at tight;")
print(f"     here 2+lambda_min=0.198 = how far C_7's min eig sits ABOVE -2 (genus-1 spectral signature).")
print()
print("EXISTENCE/M column = COVERING RADIUS (the corrected covering-min):")
print("  EVEN n: even block 2*{1..n-1} at t=1/(2n) = equally-spaced AP = the n-th roots = C_n vertices; radius 1/n.")
print("  => the covering-min IS the C_n cycle covering radius (equally-spaced optimum). ODD n: C_n unreachable by")
print("     an all-even covering set (needs odd mult of n) => covering radius > 1/n (realizability).")
print("  TWO COLUMNS, both the cycle C_*: MEASURE = C_p SPECTRAL slack (Ramanujan/Ihara-RH); M = C_n COVERING radius.")
