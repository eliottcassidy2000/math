"""Cayley transform (tournament -> SO(7) -> 7 points on the circle -> Verblunsky), Paley/QR comparison, tiling model."""
import numpy as np
from itertools import combinations
V=list(range(7))
def rot_tour(conn):
    A=np.zeros((7,7),int)
    for i in V:
        for j in V:
            if i!=j and ((j-i)%7) in conn: A[i,j]=1
    return A
def verblunsky(moments):
    N=len(moments)-1; r=np.array(moments,dtype=complex); a=np.zeros(N+1,dtype=complex); a[0]=1; E=r[0].real; al=[]
    for n in range(1,N+1):
        acc=r[n]+sum(a[j]*r[n-j] for j in range(1,n)); k=-acc/E if abs(E)>1e-14 else 0
        anew=a.copy()
        for j in range(1,n): anew[j]=a[j]+k*np.conj(a[n-j])
        anew[n]=k; a=anew; E=E*(1-abs(k)**2); al.append(-np.conj(k))
    return al
for name,conn in [("ROTATIONAL {1,2,3}",{1,2,3}),("PALEY/QR {1,2,4}",{1,2,4})]:
    A=rot_tour(conn); S=(A-A.T).astype(float)
    ev=np.linalg.eigvals(S); 
    U=np.linalg.solve(np.eye(7)+S, np.eye(7)-S)  # Cayley (I-S)(I+S)^-1 in SO(7)
    uev=np.linalg.eigvals(U); ang=np.sort(np.angle(uev))
    print(f"\n{name}: |skew eigenvalue|_max={max(abs(ev)):.4f} (sqrt7={np.sqrt(7):.4f})")
    print(f"  Cayley U eigenvalues on circle, angles/2pi = {np.round(ang/(2*np.pi),4)}")
    # are they odd/14 or k/7 ?
    frac=np.sort((ang/(2*np.pi))%1)
    print(f"  angles mod 1 = {np.round(frac,4)}   (7th roots k/7={np.round([k/7 for k in range(7)],4)})")
    # Verblunsky of uniform measure on the 7 Cayley-eigenvalue points
    mom=[np.mean(uev**k) for k in range(7)]
    al=verblunsky(mom); print(f"  Verblunsky |alpha| of Cayley spectrum = {[round(abs(x),3) for x in al]}")
# ---- TILING model: base path 0->1->...->6, tiles=chords (i,j) j>=i+2; tournament=tiling (flipped = 'down' arcs) ----
print("\n== TILING model (base Ham path 0->1->..->6; 15 tiles = chords; transitive=empty tiling) ==")
for name,conn in [("ROTATIONAL {1,2,3}",{1,2,3}),("PALEY/QR {1,2,4}",{1,2,4})]:
    A=rot_tour(conn); flipped=[(i,j) for i in V for j in V if j>=i+2 and A[i,j]==0]  # tile down = j->i (flipped from transitive i->j)
    byrange={}
    for (i,j) in flipped: byrange.setdefault(j-i,[]).append((i,j))
    print(f"  {name}: {len(flipped)} flipped tiles (Hamming dist to transitive); by chord-length { {r:len(v) for r,v in sorted(byrange.items())} }")
    print(f"      flipped = the LONG chords {sorted(set(j-i for i,j in flipped))} (>= 4 rotational / {{3,5,6}} Paley)")
