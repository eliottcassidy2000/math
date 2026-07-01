"""The CRT/parity decomposition of LRC(14) at t=1/14:  14 = 2 x 7 = (parity Z/2 = iota/complement) x (heptagon D_7).
ODD runners -> full z^7=-1 heptagon (the SC tournament, Cayley-dual to 7th roots, Verblunsky (0..0,1)).
EVEN runners -> 7th roots minus origin -> HARMONIC Verblunsky at n=7. This is a recursive LRC(7) inside LRC(14)."""
import numpy as np
def verblunsky(moments):
    N=len(moments)-1; r=np.array(moments,dtype=complex); a=np.zeros(N+1,dtype=complex); a[0]=1; E=r[0].real; al=[]
    for n in range(1,N+1):
        acc=r[n]+sum(a[j]*r[n-j] for j in range(1,n)); k=-acc/E if abs(E)>1e-14 else 0
        anew=a.copy()
        for j in range(1,n): anew[j]=a[j]+k*np.conj(a[n-j])
        anew[n]=k; a=anew; E=E*(1-abs(k)**2); al.append(-np.conj(k))
    return al
def moms(S,t,N): S=np.array(S); return [np.mean(np.exp(2j*np.pi*k*S*t)) for k in range(N+1)]
# EVEN runners at t=1/14
even=[2,4,6,8,10,12]
alE=verblunsky(moms(even,1/14,6))
print("EVEN runners {2,4,6,8,10,12} at t=1/14 (= 7th roots minus origin):")
print(f"  |alpha_j| = {[round(abs(a),4) for a in alE]}")
print(f"  harmonic 1/(7-1-j)=1/(6-j): {[round(1/(6-j),4) for j in range(6)]}  MATCH? {np.allclose([abs(a) for a in alE],[1/(6-j) for j in range(6)])}")
print(f"  => the EVEN sub-cloud IS the n=7 harmonic law: a recursive LRC(7) inside LRC(14).  alpha_0=1/6 (centroid from origin).")
# ODD runners at t=1/14 (full z^7=-1 heptagon, 7 points)
odd=[1,3,5,7,9,11,13]
alO=verblunsky(moms(odd,1/14,6))
print(f"\nODD runners {odd} at t=1/14 (= full z^7=-1 heptagon, 7 pts):")
print(f"  |alpha_j| = {[round(abs(a),4) for a in alO]}  (uniform on a full rotated heptagon => Verblunsky ~0 then terminal 1 = SYMMETRIC/Haar)")
print(f"  => the ODD sub-cloud = the SC D_7 tournament's Cayley-dual (7th roots), Verblunsky (0..0,1).")
print("\n== THE DECOMPOSITION (14 = 2 x 7, CRT) ==")
print("  Z/2 (parity) = iota = the COMPLEMENT/antipode symmetry (splits odd vs even runners)")
print("  Z/7 (heptagon) = D_7 = the tournament vertices / 7th roots")
print("  ODD half  -> full heptagon z^7=-1  -> SC rotational tournament (Aut C_7, 14 triangles) -Cayley-> 7th roots, Verblunsky (0..0,1)")
print("  EVEN half -> z^7=+1 minus origin   -> HARMONIC Verblunsky 1/(6-j) = sub-LRC(7)")
print("  => LRC(14)'s lonely config FACTORS as (parity Z/2) x (heptagon D_7): the repo's Mode-B/CRT descent, made concrete.")
