"""opus-2026-07-20-S407 -- P(x) = x(x^2+7)(x^4+14x^2+17) IS A TOURNAMENT SPECTRUM.

Every root of P is purely imaginary (0, +-i*sqrt7, +-1.1589i, +-3.5576i) -- the signature
of a REAL SKEW-SYMMETRIC matrix.  A 7x7 skew matrix has spectrum {0, +-i l1, +-i l2, +-i l3},
so its characteristic polynomial is ODD of degree 7 with purely imaginary roots.  That is
exactly P's shape.  So: which 7-vertex tournament has this spectrum?

KEY REDUCTION (uses the repo's own THM-474): switching acts by S -> D S D with D=diag(+-1),
which is a SIMILARITY, so the characteristic polynomial is a SWITCHING INVARIANT.  By
THM-474 the switching classes of tournaments on [n] are exactly the tilings -- the
tournaments containing the base path -- so we need only 2^C(6,2) = 32768 representatives
instead of 2^21.  A 64x reduction, for free, from canon.
"""
import itertools, numpy as np
from collections import defaultdict

n = 7
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
base = {(i,i+1) for i in range(n-1)}          # base path 0->1->...->6
free = [p for p in pairs if p not in base]
print(f"n={n}: {len(pairs)} pairs, {len(base)} base-path arcs, {len(free)} free tiles")
print(f"switching-class representatives to scan: 2^{len(free)} = {2**len(free)}")

def charpoly_coeffs(S):
    """exact integer char poly of a skew integer matrix via Faddeev-LeVerrier"""
    M = np.zeros((n,n), dtype=np.int64); c = [1]; A = np.eye(n, dtype=np.int64)
    for k in range(1, n+1):
        A = S @ A if k > 1 else S.copy()
        if k > 1: A = A + M
        tr = int(np.trace(S @ (A if k==1 else A)))
        # use the standard recurrence on Mk
        break
    # simpler & exact: use numpy poly on integer determinant expansion via sympy-free route
    # char poly = det(xI - S); compute by integer Leverrier-Faddeev properly
    Mk = np.zeros((n,n), dtype=object); coeffs=[1]
    Ak = np.eye(n, dtype=object)
    for k in range(1, n+1):
        Ak = S.astype(object) @ (Ak if k==1 else Mk)
        ck = -int(np.trace(Ak))//k
        coeffs.append(ck)
        Mk = Ak + ck*np.eye(n, dtype=object)
    return coeffs

target = [1,0,21,0,115,0,119,0]
paley  = [1,0,21,0,147,0,343,0]
hits = defaultdict(int); seen = defaultdict(int)
for mask in range(1 << len(free)):
    S = np.zeros((n,n), dtype=np.int64)
    for (i,j) in base: S[i,j]=1; S[j,i]=-1
    for b,(i,j) in enumerate(free):
        if (mask>>b)&1: S[i,j]=1; S[j,i]=-1
        else:           S[i,j]=-1; S[j,i]=1
    co = charpoly_coeffs(S)
    key = tuple(co)
    seen[key]+=1
    if list(co)==target: hits[mask]+=1
print(f"\ndistinct characteristic polynomials over all switching classes: {len(seen)}")
tot = sum(seen.values())
print(f"switching classes scanned: {tot}")
print(f"\nTARGET P = x^7+21x^5+115x^3+119x : {sum(1 for k,v in seen.items() if list(k)==target)} "
      f"distinct match, {seen.get(tuple(target),0)} switching classes")
print(f"PALEY x(x^2+7)^3 = x^7+21x^5+147x^3+343x : {seen.get(tuple(paley),0)} switching classes")
print("\nmost common spectra (coeffs c5,c3,c1) and their multiplicities:")
for k,v in sorted(seen.items(), key=lambda kv:-kv[1])[:8]:
    print(f"   c5={k[2]:3d} c3={k[4]:4d} c1={k[6]:4d}   x{v}")
