"""ALGEBRA<->GEOMETRY dictionary. Cayley transform S(skew tournament)->U=Cayley(S) in SO(n), eigenvalues on the
unit circle. Test: complement (S->-S) = complex conjugation/inversion of the spectrum (a 2nd geometry of the
staircase anti-diagonal reflection R=complement); SC <=> conjugation-symmetric spectrum. Also c3 & score links."""
import numpy as np, itertools
def tours(n):
    pairs=list(itertools.combinations(range(n),2))
    for bits in itertools.product([0,1],repeat=len(pairs)):
        A=np.zeros((n,n),int)
        for (i,j),b in zip(pairs,bits):
            if b==0: A[i,j]=1
            else: A[j,i]=1
        yield A
def skew(A): return A-A.T
def cayley_angles(A):
    S=skew(A).astype(float); n=A.shape[0]
    U=np.linalg.solve(np.eye(n)+S, np.eye(n)-S)   # (I+S)^{-1}(I-S) = Cayley
    ev=np.linalg.eigvals(U)
    ang=np.sort(np.round(np.angle(ev)/np.pi,4))     # angles/pi on the circle
    return ang
def canon(A,perms,n): 
    return min(tuple(A[np.ix_(p,p)].flatten()) for p in perms)
for n in [3,4,5]:
    perms=[list(p) for p in itertools.permutations(range(n))]
    # pick one rep per iso class, test the links
    seen={}
    for A in tours(n):
        c=canon(A,perms,n)
        if c not in seen: seen[c]=A.copy()
    print(f"\n=== n={n}: {len(seen)} iso classes ===")
    n_sc=0; comp_conj=True; sc_symm=True
    for c,A in seen.items():
        ang=cayley_angles(A)
        Aop=A.T                      # complement/transpose = reverse arcs
        angop=cayley_angles(Aop)
        # complement <=> conjugate spectrum (angle -> -angle)
        if not np.allclose(np.sort(-ang), angop, atol=1e-6): comp_conj=False
        # SC?
        isSC = (canon(Aop,perms,n)==c)
        if isSC:
            n_sc+=1
            # SC => spectrum symmetric under conjugation (set closed under angle->-angle)
            if not np.allclose(np.sort(-ang), ang, atol=1e-6): sc_symm=False
    print(f"  COMPLEMENT (S->-S) == complex CONJUGATION of Cayley spectrum (angle->-angle)? {comp_conj}")
    print(f"  SC classes: {n_sc}; every SC class has a CONJUGATION-SYMMETRIC spectrum (real axis mirror)? {sc_symm}")
    # example angles
    ex=list(seen.values())[0]; print(f"  example spectrum angles/pi: {cayley_angles(ex)}")
# c3 and score spectral links (one n)
print("\n=== c3 & score links (n=5) ===")
perms=[list(p) for p in itertools.permutations(range(5))]
import math
for A in list(tours(5))[:1]: pass
seen={}
for A in tours(5):
    c=canon(A,perms,5)
    if c not in seen: seen[c]=A.copy()
for c,A in list(seen.items())[:4]:
    sc=A.sum(1)                              # scores
    c3=math.comb(5,3)-sum(math.comb(int(s),2) for s in sc)   # #3-cycles (Kendall)
    S=skew(A).astype(float); tr2=-0.5*np.trace(S@S)          # = #edges = C(n,2) (skew)
    print(f"  scores={list(sc)}  c3={c3}  (c3 = C(n,3)-sum C(s_i,2) = the GRID-TRIANGLE cyclic count)")
