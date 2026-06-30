"""
M(S) vs the descent (cusp/off-cusp). Verify M>=1/14 across both. Explore the renormalization map
S -> E/2 (the odd-core trajectory) and the doublet attractor.
"""
import math, cmath
from fractions import Fraction
w=cmath.exp(2j*cmath.pi/7)
def g(res):
    return min(abs(sum(w**((k*x)%7) for x in res))**2 for k in range(1,7))
def descent(S):
    ch=[]; cur=list(S)
    while cur:
        O=[x for x in cur if x%2==1]; E=[x//2 for x in cur if x%2==0]
        if O: ch.append(O)
        cur=E
        if len(ch)>25: break
    return ch
def cusp_levels(S):
    return [j for j,O in enumerate(descent(S)) if set(x%7 for x in O)==set(range(7))]
def M_exact(S):
    # M = max over rational t=a/q of min_i ||s_i a/q||; optimal q <= some bound. Use q up to 400.
    best=Fraction(0)
    Smax=max(S)
    for q in range(2, 420):
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            m=min(min((s*a)%q, q-(s*a)%q) for s in S)
            val=Fraction(m,q)
            if val>best: best=val
        # early structural: M is rational with small q typically
    return best
sets={
 "AP {1..13} (cusp,extremal)": list(range(1,14)),
 "{1..11,13,84} (cusp,cover)": list(range(1,12))+[13,84],
 "{1..12,182} (offcusp,cover)": list(range(1,13))+[182],
 "{2..14} (offcusp,cover)": list(range(2,15)),
 "cusp-cover A": [1,3,5,7,9,11,13,2,4,6,8,10,84],   # Z7 odds + evens, 84=7.12 for q=12
 "cusp-cover B": [1,3,5,7,9,11,13,4,6,8,10,14,12],   # Z7 odds, has 14
}
print(f"{'set':>30} {'cusp?':>8} {'M=':>9} {'>=1/14?':>8} {'1/14=0.07143'}")
for name,S in sets.items():
    cl=cusp_levels(S); M=M_exact(S)
    print(f"{name:>30} {str(cl):>8} {str(M):>9} {float(M)>=1/14-1e-9!s:>8}  M={float(M):.5f}")
print()
# renormalization: the odd-core trajectory mod 7, and the doublet attractor
print("RENORMALIZATION map S -> E/2 (peel odds, halve evens). Odd-core trajectory mod 7 for AP:")
ch=descent(list(range(1,14)))
traj=[sorted(set(x%7 for x in O)) for O in ch]
for j,(O,t) in enumerate(zip(ch,traj)):
    print(f"   level {j}: |O|={len(O)} mod7={t} g={g([x%7 for x in O]):.4f}")
print("   => trajectory Z_7 -> {1,3,5} -> {1,3} -> {1}: shrinks toward the DOUBLET then SINGLETON.")
print("   the doublet {1,3} (g=0.198) is the last NONTRIVIAL core -- the attractor carrying the floor.")
print()
print("IDEA: the descent is a RENORMALIZATION. Each level halves the even sublattice (2-adic zoom).")
print("  The odd core mod 7 is the 'apex state'. Z_7 = the cusp (measure 0). Proper cores g>0 (measure>0).")
print("  Trajectory always contracts to {singleton} (g=1); the last doublet sets the per-level floor 0.198.")
