from math import gcd, lcm, ceil, floor
from sympy import nextprime

def clears_at(V, q):
    """Does V clear at modulus q? i.e. exists rotation c (gcd(c,q)=1) with all
    (v_i*c)%q in the safe band [ceil(2q/25), floor(23q/25)] -> M >= 2/25 at t=c/q."""
    lo = ceil(2*q/25); hi = floor(23*q/25)
    if lo > hi: return None  # band empty, q too small to certify 2/25
    for c in range(1, q):
        if gcd(c, q) != 1: continue
        if all(lo <= (v*c) % q <= hi for v in V):
            return c
    return False

def first_clearing_q(V, qmax):
    for q in range(6, qmax+1):
        r = clears_at(V, q)
        if r not in (None, False):
            return q, r
    return None, None

def escapes_upto(V, Q0):
    """True if V clears at NO q in [6, Q0]."""
    for q in range(6, Q0+1):
        r = clears_at(V, q)
        if r not in (None, False):
            return False, (q, r)
    return True, None

print("=== INDEPENDENT AUDIT of mac-mini-S36 escape families (opus-S130) ===\n")
AP = list(range(1,13))
print(f"AP={AP}: clears at any q<=200? ", first_clearing_q(AP, 200), " (expect None,None = tight locus)")
print()

for Q0 in [12, 25, 32, 37]:
    L = lcm(*range(2, Q0+1))
    p = nextprime(Q0)
    # varying-k pattern (1,2,1,2,...) all positive, non-translate
    k = [1 + (i % 2) for i in range(12)]
    V = [ (i+1) + L*k[i] for i in range(12) ]
    ratio = max(V)/min(V)
    esc, hit = escapes_upto(V, Q0)
    qcl, ccl = first_clearing_q(V, 3*p)
    # covering lower bound on M from the clearing modulus
    Mlb = ceil(2*qcl/25)/qcl if qcl else None
    congAP = all( (V[i] - (i+1)) % L == 0 for i in range(12) )  # V == AP mod L?
    print(f"Q0={Q0}: L=lcm(2..{Q0})={L}  nextprime={p}")
    print(f"  k(varying)={k}")
    print(f"  V (mod small, all ~L): min={min(V)} max={max(V)} ratio={ratio:.3f}  compressed(<=13)? {ratio<=13}")
    print(f"  V == AP mod L? {congAP}   (=> same residues as AP at every q|L, i.e. every q<={Q0})")
    print(f"  escapes ALL q in [6,{Q0}]? {esc}" + (f"  (clears at {hit})" if not esc else ""))
    print(f"  first clearing q (search up to {3*p}): q={qcl} (rotation c={ccl})   == nextprime? {qcl==p}")
    print(f"  covering M lower bound = ceil(2*{qcl}/25)/{qcl} = {Mlb:.5f}   > 2/25={2/25:.5f}? {Mlb>2/25 if Mlb else '-'}  => LOOSE (not a gap member)")
    print()
