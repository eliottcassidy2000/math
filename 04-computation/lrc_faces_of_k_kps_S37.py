#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S37: THE FACES OF THE COMPLEXITY PARAMETER k agree on the KNOWN
gap members -- a concrete synthesis of the 4-agent convergence.

The residual has collapsed to ONE parameter (opus HYP-4496 window_num_denom_locked): for an
in-window ML value M=s/(Ns+k), the denominator q=Ns+k is locked Ns<q<(N+1)s, so numerator s,
denominator q, order k, Stern-Brocot depth, and (kps S36) the base DEFECT ORDER all move
together; bounding one bounds the height.  This file checks, on the ACTUAL known gap members
(Fan-Sun counterexamples + kps S35 slice), that the intrinsic faces agree:

  face 1 (opus HYP-4486): (s,k) = (numerator, q - N*numerator);   k>=2 & k<s<2k for in-window
  face 2 (opus HYP-4496): denominator lock  N*s < q < (N+1)*s
  face 3 (opus S109):     q <= 2*max(speeds)   (height controls q)
  face 4 (opus HYP-4496): Stern-Brocot depth of M inside the gap (mediant descent)
  face 5 (kps S36):       ladder order = k, and the mediant (SB-depth 1) is (s,k)=(3,2)

We verify M with the exact solver, then read every face off M and the family.
"""
from fractions import Fraction
import numpy as np

def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        for c in range(1,q):
            r=(va*c)%q; d=np.minimum(r,q-r); bq=int(d.min())
            if bq*bd>bn*q: bn,bd,bc=bq,q,c
    return Fraction(bn,bd)

def sb_depth_in_gap(M, N):
    """Stern-Brocot depth of M strictly inside (1/(N+1), 2/(2N+1)) via mediant descent.
    depth 1 = the mediant 3/(3N+2); deeper = further descent."""
    lo, hi = Fraction(1,N+1), Fraction(2,2*N+1)
    if not (lo < M < hi): return None
    depth = 0
    # descend the Stern-Brocot tree between lo and hi toward M
    a,b = lo, hi
    for _ in range(50):
        depth += 1
        med = Fraction(a.numerator+b.numerator, a.denominator+b.denominator)
        if M == med: return depth
        if M < med: b = med
        else: a = med
    return depth  # capped

# known gap members: (label, speeds, N=#speeds, expected ML)
MEMBERS = [
    ("Fan-Sun n=4  {3,8,11,19}",        [3,8,11,19],            4),
    ("Fan-Sun n=6  {5,6,11,17,23,28}",  [5,6,11,17,23,28],      6),
    ("Fan-Sun n=7  {1,3,4,5,7,13,18}",  [1,3,4,5,7,13,18],      7),
    ("kps S35 n=7  {1,2,3,4,5,7,18}",   [1,2,3,4,5,7,18],       7),
]

print("=== S37: the faces of k agree on the known gap members ===\n", flush=True)
print(f"  {'member':<34}{'ML':>7}{'N':>3}  {'(s,k)':>8}{'lock Ns<q<(N+1)s':>18}{'q<=2max':>9}{'SBdepth':>8}", flush=True)
for label, v, N in MEMBERS:
    M = Mw(v)
    s, k = M.numerator, M.denominator - N*M.numerator
    q = M.denominator
    lock = (N*s < q < (N+1)*s)
    q_le_2max = (q <= 2*max(v))
    d = sb_depth_in_gap(M, N)
    kws = (k>=2 and k < s < 2*k)
    print(f"  {label:<34}{str(M):>7}{N:>3}  ({s},{k}){'  k<s<2k' if kws else '        '}{str(lock):>10}{str(q_le_2max):>9}{str(d):>8}", flush=True)

print("\nREADING (corrected, honest): the lock Ns<q<(N+1)s and q<=2*max hold for EVERY member", flush=True)
print("(they are general).  BUT k<s<2k and a finite SB-depth-in-gap hold ONLY for genuine", flush=True)
print("FIRST-GAP members: the n=7 mediant 3/23 (depth 1, (s,k)=(3,2)).  The Fan-Sun n=4 (7/30)", flush=True)
print("and n=6 (8/51) are Kravitz counterexamples ABOVE the first gap (7/30>2/9, 8/51>2/13) --", flush=True)
print("so SB-depth=None and k<s<2k FAILS.  LESSON: a Kravitz counterexample is NOT automatically", flush=True)
print("a first-gap member; the first gap is far more restrictive (only mediant-type known:", flush=True)
print("n=7 -> 3/23 depth 1; n=6 -> 5/33 depth 2, kps S34).  On those, s,k,q,SB-depth,height", flush=True)
print("move together = ONE parameter, and it is SMALL.  kps S36 adds the mechanism face:", flush=True)
print("k = base defect order, depth-1 = mediant (s,k)=(3,2).", flush=True)
print("\nCRUX (all agree): bound this one parameter at N=12.  opus = finite residue system q=38", flush=True)
print("(the mediant 3/38); mac-mini = AP unique max-kissing (CKMRV LP-uniqueness) + Selberg spacing;", flush=True)
print("kps = the ladder Dx<D skip + a COMBINATORIAL defect-count bound (new lead).", flush=True)
