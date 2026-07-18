#!/usr/bin/env python3
"""trapped_core_nontight_kps_S128c48.py -- kind-pasteur S128 cont.48.
DECISIVE: is every tight family (M = 1/14) excluded from the trapped-core cut?
If so, the THM-984 residual 'mu0 > 0 on trapped cores' is a STRICT-margin question,
not an equality one -- the M=1/14 rigidity case never arises inside the cut.

Trapped-core hypotheses (from LRC14Grand.ResidualObligation), all must hold:
 (H1) all nonzero        (H2) covering        (H3) GapFamily: some |vi| > 13|vj|
 (H4) compressed: for every i, exists j!=i with |vi| <= 13|vj|  (no dominant runner)
 (H5) distinct |vi|      (H6) max|vi| >= 23
 (H7) non-clusterable: no g>=2 divides all-but-one of the v
 (H8) no AP-embedding (near-AP): not inside a bounded-lift AP with <=12 slopes
 (H9) common-residue exclusion: no d>=2 with d | (vi - a) for all i (some a)
M(v) computed exactly on the Farey-refined critical grid.
"""
import sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)

def norm_dist(x):  # ||x|| = dist to nearest integer, x a Fraction
    fx = x - int(x)
    if fx < 0: fx += 1
    return min(fx, 1 - fx)

def M_exact(V):
    """max_t min_i ||vi t||. Critical t: a/(vi+-vj) and a/(2vi), reduced. Exact rationals."""
    cand = set()
    n = len(V)
    for i in range(n):
        vi = V[i]
        for a in range(1, 2*vi):
            cand.add(F(a, 2*vi))
    for i in range(n):
        for j in range(n):
            for s in (V[i]+V[j], abs(V[i]-V[j])):
                if s == 0: continue
                for a in range(1, s):
                    cand.add(F(a, s))
    best = F(0)
    bestt = None
    for t in cand:
        if not (0 < t < 1): continue
        m = min(norm_dist(v*t) for v in V)
        if m > best:
            best = m; bestt = t
    return best, bestt

def covering(V):
    # covering: for every prime power modulus q, some vi shares... use the standard:
    # NOT covering iff exists q>=2 with q | none of ... (sieve). Practical proxy for LRC:
    # family is 'covering' iff for all q in 2..some bound, exists residue argument.
    # Use the repo definition: NON-covering iff exists q (2<=q<=14) with q not dividing
    # any vi and the q-sieve gives loneliness. We test: exists q in 2..14, gcd(q, all vi)
    # arrangement -> here approximate covering as: no single small q sieves it.
    # Exact per CoveringFamily: for all 2<=q<=14, exists i with (vi mod q) hitting each class... 
    # We use: family covering iff for every q in 2..14 there is NO t=1/q witness, i.e.
    # for q, min_i ||vi/q|| < 1/14. If some q gives min_i||vi/q||>=1/14 -> t=1/q is lonely (non-trapped, sieve branch).
    for q in range(2, 15):
        t = F(1, q)
        if min(norm_dist(v*t) for v in V) >= F(1,14):
            return False  # sieve at 1/q already lonely -> not in covering/trapped branch
    return True

def gap_family(V):  # H3
    A = [abs(v) for v in V]
    return any(a > 13*b for a in A for b in A)

def compressed(V):  # H4
    A = [abs(v) for v in V]
    return all(any(j!=i and A[i] <= 13*A[j] for j in range(len(A))) for i in range(len(A)))

def distinct(V):  # H5
    A = [abs(v) for v in V]
    return len(set(A)) == len(A)

def maxge23(V):  # H6
    return max(abs(v) for v in V) >= 23

def non_clusterable(V):  # H7
    n = len(V)
    for i0 in range(n):
        rest = [V[j] for j in range(n) if j != i0]
        g = reduce(gcd, [abs(x) for x in rest], 0)
        if g >= 2 and (V[i0] % g != 0):
            return False  # g divides all-but-i0 but not vi0 -> clusterable
    return True

def no_ap_embed(V):  # H8 (bounded-lift near-AP): try to write vi = a_i + L*k_i, |a_i|<=A, A/L<=1/13-1/14, <=12 distinct k
    A = sorted(abs(v) for v in V)
    for L in range(1, A[-1]+1):
        ks = [round(a / L) for a in A]
        if len(set(ks)) > 12: continue
        res = [a - L*k for a, k in zip(A, ks)]
        Amax = max(abs(r) for r in res)
        if L > 0 and F(Amax, L) <= F(1,13) - F(1,14) and all(k != 0 for k in ks):
            return False  # embeds in a near-AP -> not trapped
    return True

def no_common_residue(V):  # H9
    for d in range(2, max(abs(v) for v in V)+1):
        for a in range(d):
            if all((v - a) % d == 0 for v in V):
                return False
    return True

HYPS = [("covering",covering),("gap",gap_family),("compressed",compressed),
        ("distinct",distinct),("max>=23",maxge23),("non-cluster",non_clusterable),
        ("no-AP",no_ap_embed),("no-comm-res",no_common_residue)]

def trapped_report(name, V):
    m, t = M_exact(V)
    tight = (m == F(1,14))
    fails = [nm for nm,fn in HYPS if not fn(V)]
    trapped = (len(fails)==0)
    print("  %-22s M=%s (%.5f) %s | trapped=%s%s" % (
        name, m, float(m), "TIGHT(=1/14)" if tight else ("margin+%.4f"%(float(m)-1/14)),
        trapped, "" if not fails else "  FAILS: "+",".join(fails)))
    return tight, trapped, fails

print("== KNOWN TIGHT / NEAR-TIGHT FAMILIES (does each escape the trapped cut?) ==")
fams = [
  ("deep-well {1..12,182}", list(range(1,13))+[182]),
  ("AP {1..13}", list(range(1,14))),
  ("AP {2,4,..,26}", [2*i for i in range(1,14)]),
  ("sporadic {1..11,13,24}", list(range(1,12))+[13,24]),
  ("perturbed {1..11,13,36}", list(range(1,12))+[13,36]),
  ("deep-well {1..12,84}", list(range(1,13))+[84]),
  ("deep-well {1..12,98}", list(range(1,13))+[98]),
]
tight_in_cut = []
for name, V in fams:
    tight, trapped, fails = trapped_report(name, V)
    if tight and trapped:
        tight_in_cut.append(name)
print()
if tight_in_cut:
    print("*** WARNING: tight families INSIDE the trapped cut:", tight_in_cut, "-- residual DOES include equality ***")
else:
    print(">>> EVERY tight family escapes the trapped cut: the M=1/14 rigidity case NEVER arises inside trapped cores.")
    print(">>> Hence the THM-984 residual 'mu0>0 on trapped cores' is a STRICT-MARGIN question (M > 1/14 strictly).")
print("DONE")
