#!/usr/bin/env python3
"""gmc2_nullcone_kps_S128c118.py -- kind-pasteur-2026-07-20-S128c118

GMC(2): TWO REAL GAUSSIANS = ONE COMPLEX GAUSSIAN Z.

Target: the NULLCONE
    N = { P in C[Z, conj(Z)] : E[P^m] = 0 for all m >= 1 }
and the conjecture that it is a Mathieu subspace, i.e. for every P in N and every Q,
E[Q P^m] = 0 for all large m.  Proving a STRUCTURE THEOREM for N settles GMC(2)
outright, which is why it is the stronger statement to aim at.

THE CHARGE GRADING IS THE WHOLE PICTURE.  Writing Z = r e^{i phi}, the monomial
Z^a conj(Z)^b is r^{a+b} e^{i(a-b) phi}, so a-b is a FOURIER MODE in phi.  The Gaussian
measure has phi uniform and independent of r, so E kills every monomial of non-zero
charge:  E[Z^a conj(Z)^b] = delta(a,b) a!.

That gives a free family of nullcone elements: if every monomial of P has charge > 0
(or every one has charge < 0), then every monomial of P^m has charge >= m*k_min > 0,
so E[P^m] = 0 for all m >= 1 automatically.  Call these ONE-SIDED.  For a one-sided P
the Mathieu property is immediate too: E[Q P^m] can only be non-zero if some monomial
of Q has charge <= -m*k_min, which fails once m > deg(Q)/|k_min|.  So

    STRUCTURE CONJECTURE (2D nullcone):  N = {one-sided P} u {0}.

and GMC(2) follows from it in one line.

THIS SCRIPT TESTS THE CONJECTURE, i.e. hunts for a MIXED-charge element of N.
Lesson from the previous session: a sieve at m <= 3 means nothing here, so candidates
must survive to serious depth before being believed, and any survivor must then be
checked for the SECOND half of the definition (some Q with non-vanishing moments).
"""
import sys
from math import factorial
from itertools import product

DEG = int(sys.argv[1]) if len(sys.argv) > 1 else 3
MM = int(sys.argv[2]) if len(sys.argv) > 2 else 9
COEF = [-1, 0, 1]

MONS = [(a, b) for a in range(DEG + 1) for b in range(DEG + 1) if 1 <= a + b <= DEG]
MONS = [(0, 0)] + MONS
print("monomials Z^a conj(Z)^b with a+b <= %d : %d of them" % (DEG, len(MONS)))


def mul(p, q):
    out = {}
    for (a, b), u in p.items():
        for (c, d), v in q.items():
            k = (a + c, b + d)
            out[k] = out.get(k, 0) + u * v
    return {k: v for k, v in out.items() if v}


def Eval(p):
    return sum(v * factorial(a) for (a, b), v in p.items() if a == b)


def in_nullcone(p, M):
    cur = dict(p)
    for m in range(1, M + 1):
        if Eval(cur) != 0:
            return m
        cur = mul(cur, p)
    return None


def charges(p):
    return sorted({a - b for (a, b) in p})


def one_sided(p):
    cs = charges(p)
    return all(c > 0 for c in cs) or all(c < 0 for c in cs)


print("searching coefficient vectors in %s^%d ..." % (COEF, len(MONS)))
found_mixed = []
n_null = 0
n_onesided = 0
tested = 0
for vec in product(COEF, repeat=len(MONS)):
    if all(v == 0 for v in vec):
        continue
    p = {MONS[i]: vec[i] for i in range(len(MONS)) if vec[i]}
    tested += 1
    # cheap pre-filter: E[P] and E[P^2]
    if Eval(p) != 0:
        continue
    p2 = mul(p, p)
    if Eval(p2) != 0:
        continue
    fail = in_nullcone(p, MM)
    if fail is not None:
        continue
    n_null += 1
    if one_sided(p):
        n_onesided += 1
    else:
        found_mixed.append(p)

print()
print("=" * 78)
print("RESULT")
print("=" * 78)
print("  coefficient vectors tested          : %d" % tested)
print("  survive E[P^m] = 0 for m = 1..%-3d   : %d" % (MM, n_null))
print("  of those, ONE-SIDED in charge       : %d" % n_onesided)
print("  of those, MIXED charge              : %d" % len(found_mixed))
if found_mixed:
    print()
    print("  MIXED-charge survivors (the conjecture would be FALSE if any is genuine):")
    for p in found_mixed[:10]:
        s = " + ".join("%s*Z^%d Zb^%d" % (v, a, b) for (a, b), v in sorted(p.items()))
        print("     charges %s :  %s" % (charges(p), s))
    print()
    print("  Each needs the SECOND test: is there a Q with E[Q P^m] non-zero for large m?")
    QS = [((1, 0), "Z"), ((0, 1), "Zb"), ((1, 1), "|Z|^2"), ((2, 0), "Z^2"),
          ((0, 2), "Zb^2"), ((2, 1), "Z^2 Zb"), ((0, 0), "1"), ((3, 0), "Z^3")]
    for p in found_mixed[:10]:
        print("     P with charges %s :" % charges(p))
        cur = dict(p)
        pows = []
        for m in range(1, MM + 1):
            pows.append(dict(cur))
            cur = mul(cur, p)
        for qm, qn in QS:
            q = {qm: 1}
            vals = [Eval(mul(q, pw)) for pw in pows]
            nz = [m for m, v in enumerate(vals, 1) if v != 0]
            print("        Q = %-8s E[QP^m] = %-34s last non-zero: %s"
                  % (qn, str(vals[:7]), max(nz) if nz else None))
else:
    print()
    print("  NO mixed-charge element of the nullcone exists in this box.")
    print("  Consistent with the structure conjecture N = {one-sided} u {0}, and hence")
    print("  with GMC(2) being TRUE.  This is a search, so it bounds the counterexample")
    print("  rather than proving the theorem -- the detection floor is: all coefficient")
    print("  vectors in {-1,0,1} over %d monomials (degree <= %d), depth m <= %d."
          % (len(MONS), DEG, MM))
