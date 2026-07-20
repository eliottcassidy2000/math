#!/usr/bin/env python3
"""babai_cameron_failure_set_kps_S128c116.py -- kind-pasteur-2026-07-20-S128c116

BABAI-CAMERON REMARK 7.4, IN THE RANGE THEIR STRUCTURE THEOREM DOES NOT COVER.

Babai-Cameron (EJC 7 (2000) #R38) Remark 7.4 asks to enumerate the switching classes
of tournaments admitting an automorphism that fixes NO member of the class, and says
"We cannot do this."

klein-S338 (THM-1440) settled n = 1 mod 4: every switching class then contains a
UNIQUE all-even-score tournament, any automorphism of the class must fix it, so the
failure set is EMPTY.  That argument uses n = 1 mod 4 essentially (it is where
sum s_v = C(n,2) is even), and says nothing at n = 3 mod 4 or n even.  Those are
exactly the cases left, and n = 3, 7 are the smallest.

THE REDUCTION USED HERE, which turns the question into linear algebra over F_2.
Fix the reference orientation i -> j for i < j and encode a tournament by the set
x of edges it REVERSES.  Then

  * switching at U reverses the cut delta(U), so SWITCHING CLASSES ARE THE COSETS
    x + C of the cut space C = span{delta(U)}, dim C = n - 1;
  * a permutation sigma acts AFFINELY, x |-> P x + c, where P permutes edge
    coordinates and c is the inversion set of sigma^{-1};
  * writing A := P + I, for v := A x + c one has
        sigma STABILISES the class x + C   <=>  v in C
        sigma FIXES SOME MEMBER of x + C   <=>  v in A(C)
    (sigma fixes y iff A y = c; with y = x + u this is A u = v).
  * A(C) is contained in C, since A delta(U) = delta(sigma U) + delta(U), so both
    conditions are well defined on the class.

Hence:  sigma witnesses failure on x + C  <=>  v in C \\ A(C).

That is a membership test in a 2^{n-1}-element set, evaluated over all cosets and all
sigma, and it is exact.  Computed at n = 3..7, with n = 5 as the CONTROL that must
return zero by klein's theorem.
"""
import sys
from itertools import permutations
import numpy as np

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7


def setup(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    eidx = {e: k for k, e in enumerate(edges)}
    m = len(edges)
    return edges, eidx, m


def cut_space(n, edges, eidx, m):
    """All delta(U) as bitmasks; a set of size 2^(n-1)."""
    out = set()
    for U in range(1 << n):
        d = 0
        for k, (i, j) in enumerate(edges):
            if ((U >> i) & 1) != ((U >> j) & 1):
                d |= 1 << k
        out.add(d)
    return out


def basis_of(vecs, m):
    piv, basis = [], []
    for v in vecs:
        cur = v
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
    return basis


def coset_reps(Cset, m):
    """Representatives of F_2^m / C: span of a complement basis."""
    Cb = basis_of(list(Cset), m)
    comp = []
    for k in range(m):
        v = 1 << k
        cur = v
        for b in Cb + comp:
            cur = min(cur, cur ^ b)
        if cur:
            comp.append(cur)
    reps = [0]
    for b in comp:
        reps = reps + [r ^ b for r in reps]
    return reps, comp


def perm_data(sigma, n, edges, eidx):
    """P as an edge permutation table, and c = inversion set of sigma^{-1}."""
    inv = [0] * n
    for i, s in enumerate(sigma):
        inv[s] = i
    Ptab = [0] * len(edges)
    c = 0
    for k, (a, b) in enumerate(edges):
        sa, sb = sigma[a], sigma[b]
        Ptab[k] = eidx[(min(sa, sb), max(sa, sb))]
        # c bit at edge e={a,b}, a<b : 1 iff sigma^{-1} reverses the order
        if inv[a] > inv[b]:
            c |= 1 << k
    return Ptab, c


def apply_P(x, Ptab, m):
    y = 0
    xx = x
    while xx:
        b = xx & -xx
        k = b.bit_length() - 1
        y |= 1 << Ptab[k]
        xx ^= b
    return y


print("=" * 78)
print("BABAI-CAMERON REMARK 7.4 FAILURE SET")
print("  switching classes admitting an automorphism that fixes NO member")
print("=" * 78)
print("  %-4s %-8s %-12s %-16s %-14s" % ("n", "n mod 4", "sw classes", "labelled bad", "bad up to iso"))
for n in range(3, NMAX + 1):
    edges, eidx, m = setup(n)
    Cset = cut_space(n, edges, eidx, m)
    reps, comp = coset_reps(Cset, m)
    Clist = sorted(Cset)
    nclass_lab = len(reps)
    bad = set()
    perms = [p for p in permutations(range(n)) if list(p) != list(range(n))]
    for sigma in perms:
        Ptab, c = perm_data(sigma, n, edges, eidx)
        AC = set()
        for u in Clist:
            AC.add(apply_P(u, Ptab, m) ^ u)
        diff = Cset - AC
        if not diff:
            continue
        for x in reps:
            v = apply_P(x, Ptab, m) ^ x ^ c
            if v in diff:
                bad.add(x)
    # count bad classes up to isomorphism
    canon = set()
    for x in bad:
        best = None
        for sigma in permutations(range(n)):
            Ptab, c = perm_data(sigma, n, edges, eidx)
            y = apply_P(x, Ptab, m) ^ c
            mn = min(y ^ u for u in Clist)
            if best is None or mn < best:
                best = mn
        canon.add(best)
    # switching classes up to iso, for context
    seen = set()
    for x in reps:
        best = None
        for sigma in permutations(range(n)):
            Ptab, c = perm_data(sigma, n, edges, eidx)
            y = apply_P(x, Ptab, m) ^ c
            mn = min(y ^ u for u in Clist)
            if best is None or mn < best:
                best = mn
        seen.add(best)
    print("  %-4d %-8d %-12d %-16d %-14d" % (n, n % 4, len(seen), len(bad), len(canon)))

print()
print("  A049313 (Babai-Cameron, switching classes of tournaments up to iso):")
print("     n=3..7 should read 1, 2, 2, 6, 12")
print()
print("  CONTROL: n = 5 is 1 mod 4, so klein-S338's theorem forces the failure set")
print("  to be EMPTY there.  A non-zero entry at n = 5 would refute that theorem and")
print("  would mean this script is wrong -- it is the instrument check, not a result.")
print()
print("  The open cases are n = 3 mod 4 (n = 3, 7) and n even, which is exactly where")
print("  klein's parity argument gives nothing: at n = 3 mod 4 the score-parity vector")
print("  has ODD weight, so there is no all-even member to be forced to stay put --")
print("  instead there are n members with a single odd score, and an automorphism is")
print("  free to permute them without fixing any.")
