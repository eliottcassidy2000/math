#!/usr/bin/env python3
"""babai_cameron_n7_kps_S128c116.py -- kind-pasteur-2026-07-20-S128c116

n = 7, the decisive case for Babai-Cameron Remark 7.4.

klein-S338 settled n = 1 mod 4 (failure set EMPTY, via the unique all-even-score
member).  My n <= 6 run gives  0, 2, 0, 4  at n = 3,4,5,6 -- so BOTH odd cases so far
are empty, and the even ones are not.  If n = 7 is also empty, the conjecture is that
the failure set is empty for EVERY odd n, which is strictly stronger than klein's
theorem and would need a different argument at 3 mod 4.

TWO STRUCTURAL FACTS THAT FRAME THE ANSWER.

  * MOON: the automorphism group of a tournament has ODD order.  So if sigma has EVEN
    order it lies in no Aut(T), hence fixes no member of any class automatically.
    Therefore any even-order sigma that STABILISES a class is an instant witness, and
    the whole question at odd n is whether such a sigma can stabilise anything.
  * At n = 3 mod 4 each class has exactly n members of parity-weight one, and sigma
    sends the member with odd score at v to the one with odd score at sigma(v).  So
    sigma fixes one of those IFF sigma fixes a VERTEX.  A fixed-point-free sigma must
    therefore be caught, if at all, on the heavier parity strata.

Same exact reduction as before (cosets of the cut space; sigma stabilises x+C iff
v := Ax + c lies in C, and fixes a member iff v lies in A(C)), but vectorised with
byte tables so that 32768 cosets x 5040 permutations is seconds rather than hours.
"""
import sys
from itertools import permutations
import numpy as np

n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
eidx = {e: k for k, e in enumerate(edges)}
m = len(edges)
print("n = %d, %d edges, cosets = 2^%d" % (n, m, m - (n - 1)))

Cset = set()
for U in range(1 << n):
    d = 0
    for k, (i, j) in enumerate(edges):
        if ((U >> i) & 1) != ((U >> j) & 1):
            d |= 1 << k
    Cset.add(d)
Clist = sorted(Cset)
print("cut space |C| = %d (expect 2^%d = %d)" % (len(Clist), n - 1, 1 << (n - 1)))


def basis_of(vecs):
    basis = []
    for v in vecs:
        cur = v
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
    return basis


Cb = basis_of(Clist)
comp = []
for k in range(m):
    cur = 1 << k
    for b in Cb + comp:
        cur = min(cur, cur ^ b)
    if cur:
        comp.append(cur)
reps = np.zeros(1, dtype=np.int64)
for b in comp:
    reps = np.concatenate([reps, reps ^ b])
print("coset representatives: %d" % len(reps))

NB = (m + 7) // 8
Carr = np.array(Clist, dtype=np.int64)


def perm_tables(sigma):
    inv = [0] * n
    for i, s in enumerate(sigma):
        inv[s] = i
    Ptab = [0] * m
    c = 0
    for k, (a, b) in enumerate(edges):
        sa, sb = sigma[a], sigma[b]
        Ptab[k] = eidx[(min(sa, sb), max(sa, sb))]
        if inv[a] > inv[b]:
            c |= 1 << k
    tabs = []
    for byte in range(NB):
        t = np.zeros(256, dtype=np.int64)
        for val in range(256):
            y = 0
            for bit in range(8):
                k = byte * 8 + bit
                if k < m and (val >> bit) & 1:
                    y |= 1 << Ptab[k]
            t[val] = y
        tabs.append(t)
    return tabs, c, Ptab


def applyP_vec(x, tabs):
    out = np.zeros_like(x)
    for byte in range(NB):
        out ^= tabs[byte][(x >> (8 * byte)) & 255]
    return out


bad = set()
witness = {}
cnt_even_order_stab = 0
for sigma in permutations(range(n)):
    if list(sigma) == list(range(n)):
        continue
    tabs, c, Ptab = perm_tables(sigma)
    AC = set()
    for u in Clist:
        y = 0
        uu = u
        while uu:
            b = uu & -uu
            y |= 1 << Ptab[b.bit_length() - 1]
            uu ^= b
        AC.add(y ^ u)
    diff = Cset - AC
    if not diff:
        continue
    v = applyP_vec(reps, tabs) ^ reps ^ c
    darr = np.array(sorted(diff), dtype=np.int64)
    hit = np.isin(v, darr)
    if hit.any():
        # order of sigma
        o, p = 1, list(sigma)
        cur = list(sigma)
        while cur != list(range(n)):
            cur = [sigma[i] for i in cur]
            o += 1
        for x in reps[hit].tolist():
            bad.add(x)
            witness.setdefault(x, (sigma, o))

print()
print("=" * 70)
print("RESULT at n = %d" % n)
print("=" * 70)
print("  labelled switching classes that FAIL (some automorphism fixes no member): %d"
      % len(bad))
if bad:
    ex = sorted(bad)[:3]
    for x in ex:
        s, o = witness[x]
        print("     coset rep %-10d witness sigma = %s (order %d)" % (x, s, o))
    orders = sorted({witness[x][1] for x in bad})
    print("  witness permutation orders seen : %s" % orders)
    print("  any EVEN-order witness (Moon-automatic) : %s"
          % any(o % 2 == 0 for o in orders))
else:
    print("  EMPTY.")
    print("  So n = %d joins n = 3 and n = 5: the failure set is empty at every odd n" % n)
    print("  computed so far, which is STRICTLY STRONGER than klein-S338's theorem")
    print("  (that covers only n = 1 mod 4) and needs a separate argument at 3 mod 4.")
    print()
    print("  Moon's theorem supplies half of it: an even-order permutation lies in no")
    print("  Aut(T), so it fixes no member automatically -- hence emptiness says NO")
    print("  even-order permutation stabilises any switching class at odd n.  That is")
    print("  the statement to prove.")
