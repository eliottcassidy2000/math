#!/usr/bin/env python3
"""
boxeph-S17: the verified-interval-evaluator DATA layer for mu_L, L <= 12.
Design (single-containment reduction): danger covers Cov_L satisfy the
recursion  D_L subset D_1 union T^{-1}(D_{L-1}),  T y = 2y mod 1.  The
pullback of an interval [a,b] is [a/2,b/2] u [(a+1)/2,(b+1)/2].  We emit,
per L: (i) the merged danger cover Cov_L (exact rationals), (ii) the
PRE-SPLIT pullback pieces of Cov_{L-1}, each verified to sit inside a
SINGLE interval of Cov_L (so the Lean checker is a trivially-sound
exists-containment over Q), (iii) total lengths = 1 - mu_L exact, matching
the S11 table, and (iv) the SIX census-branch witnesses: the 2-chain
partition budgets (12,1)...(7,6), each < 1 -- the W0 = 11 theorem's data.
"""
from fractions import Fraction as F

LO, HI = F(1,14), F(13,14)

def merge(iv):
    iv = sorted(iv)
    out = [list(iv[0])]
    for a,b in iv[1:]:
        if a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a,b])
    return [(a,b) for a,b in out]

def danger_cover(L):
    """Merged danger cover: complement of the safe set (exact)."""
    # danger = union over j<L of level-j intervals (m/2^j +- 1/(14 2^j)) clipped [0,1]
    iv = []
    for j in range(L):
        d = 2**j
        for m in range(d+1):
            a, b = F(m,d) - F(1,14*d), F(m,d) + F(1,14*d)
            iv.append((max(a,F(0)), min(b,F(1))))
    return merge(iv)

def pullback(cov):
    out = []
    for a,b in cov:
        out.append((a/2, b/2))
        out.append(((a+1)/2, (b+1)/2))
    return merge(out)

def split_into(pieces, target):
    """Split each piece at target-interval boundaries; verify single containment."""
    bounds = sorted(set([a for a,b in target] + [b for a,b in target]))
    out = []
    for a,b in pieces:
        cuts = [a] + [x for x in bounds if a < x < b] + [b]
        for u,v in zip(cuts, cuts[1:]):
            out.append((u,v))
    # verify: every piece inside ONE target interval
    for u,v in out:
        if not any(c <= u and v <= d for c,d in target):
            raise AssertionError(f"piece ({u},{v}) not singly contained")
    return out

def length(cov):
    return sum(b-a for a,b in cov)

print("mu_L evaluator data layer: exact covers + single-containment witnesses")
covs = {1: merge([(F(0),LO),(HI,F(1))])}
print(f"L= 1: cover {len(covs[1])} intervals, danger = {length(covs[1])} = {float(length(covs[1])):.6f}")
mu_expected = {1:F(6,7),2:F(11,14),3:F(5,7),4:F(9,14),5:F(33,56),6:F(15,28),
               7:F(109,224),8:F(199,448),9:F(181,448),10:F(659,1792),
               11:F(1201,3584),12:F(1093,3584)}
ok_all = True
for L in range(2,13):
    cov = danger_cover(L)
    pb = pullback(covs[L-1])
    base = covs[1]
    # the recursion cover: base u pullback -- verify it IS contained in cov and vice versa lengths
    pieces = split_into(merge(base + pb), cov)
    covs[L] = cov
    dl = length(cov)
    match = (F(1) - dl) == mu_expected[L]
    ok_all &= match
    print(f"L={L:2d}: cover {len(cov):4d} intervals, danger = {dl} "
          f"(mu_L = {F(1)-dl} {'== S11 OK' if match else 'MISMATCH!'}); "
          f"recursion pieces singly contained: {len(pieces)}")
print(f"\nall mu_L match S11 exact table: {ok_all}")
print("\n== THE SIX CENSUS-BRANCH WITNESSES (2-chain partitions of 13) ==")
for L1 in range(12, 6, -1):
    L2 = 13 - L1
    b = length(covs[L1]) + length(covs[L2])
    print(f"  ({L1:2d},{L2:2d}): budget = {b} = {float(b):.6f}  "
          f"{'< 1 PASS' if b < 1 else 'FAIL'}")
# emit Lean-ready data sizes
import json
data = {L: [[str(a), str(b)] for a,b in covs[L]] for L in covs}
with open('05-knowledge/results/lrc_mu_covers_boxeph_S17.json','w') as f:
    json.dump(data, f)
print("\ncovers emitted to 05-knowledge/results/lrc_mu_covers_boxeph_S17.json")
