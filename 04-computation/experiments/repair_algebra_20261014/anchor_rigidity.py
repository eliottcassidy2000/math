#!/usr/bin/env python3
"""Positivity of canonical numerators at -1 and rigidity of the anchor row (THM-4640 (1), last bullets).

K(()) = -1, K(t.e) = 3 K(t) + 2^e.  A representation s = (e_1 < ... < e_m) is canonical if e_1 >= 2.

Lemma (PROVED in THM-4640):  for canonical s with m >= 1,
    K(s) >= 3^(m+1) - 2^(m+2) >= 1     and     2^(e_m) <= K(s) + 3, with equality iff m = 1.
So every K has finitely many canonical representations, all with exponents <= log2(K + 3).

Part A checks the lemma on every canonical s with exponents <= EMAX_A.
Part B lists, for 3 <= a <= AMAX, every canonical representation of the anchor 2^a - 3 (by the lemma, all have
exponents <= a, so the search is complete) and reports the anchors with a second representation.
Part C reports the canonical representations of the Wilson primes 5, 13, 563 and of a few small odd K.
Part D cross-checks completeness: for every odd K < 2^(EMAX_D - 1), the representations found by the unbounded
backward recursion (bounded only by the lemma) equal those found by brute force with exponents <= EMAX_D.

Usage: python3 anchor_rigidity.py [EMAX_A=20] [AMAX=31] [EMAX_D=16]
"""
import sys
from itertools import combinations

EMAX_A = int(sys.argv[1]) if len(sys.argv) > 1 else 20
AMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 31
EMAX_D = int(sys.argv[3]) if len(sys.argv) > 3 else 16
checks = 0


def K_of(s):
    K = -1
    for e in s:
        K = 3 * K + (1 << e)
    return K


def canonical_reps(K):
    """All canonical increasing s with K(s) = K; complete by the lemma (exponents <= log2(K + 3))."""
    out = []

    def back(K, emax, tail):
        # tail: exponents already chosen (to the right), all > emax
        if K == -1:
            out.append(tuple(reversed(tail)))
        for e in range(2, emax + 1):
            r = K - (1 << e)
            if r % 3:
                continue
            P = r // 3
            # P = K(prefix) with prefix canonical (exponents >= 2, < e); P = -1 (empty) or P >= 1
            if P == -1 or P >= 1:
                back(P, e - 1, tail + [e])

    if K < -1 or K % 2 == 0:
        return out
    if K == -1:
        return [()]
    back(K, (K + 3).bit_length() - 1, [])
    return sorted(out, key=lambda s: (len(s), s))


# Part A: the lemma on all canonical s with exponents <= EMAX_A
cnt = 0
for m in range(1, EMAX_A):
    lo = 3 ** (m + 1) - 2 ** (m + 2)
    for s in combinations(range(2, EMAX_A + 1), m):
        K = K_of(s)
        assert K >= lo >= 1, (s, K, lo)
        assert (1 << s[-1]) <= K + 3, s
        assert ((1 << s[-1]) == K + 3) == (m == 1), s
        cnt += 1
checks += cnt
print(f"A  lemma K >= 3^(m+1) - 2^(m+2) >= 1 and 2^(e_m) <= K + 3 (= iff m = 1): "
      f"{cnt} canonical s with exponents <= {EMAX_A}: OK")

# Part B: anchors 2^a - 3
second = []
for a in range(2, AMAX + 1):
    K = (1 << a) - 3
    reps = canonical_reps(K)
    assert (a,) in reps
    checks += 1
    others = [s for s in reps if s != (a,)]
    for s in others:
        assert K_of(s) == K and s[-1] < a
        checks += 1
    if others:
        second.append((a, K, others))
print(f"B  anchors 2^a - 3, 2 <= a <= {AMAX}: second canonical representation exactly for a in "
      f"{[a for a, _, _ in second]}")
for a, K, others in second:
    for s in others:
        print(f"     a = {a:2d}  K = {K}  (a) ~ {s}   last exponents {a} vs {s[-1]}, gap {a - s[-1]}, "
              f"length difference {len(s) - 1}")

# Part C: Wilson primes and small K
for K in (1, 5, 13, 563):
    reps = canonical_reps(K)
    print(f"C  K = {K:3d}: canonical representations {reps}")
checks += 4
small = [K for K in range(1, 2001, 2) if K % 3 and not canonical_reps(K)]
print(f"C  odd K < 2001, K not divisible by 3, with no canonical representation: {len(small)} values; "
      f"first ones {small[:12]}")

# Part D: completeness cross-check against brute force with exponents <= EMAX_D
brute = {}
for m in range(0, EMAX_D):
    for s in combinations(range(2, EMAX_D + 1), m):
        brute.setdefault(K_of(s), []).append(s)
bound = 1 << (EMAX_D - 1)
nk = 0
for K in range(-1, bound, 2):
    a = sorted(canonical_reps(K), key=lambda s: (len(s), s))
    b = sorted(brute.get(K, []), key=lambda s: (len(s), s))
    assert a == b, (K, a, b)
    nk += 1
checks += nk
print(f"D  backward recursion = brute force (exponents <= {EMAX_D}) for all {nk} odd K in [-1, 2^{EMAX_D - 1}): OK")
print(f"total exact checks: {checks}")
