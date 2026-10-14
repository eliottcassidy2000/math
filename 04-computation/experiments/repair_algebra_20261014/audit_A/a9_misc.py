#!/usr/bin/env python3
"""Audit A: small residual checks for THM-4640.
(1) repair-term form of statement 3 on every canonical primitive pair with exponents <= E:
    last exponents e < e' = e + 2i  =>  K(s_-) - K(s'_-) = 2^e (4^i - 1)/3 (s = the one with the smaller last exponent);
(2) the 2r signed terms of a padded difference are pairwise distinct (no hidden cancellation), statement 6;
(3) reps of every primitive value with exponents <= E agree with brute force for E <= 16 (completeness of the
    pruned backward search used in a3/a7).
Usage: python3 a9_misc.py [E=20]"""
import sys, itertools
from collections import defaultdict
import numpy as np

E = int(sys.argv[1]) if len(sys.argv) > 1 else 20


def K_of(s):
    K = -1
    for e in s:
        K = 3 * K + (1 << e)
    return K


def canon_reps(K, emax):
    out = []

    def back(K, ehi, tail):
        if K == -1:
            out.append(tuple(reversed(tail)))
            return
        for e in range(2, ehi + 1):
            r = K - (1 << e)
            if r % 3:
                continue
            P = r // 3
            if P == -1 or 1 <= P <= 3 ** (e - 1) - 2 ** e:
                back(P, e - 1, tail + [e])

    back(K, emax, [])
    return out


S = np.array([-1], dtype=np.int64); P = np.array([], dtype=np.int64)
for e in range(1, E + 1):
    T = 3 * S + (1 << e)
    idx = np.searchsorted(S, T); idx[idx == len(S)] = len(S) - 1
    P = np.union1d(P, T[S[idx] == T]); S = np.union1d(S, T)
n1 = n2 = 0
for k in P.tolist():
    if k == -1:
        continue
    reps = canon_reps(k, E)
    for s, t in itertools.combinations(reps, 2):
        if s[-1] == t[-1]:
            continue
        if s[-1] > t[-1]:
            s, t = t, s
        e, i2 = s[-1], t[-1] - s[-1]
        assert i2 % 2 == 0
        i = i2 // 2
        assert K_of(s[:-1]) - K_of(t[:-1]) == (2 ** e) * (4 ** i - 1) // 3
        n1 += 1
        a, b = (s, t) if len(s) < len(t) else (t, s)
        L = len(b)
        ap = (1,) * (L - len(a)) + a
        terms = []
        for j, (x, y) in enumerate(zip(ap, b)):
            if x != y:
                terms += [3 ** (L - 1 - j) * 2 ** x, -(3 ** (L - 1 - j)) * 2 ** y]
        assert len(set(terms)) == len(terms) and len({abs(x) for x in terms}) == len(terms)
        n2 += 1
print(f"E <= {E}: repair-term identity on {n1} canonical primitive pairs OK; 2r terms pairwise distinct (also in "
      f"absolute value) on {n2} pairs OK")
# (3) completeness of the backward search against brute force, E <= 16
for EE in (12, 14, 16):
    brute = defaultdict(list)
    for n in range(EE):
        for s in itertools.combinations(range(2, EE + 1), n):
            brute[K_of(s)].append(s)
    bad = 0
    for k, rs in brute.items():
        if sorted(rs) != sorted(canon_reps(k, EE)):
            bad += 1
    print(f"backward search = brute force for every value of a canonical rep with exponents <= {EE}: "
          f"{len(brute)} values, {bad} mismatches")
