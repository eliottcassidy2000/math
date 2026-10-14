#!/usr/bin/env python3
"""Audit A (independent), THM-4640 (1) positivity/rigidity, (4) parity, (5c) padded differences, (7) gap census,
and a degeneracy test of the vanishing {2,3}-sums behind primitive collisions.

Representations: increasing exponent tuples s, K(()) = -1, K(t.e) = 3K(t) + 2^e; canonical = first exponent >= 2.
Usage: python3 a3_reps.py [E_GAPS=22] [E_PAD=14] [E_DEG=18] [E_LEMMA=20] [AMAX=40]
"""
import sys, time, itertools
from collections import defaultdict, Counter
import numpy as np

E_GAPS = int(sys.argv[1]) if len(sys.argv) > 1 else 22
E_PAD = int(sys.argv[2]) if len(sys.argv) > 2 else 14
E_DEG = int(sys.argv[3]) if len(sys.argv) > 3 else 18
E_LEMMA = int(sys.argv[4]) if len(sys.argv) > 4 else 20
AMAX = int(sys.argv[5]) if len(sys.argv) > 5 else 40
t0 = time.time()


def K_of(s):
    K = -1
    for e in s:
        K = 3 * K + (1 << e)
    return K


def canon_reps(K, emax):
    """all canonical increasing s with K(s) = K and exponents <= emax (exact; prunes with -1 or
    1 <= P <= max K over canonical exponents < e, which is K((2..e-1)) = 3^(e-1) - 2^e)"""
    out = []

    def back(K, ehi, tail):
        if K == -1:
            out.append(tuple(reversed(tail)))
            return  # a canonical non-empty prefix has K >= 1, so -1 can only be the empty prefix
        for e in range(2, ehi + 1):
            r = K - (1 << e)
            if r % 3:
                continue
            P = r // 3
            if P == -1 or 1 <= P <= 3 ** (e - 1) - 2 ** e:
                back(P, e - 1, tail + [e])

    if K == -1:
        return [()]
    back(K, emax, [])
    return out


# ---------------------------------------------------------------- positivity lemma (own enumeration, numpy)
def all_canonical(E):
    K = np.array([-1], dtype=np.int64); m = np.array([0], dtype=np.int64)
    first = np.array([0], dtype=np.int64); last = np.array([0], dtype=np.int64)
    for e in range(2, E + 1):
        K = np.concatenate([K, 3 * K + (1 << e)])
        first = np.concatenate([first, np.where(m == 0, e, first)])
        m = np.concatenate([m, m + 1])
        last = np.concatenate([last, np.full(len(last), e, dtype=np.int64)])
    return K, m, first, last


K, m, first, last = all_canonical(E_LEMMA)
ne = m > 0
lo = 3 ** (m[ne] + 1) - 2 ** (m[ne] + 2)
assert (K[ne] >= lo).all() and (lo >= 1).all()
assert ((1 << last[ne]) <= K[ne] + 3).all()
assert (((1 << last[ne]) == K[ne] + 3) == (m[ne] == 1)).all()
# parity / mod 8 (statement 4): K mod 8 = 3^(m-1) if e_1 = 2, -3^m if e_1 >= 3 (m = 0 counts as the second case)
pred = np.where(first == 2, pow3mod8 := np.where((m - 1) % 2 == 0, 1, 3), np.where(m % 2 == 0, 7, 5))
assert (K % 8 == pred).all()
# equal length -> equal representation (all increasing reps, canonical here; non-canonical = (1).t)
km = np.unique(np.stack([K, m]), axis=1)
assert km.shape[1] == len(K), "two distinct canonical reps of equal length share K"
print(f"lemma (K >= 3^(m+1) - 2^(m+2) >= 1, 2^e_m <= K+3, = iff m = 1), K mod 8 rule, and 'equal length => equal "
      f"representation' hold on all {len(K)} canonical s with exponents <= {E_LEMMA}  [{time.time()-t0:.1f}s]", flush=True)
del K, m, first, last, km, lo, pred

# ---------------------------------------------------------------- anchor rigidity
sec = []
for a in range(2, AMAX + 1):
    reps = canon_reps((1 << a) - 3, a)
    assert (a,) in reps
    others = sorted(r for r in reps if r != (a,))
    if others:
        sec.append((a, others))
print(f"anchors 2^a - 3, 2 <= a <= {AMAX}: second canonical representation for a in {[a for a, _ in sec]}; "
      f"counts {[len(o) for _, o in sec]}  [{time.time()-t0:.1f}s]", flush=True)
for a, o in sec:
    print(f"   a={a}: {o}")
# smallest value with two canonical representations: all reps of K <= 2^12 - 3 have exponents <= 12
two = [k for k in range(-1, 4093, 2) if len(canon_reps(k, 12)) >= 2]
print(f"smallest values with >= 2 canonical representations: {two[:6]}", flush=True)

# ---------------------------------------------------------------- primitive values and their canonical reps (E_GAPS)
# primitive K at level E: from the set recursion (same as a2_census.py)
def primitive_values(E):
    S = np.array([-1], dtype=np.int64); P = np.array([], dtype=np.int64)
    for e in range(1, E + 1):
        T = 3 * S + (1 << e)
        idx = np.searchsorted(S, T); idx[idx == len(S)] = len(S) - 1
        P = np.union1d(P, T[S[idx] == T])
        S = np.union1d(S, T)
    return P


Pv = primitive_values(E_GAPS)
nrep = Counter(); gaps = Counter(); lend = Counter(); padhist_c = Counter(); npairs = 0
pairs_by_E = []
for k in Pv.tolist():
    if k == -1:
        continue
    reps = canon_reps(k, E_GAPS)
    nrep[len(reps)] += 1
    for s, t in itertools.combinations(reps, 2):
        if s[-1] == t[-1]:
            continue
        assert len(s) != len(t) and (len(s) - len(t)) % 2 == 0 and (s[-1] - t[-1]) % 2 == 0
        gaps[abs(s[-1] - t[-1])] += 1; lend[abs(len(s) - len(t))] += 1; npairs += 1
        if len(s) > len(t):
            s, t = t, s
        sp = (1,) * (len(t) - len(s)) + s
        r = sum(1 for x, y in zip(sp, t) if x != y)
        assert r >= 3
        pairs_by_E.append((max(s[-1], t[-1]), k, s, t, r))
print(f"E <= {E_GAPS}: {len(Pv) - 1} non-root primitive values; canonical reps per value {dict(sorted(nrep.items()))}; "
      f"{npairs} primitive canonical pairs; gaps {dict(sorted(gaps.items()))}; length differences "
      f"{dict(sorted(lend.items()))}; min padded difference >= 3 for all  [{time.time()-t0:.1f}s]", flush=True)

# padded-difference histogram at E_PAD (canonical pairs; x4 for twin variants, +1 root)
hc = Counter(r for (E, k, s, t, r) in pairs_by_E if E <= E_PAD)
print(f"E <= {E_PAD}: {sum(hc.values())} canonical primitive pairs, padded differences {dict(sorted(hc.items()))}; "
      f"with twin variants x4 + root: {4 * sum(hc.values()) + 1}", flush=True)

# ---------------------------------------------------------------- degeneracy of the vanishing sums
def has_vanishing_proper_subsum(terms):
    n = len(terms)
    h = n // 2
    A, B = terms[:h], terms[h:]
    sa = defaultdict(list)
    for mask in range(1 << len(A)):
        sa[sum(A[i] for i in range(len(A)) if mask >> i & 1)].append(mask)
    full = (1 << n) - 1
    for mb in range(1 << len(B)):
        sb = sum(B[i] for i in range(len(B)) if mb >> i & 1)
        for ma in sa.get(-sb, ()):
            tot = ma | (mb << h)
            if tot != 0 and tot != full:
                return True
    return False


deg = Counter(); ndeg_ex = []
for (E, k, s, t, r) in pairs_by_E:
    if E > E_DEG:
        continue
    sp = (1,) * (len(t) - len(s)) + s
    L = len(t)
    terms = []
    for j, (x, y) in enumerate(zip(sp, t)):
        if x != y:
            terms.append(3 ** (L - 1 - j) * 2 ** x); terms.append(-(3 ** (L - 1 - j)) * 2 ** y)
    assert sum(terms) == 0
    d = has_vanishing_proper_subsum(terms)
    deg[d] += 1
    if d and len(ndeg_ex) < 5:
        ndeg_ex.append((k, s, t, r))
print(f"E <= {E_DEG}: vanishing 2r-term sums of primitive canonical pairs: degenerate (a vanishing proper subsum) "
      f"{deg[True]}, non-degenerate {deg[False]}; examples {ndeg_ex}  [{time.time()-t0:.1f}s]", flush=True)
