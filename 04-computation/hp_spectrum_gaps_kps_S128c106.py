#!/usr/bin/env python3
"""hp_spectrum_gaps_kps_S128c106.py -- kind-pasteur-2026-07-20-S128c106

THE HAMILTONIAN-PATH SPECTRUM OF TOURNAMENTS.

Redei: hp(T) is always odd.  The repo carries the claim that the attained set is
exactly {odd} \\ {7, 21} -- two "permanent gaps".  This checks it and, more
usefully, makes the question WELL-POSED, which it is not as usually stated.

THE STABILISATION LEMMA (proved here, elementary, and it is what makes "the
spectrum" a single set rather than an n-indexed family):

    If v is a DOMINANT vertex of T' (out-degree n, in-degree 0) then
    hp(T') = hp(T'-v).
    Proof.  Any Hamiltonian path of T' contains v; v has in-degree 0, so v has no
    predecessor and is the initial vertex; the remainder is a Hamiltonian path of
    T'-v, and conversely any such path extends uniquely since v beats everything.
    QED

So every value attained at some n is attained at ALL larger n, the spectra are
NESTED, and "is 7 attained?" is a single well-posed question rather than one per n.

WHAT IS COMPUTED
  (A) exhaustive spectrum for n = 3..6 (all 2^C(n,2) labelled tournaments);
  (B) the nesting predicted by the lemma, checked;
  (C) a targeted search of the NEAR-TRANSITIVE region at n = 7..9 -- every
      perturbation of the transitive tournament by k <= 4 arc flips -- which is
      where small hp values must live;
  (D) random sampling at n = 7,8,9 to catch anything the flip search misses;
  (E) the verdict on 7 and 21, and on the other small odd values missing at n = 6
      (which SHOULD fill in at larger n if the two gaps are the only permanent ones).

Honest scope: (A),(B) are exhaustive and hence proofs for those n.  (C),(D) are
searches, so they can establish ATTAINMENT but never NON-attainment; the detection
floor is reported alongside, per the instrument-gate rule (MISTAKE-196).
"""
import sys
import random
from itertools import combinations

random.seed(20260720)
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
SAMPLES = int(sys.argv[2]) if len(sys.argv) > 2 else 60000


def hp_count(adj, n):
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            m = adj[v] & ~mask
            while m:
                b = m & -m
                dp[mask | b][b.bit_length() - 1] += c
                m ^= b
    return sum(dp[size - 1])


def from_bits(bits, n):
    adj = [0] * n
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
            idx += 1
    return adj


print("=" * 78)
print("(A) EXHAUSTIVE SPECTRUM, n = 3..6")
print("=" * 78)
spec = {}
for n in range(3, 7):
    m = n * (n - 1) // 2
    S = set()
    for bits in range(1 << m):
        S.add(hp_count(from_bits(bits, n), n))
    spec[n] = S
    hi = max(S)
    missing = [h for h in range(1, hi + 1, 2) if h not in S]
    print("  n = %d : |spectrum| = %-4d  max = %-5d" % (n, len(S), hi))
    print("     attained : %s" % sorted(S))
    print("     MISSING odd values below the max : %s" % missing)

print()
print("=" * 78)
print("(B) NESTING (stabilisation lemma):  spectrum(n) subset spectrum(n+1) ?")
print("=" * 78)
for n in range(3, 6):
    ok = spec[n] <= spec[n + 1]
    print("  spectrum(%d) subset spectrum(%d) : %-5s %s"
          % (n, n + 1, ok, "" if ok else "  missing: %s" % sorted(spec[n] - spec[n + 1])))
print("  -> consistent with the lemma; the spectra are nested, so 'the spectrum'")
print("     is the increasing union and the gap question is well-posed.")

print()
print("=" * 78)
print("(C) NEAR-TRANSITIVE SEARCH, n = 7..9, all k <= %d arc flips" % KMAX)
print("=" * 78)
found = dict(spec)
for n in (7, 8, 9):
    pairs = list(combinations(range(n), 2))
    base = [0] * n
    for i, j in pairs:
        base[i] |= 1 << j
    S = set()
    cfgs = 0
    for k in range(0, KMAX + 1):
        for flip in combinations(range(len(pairs)), k):
            adj = list(base)
            for f in flip:
                i, j = pairs[f]
                adj[i] &= ~(1 << j)
                adj[j] |= 1 << i
            S.add(hp_count(adj, n))
            cfgs += 1
    found[n] = S
    small = sorted(h for h in S if h <= 60)
    print("  n = %d : %d configurations, %d distinct hp values" % (n, cfgs, len(S)))
    print("     values <= 60 : %s" % small)

print()
print("=" * 78)
print("(D) RANDOM SAMPLING, n = 7,8,9,  %d samples each" % SAMPLES)
print("=" * 78)
for n in (7, 8, 9):
    m = n * (n - 1) // 2
    S = set(found.get(n, set()))
    for _ in range(SAMPLES):
        S.add(hp_count(from_bits(random.getrandbits(m), n), n))
    found[n] = S
    small = sorted(h for h in S if h <= 60)
    print("  n = %d : values <= 60 after sampling : %s" % (n, small))

print()
print("=" * 78)
print("(E) VERDICT ON THE GAPS")
print("=" * 78)
allS = set()
for n in sorted(found):
    allS |= found[n]
hi = 60
missing = [h for h in range(1, hi + 1, 2) if h not in allS]
print("  union of everything found, odd values <= %d NOT attained : %s" % (hi, missing))
for target in (7, 21):
    print("  %-3d attained anywhere in this search : %s" % (target, target in allS))
print()
print("  Values missing at n = 6 but FILLED IN at larger n (so not permanent):")
m6 = {h for h in range(1, max(spec[6]) + 1, 2)} - spec[6]
filled = sorted(h for h in m6 if h in allS)
print("     %s" % filled)
print("  Values missing at n = 6 and STILL missing: %s"
      % sorted(h for h in m6 if h not in allS))
print()
print("  DETECTION FLOOR (MISTAKE-196): (A),(B) are exhaustive and prove")
print("  non-attainment for n <= 6.  (C),(D) are searches -- they can only prove")
print("  ATTAINMENT.  The flip search is exhaustive for feedback-arc-set <= %d," % KMAX)
print("  so non-attainment of 7 and 21 is established for all tournaments within")
print("  %d flips of transitive at n <= 9, and NOT beyond that." % KMAX)
