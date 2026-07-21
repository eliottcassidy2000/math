#!/usr/bin/env python3
"""unstable_is_transitive_proof_kps_S128c132.py -- kind-pasteur-2026-07-20-S128c132

Toward a PROOF that the GIT-unstable tournaments (char_A with a root of multiplicity > n/2)
are EXACTLY the transitive ones (THM-1810, was THM-1805).  Two clean lemmas + the gap.

LEMMA A (mult > n/2 => integer eigenvalue).  If lambda is a root of char_A of multiplicity
  mu > n/2, its minimal polynomial f over Q is irreducible with f^mu | char_A, so
  deg(f)*mu <= n, forcing deg(f) = 1, so lambda in Q; being an eigenvalue of an integer
  matrix it is an algebraic integer, hence lambda in Z.  (Verified: in fact lambda = 0.)

LEMMA B (geometric multiplicity <= floor(n/2) off Perron & -1/2).  Since
  (A - lambda I) + (A - lambda I)^T = (J - I) - 2 lambda I = J - (1+2lambda) I, and
  rank(M) >= (1/2) rank(M + M^T), we get rank(A - lambda I) >= (1/2) rank(J - (1+2lambda)I).
  J - (1+2lambda)I has full rank n unless lambda = (n-1)/2 (Perron) or lambda = -1/2, so for
  every other lambda, rank(A - lambda I) >= n/2 and the geometric multiplicity
  g(lambda) = n - rank(A - lambda I) <= floor(n/2).  In particular rank(A) >= ceil(n/2).

CONSEQUENCE: mu > n/2 forces lambda in Z (Lemma A, so lambda != -1/2), and unless lambda is
the Perron value the geometric multiplicity is <= floor(n/2) < mu (Lemma B), so A is
NON-diagonalizable at lambda -- the excess multiplicity is Jordan structure.  Transitive
(lambda = 0, one nilpotent Jordan block of size n) realizes this.  The residual claim
"mu > n/2 => transitive" is thus reduced to the nilpotent/Jordan case.

This script VERIFIES both lemmas and the multiplicity GAP (algebraic mult lands in
{1..floor(n/2)} u {n}, never in between) exhaustively n<=6, and samples n=7.
"""
import sys
import numpy as np
import sympy as sp
from collections import Counter

x = sp.Symbol('x')


def from_bits(bits, n):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


print("=" * 86)
print("EXHAUSTIVE n<=6: Lemma A (mult>n/2 => integer), gap, rank(A)>=ceil(n/2)")
print("=" * 86)
for n in (3, 4, 5, 6):
    m = n * (n - 1) // 2
    lemA_ok = True
    gap_ok = True
    rank_ok = True
    unstable_lams = Counter()
    multset = set()
    minrank = n
    for bits in range(1 << m):
        A = from_bits(bits, n)
        M = sp.Matrix(A)
        p = sp.Poly(M.charpoly(x).as_expr(), x)
        fl = sp.factor_list(p.as_expr())
        # max root multiplicity and the responsible factor
        mu = max((e for _, e in fl[1]), default=1)
        multset.add(mu)
        rk = M.rank()
        minrank = min(minrank, rk)
        if rk < (n + 1) // 2:
            rank_ok = False
        if mu > n / 2:
            # find the factor(s) with that multiplicity; check its roots are integers
            for fac, e in fl[1]:
                if e == mu:
                    if sp.Poly(fac, x).degree() != 1:
                        lemA_ok = False
                    else:
                        root = sp.Poly(fac, x).all_roots()[0]
                        unstable_lams[sp.nsimplify(root)] += 1
                        if not sp.nsimplify(root).is_integer:
                            lemA_ok = False
            # gap: mu>n/2 should mean mu==n
            if mu != n:
                gap_ok = False
    print("  n=%d : max-mult set observed = %s   (gap (floor(n/2),n) empty: %s)"
          % (n, sorted(multset), gap_ok))
    print("        Lemma A (unstable eigenvalue is an integer): %s   unstable lambdas: %s"
          % (lemA_ok, dict(unstable_lams)))
    print("        min rank(A) over all T = %d   (>= ceil(n/2)=%d : %s)"
          % (minrank, (n + 1) // 2, rank_ok))
    sys.stdout.flush()

print()
print("=" * 86)
print("LEMMA B check (geometric mult <= floor(n/2) off Perron & -1/2), n=5,6 exhaustive")
print("=" * 86)
for n in (5, 6):
    m = n * (n - 1) // 2
    worst = 0
    bad = 0
    for bits in range(1 << m):
        A = np.array(from_bits(bits, n))
        ev = np.linalg.eigvals(A.astype(float))
        # cluster eigenvalues
        used = [False] * n
        for i in range(n):
            if used[i]:
                continue
            grp = [i]
            used[i] = True
            for j in range(i + 1, n):
                if not used[j] and abs(ev[i] - ev[j]) < 1e-6:
                    grp.append(j)
                    used[j] = True
            lam = ev[i]
            gmult = np.linalg.matrix_rank(A - lam * np.eye(n), tol=1e-6)
            gmult = n - gmult                     # geometric multiplicity
            perron = abs(lam - (n - 1) / 2) < 1e-6
            half = abs(lam + 0.5) < 1e-6
            if not perron and not half:
                worst = max(worst, gmult)
                if gmult > n // 2:
                    bad += 1
    print("  n=%d : max geometric mult off {Perron,-1/2} = %d  (<= floor(n/2)=%d : %s)"
          % (n, worst, n // 2, worst <= n // 2))
    sys.stdout.flush()

print()
print("=" * 86)
print("n=7 SAMPLE: is the gap (3.5, 7) empty? (no non-transitive with mult 4,5,6)")
print("=" * 86)
import random
random.seed(1)
n = 7
maxmult_nontrans = 0
worst_example = None
NS = 400000
for _ in range(NS):
    bits = random.getrandbits(n * (n - 1) // 2)
    A = np.array(from_bits(bits, n))
    ev = np.linalg.eigvals(A.astype(float))
    # multiplicity by clustering
    order = np.argsort(ev.real)
    evs = ev[order]
    best = 1
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(evs[j + 1] - evs[i]) < 1e-6:
            j += 1
        best = max(best, j - i + 1)
        i = j + 1
    c3 = int(round(np.trace(np.linalg.matrix_power(A, 3)) / 3))
    if c3 > 0 and best > maxmult_nontrans:      # non-transitive (has a 3-cycle)
        maxmult_nontrans = best
        worst_example = bits
# also the transitive one
Atr = np.array(from_bits(0, n))
evtr = np.linalg.eigvals(Atr.astype(float))
print("  transitive char_A max |eig| = %.3e (should be ~0, nilpotent) -> mult 7" % max(abs(evtr)))
print("  over %d random n=7 tournaments, max multiplicity among NON-transitive = %d"
      % (NS, maxmult_nontrans))
print("  floor(n/2) = 3, so gap forbids mult in {4,5,6}; observed max non-transitive mult %s 3"
      % ("<=" if maxmult_nontrans <= 3 else ">"))
print("  (SAMPLE, not exhaustive: detection floor is %d random tournaments.)" % NS)
