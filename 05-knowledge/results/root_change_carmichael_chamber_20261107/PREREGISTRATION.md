# Pre-registration: two root changes to 10^9, and the Carmichael chamber of root-change words

kind-pasteur-2026-11-07-S14. Reserves **THM-4904** and **HYP-9494**.

Owner prompt (2026-11-07): merge in Smith numbers and repunit primes in all bases, investigate prior work on the
Carmichael conjecture and the scaling of 14, 26, 34, 38, ..., while working on the two-root-change family.

Concurrent work, credited and not duplicated: Session218 (`05-knowledge/results/session218_smith_repunit_totient_20261107/`,
pushed 10:44 and 10:52 today) already covers the all-base Smith classification, the binary repunit multiplier
obstruction, the Mersenne-Wagstaff to nontotient bridge, the identification of 14, 26, 34, 38 (even nontotients versus
the repo's Camion exponents, which diverge at 50), dyadic repair of nontotients, and the THM-4884 multiplicity scope
repair. Session200 (`session200/carmichael-fibers/`) owns the source-profile Carmichael import (the O/J twin, singleton
profiles {3, 11} to 50,001); HYP-9477 (A) owns the strong-class exceptional set {3, 5, 11, 13}. This session works on
the root-change word grammar of HYP-9480 (D) (kind-pasteur S13), which none of these covers.

## Objects

Words over {O, Z, J} act on the rooted edge t = (D,E,U,d) = (1,0,2,1) by the physical moves
O(D,E,U,d) = (D,0,2U,d), Z = (D-E+U,d,2U-d,d), J = (D-E+U,D-E,D+U,D). The value of a word is the first coordinate
D of the final state. For K >= 0, M_K(N) is the number of words with at most K letters J, ending in Z or J (or empty),
whose value is N. This is the multiplicity counted by `rc_census.c` (the S13 enumerator `root_changes.c` with an
OpenMP task split and an 8-bit saturating counter). F_<=K is the set of values with M_K(N) > 0.

## Status of the run at registration

`rc_census.c` reproduced S13's `root_changes_1e8.txt` byte for byte on its coverage lines (X = 10^8, K = 2, 101 s on
12 threads). The run `rc_census 1000000000 2 512` was started before this file was written. **No line of its output
has been read.** This file is pushed before the output is opened.

Already seen (exploratory, X = 10^6 and 10^8 runs): M_K is even for almost all N; for K = 2 no N > 1 below 10^8 has
M_2(N) = 1; minima per dyadic window over N = 1 (mod 3) were below the minima over the other classes for windows
2^22 to 2^26.

## Predictions for X = 10^9, K = 2

- **P1 (coverage, HYP-9480 (D)).** The odd N <= 10^9 missed by F_<=2 are exactly the 140 misses below 10^8 (largest
  3,780,967). No odd N in (10^8, 10^9] is missed.
- **P2 (minimal multiplicity).** In each of the windows [2^26, 2^27), [2^27, 2^28), [2^28, 2^29) and [2^29, 10^9], the
  minimum of M_2(N) over odd N is at least 10.
- **P3 (no singleton words).** No odd N with 1 < N <= 10^9 has M_2(N) = 1.
- **P4 (mod-3 asymmetry).** In each of the windows [2^26, 2^27), [2^27, 2^28), [2^28, 2^29) and [2^29, 10^9], the minimum
  of M_2 over N = 1 (mod 3) is strictly smaller than both minima over N = 0 and N = 2 (mod 3).

Any failure is reported as a failure. P1 is the decisive test; P2-P4 are secondary.

## Planned proof content (not a prediction)

THM-4904 (to be proved and independently audited):

1. **Parity law.** M_K(N) = [N = 1] + 2 A_(K-1)(N) + B_K(N), where A_j(N) counts words w with at most j letters J
   and F(t.w) = N, and B_K(N) counts words with exactly K letters J. With Session199's adjoint F(t.w) = F(t.w^dagger)
   (reverse, swap O and Z, fix J), M_K(N) is congruent mod 2 to [N = 1] plus the number of dagger-fixed words with
   exactly K letters J and value N.
2. **Chamber.** Every word realizing a singleton value N > 1 (M_K(N) = 1) is dagger-fixed, has exactly K letters J,
   and neither begins nor ends with J (Session200's O/J twin).

Exact finite checks of (1) for K = 0, 1, 2 below a stated bound will be run by a separate verifier.

## Exploratory, not registered

Singleton values of M_0 and M_1, the comparison of chamber exceptional sets with Session200's profile set and
HYP-9477's class set, and scaling comparisons of the F_<=K miss lists with the even nontotients.
