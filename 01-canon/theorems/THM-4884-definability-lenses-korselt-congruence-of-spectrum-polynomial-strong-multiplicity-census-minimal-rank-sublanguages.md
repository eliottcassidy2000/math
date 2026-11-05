---
id: THM-4884
title: "Definability lenses on the Hamiltonian-path spectrum. (1) Korselt congruence (PROVED): the spectrum polynomial P_n(x) = sum over labelled tournaments of x^H satisfies P_n = F_{n,p} (mod p) for every prime p | n, F_{n,p} the sum over tournaments invariant under the rotation by n/p; for prime p the invariant tournaments are the circulants and their counts are multiples of p, so the number of labelled p-tournaments with a count not divisible by p is divisible by p (checked exactly for n <= 7: at n = 7 only h = 175 and 189 escape). (2) Carmichael multiplicity (FINITE-EXACT from THM-4744's inventories): m(h) = number of strong classes with H = h satisfies m = 0 only at 7, 21, m = 1 only at 3, 5, 11, 13, m = 2 only at 9, 17, 19, 35, 39, 47, 63, m = 3 only at 25, m >= 6 for every attained h in (63, 625), and the multiplicities 9, 13, 21, 22, 34, 37 do not occur below 625; Session 200's profile singletons {3, 11} differ from the class singletons exactly by the collision at 5 and the absence of a uniquely Hamiltonian realization of 13, and the C1 readout misses exactly Session 197's 29 values through 50,001. (3) Minimal rank (FINITE-EXACT): the sub-languages of G(0,2) without double events (rank 14), with first-in-first-out closing (rank 14) and with last-in-first-out closing (rank 10) are not complete to 2047 (misses 43, 47, 63, 79, 103, 113, 235, 709; 24 values; 582 values), so shared-endpoint configurations are indispensable and the smallest complete candidate found keeps rank 17"
status: >
  PROVED: (1) (orbit counting; the free action on Hamiltonian paths). FINITE-EXACT: (2) read off THM-4744's complete
  inventories (strong_ears_375.out, strong_ears_625.out); (3) the automata, ranks (modular, two primes; empty
  tournament = one path) and censuses to 2047 (S79 engine with a letter filter). KNOWN / CITED: the orbit-counting
  mechanism of THM-4746 (Camion pseudoprimes); converse pairing P_d = SC_d (mod 2) (THM-4744 (3)); the grammars,
  ranks and censuses of THM-4878; Session 200's Diophantine source compiler (concurrent, same directive); the
  attached "Z is Diophantine in Q" paper is used only as an architecture and is not verified.
  FINITE-EXACT (2b, merging Session 200's carmichael-fibers packet): its profile inventory replayed (5,757 / 126,910 /
  3,816,297 profiles at U <= 624 / 5,000 / 50,000, agreeing with that lane and its auditor); profile singletons {3, 11}
  versus class singletons {3, 5, 11, 13}; the C1 readout misses exactly Session 197's 29 values through 50,001 (the
  list through 5,001 is Session 197's C1_COMPLETE.md, KNOWN; the range (5001, 50001] is new).
  Pre-registered Q1 (confirmed in its containment, size and membership clauses; "41" and "every multiplicity 1..10"
  refuted), Q2 (completeness of modes 1 and 2 refuted; mode 3 and the rank ordering mode 3 < mode 1 confirmed; the
  ordering mode 2 < mode 1 refuted), Q3 (a check).
  NOT CLAIMED: Chan-Pak coverage; any statement about the paper's correctness.
  AUDIT: PENDING.
session: opus-2026-11-05-S81
source: 05-knowledge/results/definability_lenses_spectrum_diophantine_carmichael_20261105.md
related:
  - THM-4878
  - THM-4882
  - THM-4744
  - THM-4746
  - THM-4747
  - THM-4748
  - THM-4719
  - HYP-9470
  - HYP-9477
  - session200/diophantine-source-compiler/PROOF.md
  - session200/height-bounded-source-lifting/PROOF.md
  - session200/sharp-dyadic-fibers/PROOF.md
  - session200/carmichael-fibers/PROOF.md
  - session197/lowcycle/C1_COMPLETE.md
---

# THM-4884 - definability lenses: Korselt congruence, Carmichael multiplicity, minimal rank

## Statements

1. **Korselt congruence (PROVED).** Let P_n(x) = sum over the labelled tournaments on Z_n of x^(H(T)). For every
   prime p dividing n, P_n(x) = F_{n,p}(x) modulo p coefficientwise, where F_{n,p} sums over the tournaments invariant
   under i -> i + n/p. *Proof.* The rotation by n/p generates a group of order p acting on labelled tournaments by
   relabelling, H is invariant, orbits have size 1 or p. For n = p prime the invariant tournaments are the circulants
   (arc i -> j iff j - i lies in a half-system S of Z_p^*), and the rotation acts freely on their Hamiltonian paths (the
   first vertex moves), so each circulant count is a multiple of p. Hence: for prime p and p not dividing h, p divides
   N(p, h), the number of labelled p-tournaments with exactly h Hamiltonian paths. Exact check for n <= 7 in
   `spectrum_polynomial_congruence_20261105.cpp` (at n = 7 the residues are nonzero only at h = 175 (6 mod 7) and
   h = 189 (2 mod 7), the eight circulant counts).

2. **Carmichael multiplicity (FINITE-EXACT).** With m(h) the number of strong tournament classes of all orders
   having H = h (finite by Moon's floor), THM-4744's inventories give for odd h < 625: m(h) = 0 iff h in {7, 21};
   m(h) = 1 iff h in {3, 5, 11, 13}, each realised by one self-converse class (forced by converse pairing for an odd
   multiplicity); m(h) = 2 iff h in {9, 17, 19, 35, 39, 47, 63}; m(h) = 3 iff h = 25; m(h) >= 6 for every attained h
   with 63 < h < 625 (minimum 6 at h = 67); m(h) >= 73 for odd h in [375, 623]. The values 9, 13, 21, 22, 34, 37 are
   not multiplicities below 625, while every
   other value up to 38 is. So the analogue of Carmichael's conjecture fails exactly at the four smallest primes other
   than 7, and the analogue of Ford's theorem (every multiplicity >= 2 occurs) fails.

   *Profile fibers versus class multiplicity (FINITE-EXACT, merging Session 200).* Let pf(n) be Session 200's number
   of numerical source profiles (D, E, U, d) of the physical O/Z/J grammar from (1, 0, 2, 1) with C1 readout
   D - E + U = n. Replaying that inventory (`profile_vs_class_multiplicity_20261105.py`; profile counts agree with
   Session 200 and its auditor at U <= 5,000 and U <= 50,000): pf(h) = 1 and m(h) >= 2 never happens below 625;
   m(h) = 1 and pf(h) >= 2 only at h = 5 (two profiles, one class); 13 has m = 1 and pf = 0. Through 50,001 the odd
   n with pf(n) = 0 are exactly Session 197's 29 values 7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247,
   251, 277, 299, 421, 449, 503, 533, 569, 583, 637, 821, 1067, 1339, 1861, 2647 (KNOWN through 5,001), so no odd
   value in (2647, 50001] is missed by the C1 readout. The two multiplicities are incomparable (pf > m at 5, 9, 17,
   19, 27, 35, 39, 53, 61, 63; pf < m at 23, 25, 29, 33, 37, 41, 43, 45, 49, 51, 55, 57, 59 below 64).

3. **Minimal rank of a complete sub-grammar (FINITE-EXACT).** Sub-languages of the cutwidth-2 grammar G(0,2)
   (THM-4878; rank 17 with the empty tournament counting one path): no simultaneous double events (every chain vertex
   opens at most one reversed arc and terminates at most one): 64 cut states, Hankel rank 14, seam rank 10 at two
   pending arcs, attains 1014 of the 1024 odd numbers <= 2047 and misses 7, 21, 43, 47, 63, 79, 103, 113, 235, 709;
   crossing chains (additionally first-in-first-out closing): rank 14, 998 of 1024; nested chains (last-in-first-out):
   58 states, rank 10, 440 of 1024. Consequently the shared-endpoint configurations (a vertex beaten by two later
   vertices or beating two earlier ones) are indispensable for completeness through cutwidth 2 (S79's exact
   minimal-arc witnesses of 43, 47, 63, 79, 103, 113, 235 each contain a shared tail or a shared head), and the
   smallest complete candidate found keeps rank 17.

4. **Architecture map (reasoning, no claim).** The attached paper's assembly (indices on a rank-one group, guarded
   multiplication, height estimate, local conditions) corresponds to the one-apex orbit, ordered joins, Moon's floor
   with the counting lemma, and modular synchronization; the missing ingredient for the positive spectrum is a group
   law on counts, which THM-4882 (4) excludes for repetition and which the signed theorem of Session 198 supplies only
   for signed responses. The source note, sections 1, 2, 6 and 7, records this and the three candidate numeration
   steps.

## Programs

`04-computation/experiments/spectrum_polynomial_congruence_20261105.cpp`, `spine_grammar_sublanguages_20261105.cpp`
(S79's engine with the SPINE_MODE letter filter), `spine_grammar_ranks_mod_20261104.py --empty1`,
`definability_lenses_prereg_20261105.md`, `profile_vs_class_multiplicity_20261105.py` (Session 200's physical
inventory replayed, pf(h) beside m(h), the C1 missed list through 50,001); the multiplicity census is read from
`04-computation/experiments/h_spectrum_multiplicities_20261021/strong_ears_375.out` and `strong_ears_625.out`.
