---
id: THM-4985
title: "Generating every tournament from first creators beyond order 11. (1) FINITE-EXACT: THM-4942's five-operation grammar (first creator, first creator plus a vertex, ordinal sum, substitution, ordinal sum re-connected by at most two arcs), exhaustive through order 11, fails at order 12: one-vertex extensions of S17's seven 3-arc-strong order-11 non-creators give 12 escaper classes (6 converse pairs, H = 48,997 .. 52,303), prime, 3-arc-strong, with every vertex deletion a non-creator. (2) FINITE-EXACT, found by annealing and verified exactly: 181 further classes (193 in all; H from 32,897 to 58,017; 17 self-converse), including two self-converse classes of decomposability defect 3 (H = 48,537 and 48,849), which are not within two arc reversals of any non-strong or non-prime tournament; every escaper found has lambda = 3 and creator depth 2, except one of depth 3 (H = 33,021). (3) FINITE-EXACT, witness-certified: the band tournaments band_3(13) and band_4(16) are escapers with lambda = 3 and 4. (4) PROVED: band_s(n) is prime for n >= 2s + 3, has lambda = s for n >= 2s + 1, and every tournament whose back arcs have span <= s has at most (s + 3)^(2n) Hamiltonian paths. (5) PROVED, conditional on a superexponential solid prefix (implied by HYP-9481 (a)): for all k and r some tournament escapes the grammar with 'plus at most k vertices' and 're-connected by at most r arcs', and k*(n) -> infinity. (6) FINITE-EXACT: first creators are 92.0% of the classes of order 11; PROVED conditional on max H(m) = o(m·m!/2^(m-1)): almost every tournament is a first creator, and every tournament is an induced subtournament of one"
status: >
  PROVED + FINITE-EXACT + PROVED-CONDITIONAL.
  - (1) FINITE-EXACT: exhaustive over the 7 × 2^11 one-vertex extensions of the seven 3-arc-strong order-11
    non-creators (S17's list). Exact H for T and every deletion; the exact spectra S(10), S(11) (klein); lambda by cut
    enumeration; primality by module enumeration; canonical forms. Re-verified by independent pure-Python code
    (verify_escapers.py). With THM-4942 (3), order 12 is the first failure. HYP-9532 (F) is answered: NO.
  - (2) FINITE-EXACT examples. The search is EXPLORATORY, so no completeness is claimed; every reported class is
    re-verified exactly (analyze_escapers.py).
  - (3) FINITE-EXACT: memberships in S(12), S(14), S(15) are certified by explicit witness tournaments; every witness H
    is recomputed in pure Python (verify_bands.py).
  - (4) PROVED. Lemma B2 uses THM-4811 (B) for n <= 4s + 1 and a coprime bijection {1..s} -> {s + 1..2s} for n >= 3s
    (Daykin-Baines 1963; in general Pomerance-Selfridge 1980; CITED). The elementary choice sigma(r) = rs + 1 covers
    n >= s^2 + s + 1, which is all that (5) needs.
  - (5) PROVED, conditional on (S-inf): S(m) contains every odd h <= g(m) other than 7 and 21, with
    g(m)^(1/m) -> infinity.
  - (6) FINITE-EXACT shares from klein's census. The density statement is PROVED conditional on (M_eps): max H(m) <=
    (1 - eps)((m + 1)/2) m!/2^(m - 1) for large m. (M_eps) is OPEN; Alon's 1990 upper bound carries m^(3/2).
  KNOWN and credited:
  - ADDED 2026-11-09 after audit P1 of mac-mini-2026-11-09-paleyeven (MISTAKE-746): the "decomposability defect" dec
    (= the module defect for n >= 3) is the classical SIMPLICITY INDEX s(T). See J. Müller, J. Pelant, Czechoslovak
    Math. J. 24 (1974) 378-391 (s(T) <= (n - 1)/2, with equality iff T is doubly regular) and H. Boussaïri,
    A. Lakhlifi, M. Talbaoui, Contrib. Discrete Math. 17 (2022) 172-186, arXiv:1907.11777 (s(T) <=
    min(delta0, twin distance)). The statements about escapers are unaffected. THM-4995 continues this lane
    (escapers of simplicity index 4 at order 16 and 5 at order 18; unbounded index under (S_c));
  - THM-4942 / HYP-9532 (kind-pasteur S17): the grammar, the seven order-11 parents, the connectivity floors;
  - THM-4922 (S16): first creators;
  - klein-2026-11-06: the exact spectra;
  - THM-4811 (B) and Session230 §5: lambda = minimum semidegree for n <= 4k + 1;
  - HYP-2879, THM-4097, THM-4114: the attachment form behind wit.c and ext13.c;
  - THM-1370: the insertion lemma and the 7/21 omission;
  - THM-589 (both legacy files): Var H;
  - Szele (1943), Alon (1990), Friedgut-Kahn, Adler-Alon-Ross: the bounds on max H;
  - THM-4880 (mac-mini-2026-11-04-ribbon): the band tournaments ARE its ribbon tournaments R(a^n) with constant letter
    a = {1..s}. Its §F already lists H(band_2(n)) and H(band_3(n)) for n <= 10. Lemma B3 is its Theorem A window
    argument made explicit. New here: primality and lambda (B1, B2), the escaper certificates, and (5). THM-4952 gives
    the packing-polynomial transfer. Credit added after audit G1 (M2).
  - Session230 (primes/REPORT.md §5, SYNTHESIS §6) posed the bounded-defect question that (2) and HYP-9565 (A) make
    precise. Credit added after audit G1 (M3).
  - THM-329 (opus-2026-05-27-S6): the search values 531,205 (order 12) and 3,719,831 (order 13), used for R_H(13), R_H(14)
    in THM-4986 (5).
  NOT CLAIMED:
  - the complete list of order-12 escapers;
  - k*(12) = 3 exactly (HYP-9565 (C));
  - any statement about Chan-Pak completeness.
  AUDIT: independent audit G1 (05-knowledge/results/genesis_20261109/AUDIT_G1.md), PASS WITH CORRECTIONS: 3 MAJOR
    (an unsourced max H(12) value; THM-4880 and Session230 priors uncredited) and 23 MINOR, all applied; MISTAKE-743.
    Every computational claim reproduced with independent code.
session: mac-mini-2026-11-09-genesis
date: 2026-11-09
source: 05-knowledge/results/genesis_20261109/PROOF.md
related:
  - THM-4942
  - HYP-9532
  - THM-4922
  - HYP-9512
  - THM-4930
  - THM-4811
  - THM-1370
  - THM-589
  - THM-4744
  - THM-4880
  - THM-4952
  - HYP-9481
  - THM-4986
  - HYP-9565
  - 05-knowledge/results/session230_creator_grammar_20261109/README.md
  - 05-knowledge/results/path_spectrum_two_angles_metagraph_20261106/README.md
script: 04-computation/experiments/genesis_20261109/ext12.py; 04-computation/experiments/genesis_20261109/tprop.c; 04-computation/experiments/genesis_20261109/canon.c; 04-computation/experiments/genesis_20261109/sa12.c; 04-computation/experiments/genesis_20261109/analyze_escapers.py; 04-computation/experiments/genesis_20261109/verify_escapers.py; 04-computation/experiments/genesis_20261109/wit.c; 04-computation/experiments/genesis_20261109/verify_bands.py; 04-computation/experiments/genesis_20261109/band_lemmas.py; 04-computation/experiments/genesis_20261109/fc_share.py
---

# THM-4985 - generation from first creators fails at order 12; band escapers; first creators are typical

Proofs and data: [PROOF.md](../../05-knowledge/results/genesis_20261109/PROOF.md), and the
[note](../../05-knowledge/results/genesis_20261109/README.md).

**Definitions.**
- A **first creator** of order n is a class with H not in S(n - 1).
- An **escaper** is a strong, prime tournament with lambda >= 3 that is not a first creator and has no first-creator
  deletion. Escapers are exactly the tournaments that THM-4942's five-operation grammar does not cover: the lambda >= 3
  part of S17's REST category.
- The **decomposability defect** is the fewest arc reversals reaching a non-strong or non-prime tournament. For n >= 3 it
  equals the module defect, since a non-strong tournament is never prime (audit G1, m1).
- **band_s(n)** has vertices 0..n - 1. For i < j it has j -> i iff j - i <= s, and i -> j otherwise. It is THM-4880's
  ribbon tournament with the constant letter {1..s}.

1. **(FINITE-EXACT) The grammar fails at order 12.** Extend the seven 3-arc-strong order-11 non-creators of THM-4942
   in all 2^11 ways.
   - Exactly 12 escaper classes appear: 6 converse pairs, all rigid.
   - The H values are 48,997, 49,215, 49,935, 50,221, 51,451 and 52,303, one converse pair each.
   - All 12 have lambda = 3 and are one arc reversal from a substitution.
   - The grammar holds through order 11 (THM-4942 (3)), so order 12 is its first failure.
   - k*(12) >= 3.
2. **(FINITE-EXACT examples, EXPLORATORY search) More escapers at order 12.** Annealing finds 181 further classes, none among the 12 of (1), so 193 in all.
   - H runs from 32,897 to 58,017; 17 classes are self-converse.
   - 190 classes are one reversal from a substitution and one class is two reversals away.
   - **Two self-converse classes, with H = 48,537 and 48,849, have decomposability defect 3.**
   - Every escaper found has lambda = 3. Every one contains a first creator on 10 vertices (creator depth 2) except one,
     with H = 33,021, of depth 3.
3. **(FINITE-EXACT, witness-certified) Bands.**
   - band_3(13) has H = 237,817, which is in S(12), and every deletion has H <= 72,073.
   - band_4(16) has H = 64,053,917, which is in S(15), and every deletion has H in S(14).
   - Both are escapers, with lambda = 3 and 4 respectively. So the grammar with "ordinal sum re-connected by at most 3
     arcs" still fails at order 16.
   - band_4(15) is a non-creator with lambda = 4 (H = 15,927,743, in S(14) by the same witness file), so
     k*(15) >= 4 and k*(16) >= 4.
4. **(PROVED) Band lemmas.**
   - band_s(n) is prime for n >= 2s + 3, and not prime for n = 2s + 2.
   - lambda(band_s(n)) = s for n >= 2s + 1.
   - If all back arcs of T (with respect to some vertex order) have span <= s, then H(T) <= (s + 3)^(2n).
5. **(PROVED, conditional on (S-inf)) No bounded grammar of this kind.** For every s >= 3 and k >= 0, and all large n,
   band_s(n):
   - is not a first creator;
   - is not a first creator plus at most k vertices;
   - is strong and prime, with lambda = s.

   So no grammar {FC; FC plus at most k vertices; ordinal sum; substitution; ordinal sum re-connected by at most r arcs}
   with fixed (k, r) generates every tournament, and HYP-9532 (D)'s k*(n) is unbounded.
6. **First creators are typical.**
   - (FINITE-EXACT) The share of first-creator classes is 0.719, 0.793, 0.8545, 0.8884 and 0.9200 at orders 7-11. Order
     11 is new: 831,487,872 of 903,753,248 classes.
   - (FINITE-EXACT) max H(n - 1)/E[H(n)] is 0.40 at n = 11.
   - (PROVED, conditional on (M_eps)) A random labelled n-tournament is a first creator with probability
     >= 1 - (W(n)/n! - 1)/eps^2 = 1 - O(1/n).
   - (PROVED, conditional on (M_eps)) Every tournament is an induced subtournament of a first creator.
