---
id: THM-5071
title: "Weighted cycle formula for the Hamiltonian path total of every digraph. (1) For a loopless digraph D, H(D) = sum over permutations whose nontrivial cycles are directed cycles of D (sign +1) or of its complement (sign (-1)^(length-1)) - Grinberg-Stanley's p-expansion of the Redei-Berge function at (1,0,0,...). (2) Equivalently H(D) = sum over packings of vertex-disjoint admissible cycles of prod w(C): w = 2 for an odd cycle through tournament pairs only, 0 for an even one, 1 for any cycle of D through a reciprocal (11) pair, (-1)^(length-1) for any cycle of the complement through a missing (00) pair; loops never enter. (3) Parity: H(D) = the number of packings of defect cycles (mod 2), refining Berge's parity theorem (Redei is the tournament case). (4) Every positive integer is the total of a semicomplete digraph (TT_n plus sink back arcs at the binary digits of h - 1); with exactly one reciprocal pair and no missing pair the total is never 1, 4 or 7, at any order (1 and 7 KNOWN, Session234; 4 new); with exactly one missing pair every total 0..129 occurs at orders <= 7 (7 at order 5, already in Session234's census; 21 at order 6). (5) Every strong semicomplete core with H = 7 has no odd tournament cycle: 7 = 1 + 6 x 1 (two classes) or 1 + 5 x 1 + 1 x 1 (one class)."
status: >
  KNOWN / CITED: (1) is Grinberg-Stanley, arXiv:2307.05569v1, Theorem 1.30 (Section 1.8.2), specialised at
  (1, 0, 0, ...), where F_S(1, 0, ...) = [S empty]; it is restated as Irving-Omar, arXiv:2412.10572v3, Corollary 14.
  The sign convention (phi sums (length - 1) over the D-cycles of sigma) agrees with theirs. Checked here on 232,263
  digraphs (pre-registered gate P5) and by the audit on 1,372,741 more (0 failures).
  PROVED:
  - (2), a regrouping of (1). Pairing a D-cycle with its reversed complement cycle is GS's own device (their proof of
    the tournament p-expansion, Section 3).
  - (3), a refinement of Berge's theorem ham(D) = ham(Dbar) (mod 2). The defect cycles of D and of Dbar are the same
    set, so (3) implies Berge. Whether GS Section 6 already states (3) was not checked, so no novelty is claimed for
    (3).
  - (4a).
  - (4b), by component factorisation (Session234 DG3), a strong spanning resolution (Session234 DG4; the four-cell
    lemma of PROOF.md is its b = 1 case), Moon's nondecreasing strong floor with f(5) = 9 (THM-1370 (3), THM-4744),
    and an exhaustive check of orders <= 4. The cases 1 and 7 are KNOWN: Session234 DG3/DG4 and its DG5 refinement
    ("H = 7 forces exactly two reciprocal pairs"). New: the exclusion of 4, by the same route.
  FINITE-EXACT: (4c) at orders <= 7 (pre-registered P4 CONFIRMED). The H = 7 one-missing-pair witness at order 5 was
  already in Session234's census (general_census.csv, row 5,7,1,0,1,600); the order-6 H = 21 witness and the 0..129
  census are new.
  (5): the classes are enumerated at orders <= 5. The all-order form holds because every strong H = 7 core has order 4
  (Session234 DG4/DG5).
  AUDIT: PASS WITH CORRECTIONS (05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/audit/AUDIT.md;
  M1 credit for (4b), m1, m2, m8, m9, m12 applied; MISTAKE-748).
  KNOWN and credited:
  - Session234 refuted the owner's conjecture first and classified the semicomplete 7/21 families (README, PROOFS);
    the core counts in (5) agree with it, as does 21 with one reciprocal pair at order 6. Its DG1 notes that loops
    never enter.
  - OCF and Redei are the tournament cases.
session: kind-pasteur-2026-11-09-S18
source: 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/README.md
proof: 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/PROOF.md
related:
  - THM-1370 (h-spectrum-omits-7-21-all-n)
  - THM-4744
  - HYP-9640
  - THM-5070
  - 05-knowledge/results/session234_digraphs_20261109/README.md
verification:
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/digraph_weights.c
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/core_weights.py
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/one_defect_witness.py
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/out/digraph_weights.out
  - 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/out/core_weights.out
---

# THM-5071 - the weighted cycle formula for digraph path totals

## Setting

D is a digraph on V, written as a square binary matrix. Loops are allowed, but no Hamiltonian path uses them, so
remove them. Each unordered pair {i, j} is:
- a **T-pair** if exactly one of A_ij, A_ji is 1;
- a **B-pair** (reciprocal, "11") if both are 1;
- a **Z-pair** (missing, "00") if both are 0.

Write Dbar for the loopless complement. A tournament has only T-pairs; a semicomplete digraph has no Z-pair; an
oriented graph has no B-pair.

## Statement

**(1)** H(D) = sum over sigma in Sym(V), each nontrivial cycle of sigma a directed cycle of D or of Dbar, of
prod over cycles [+1 if the cycle is in D; (-1)^(len - 1) if it is in Dbar].

**(2) Weighted packing form.**

    H(D) = sum over packings of vertex-disjoint admissible cycles of prod w(C),

with:

| admissible cycle C | w(C) |
|---|---|
| odd cycle of D using only T-pairs | 2 |
| even cycle of D using only T-pairs | 0 |
| cycle of D (length >= 2) using at least one B-pair | 1 |
| cycle of Dbar (length >= 2) using at least one Z-pair | (-1)^(len - 1) |

- For tournaments this is OCF, H = I(Omega, 2).
- For semicomplete D every weight is 0, 1 or 2.

**(3) Parity.** H(D) is congruent, mod 2, to the number of packings (the empty one included) of cycles of the last
two kinds. This refines Berge's theorem H(D) = H(Dbar) (mod 2); Redei's theorem is the tournament case.

**(4) Spectra with few defects.**

- **(4a) Sink back arcs.** Let D be TT_n (i -> j for i < j) plus arcs from the sink n-1 to each x in S. Then
  H(D) = 1 + sum_{x in S} 2^(n-2-x). So every positive integer h is the total of a semicomplete digraph with
  popcount(h - 1) reciprocal pairs.
- **(4b) One reciprocal pair.** If D is loopless, semicomplete and has exactly one B-pair, then H(D) is not 1, 4 or 7.
  This holds at every order.
- **(4c) Census, orders <= 7.**
  - With exactly one B-pair (and no Z-pair), the totals below 198 are exactly {2, ..., 197} minus {4, 7}.
  - With exactly one Z-pair (and no B-pair), every total 0..129 occurs, including 7 (order 5) and 21 (order 6).

**(5) Cores.** The strong semicomplete digraphs with H = 7 form 3 classes, all of order 4 with b = 2. None has an odd
T-only cycle. The enumeration covers orders <= 5. It is complete at all orders because every strong H = 7 core has
order 4 (Session234 DG4/DG5).

| classes | decomposition |
|---|---|
| two | 7 = 1 + 6 x 1: six pairwise intersecting B-cycles |
| one | 7 = 1 + 5 x 1 + 1 x 1 |

The 16 classes with H = 21 at order 5 decompose with at most three odd T-only cycles; see `out/core_weights.out`.

## Proof

**(1).** Grinberg and Stanley (arXiv:2307.05569v1, Theorem 1.30) expand the Rédei-Berge function as
U_X = sum_sigma (-1)^phi(sigma) p_type(sigma). Here sigma runs over the permutations whose nontrivial cycles are
cycles of X or of Xbar, and phi(sigma) sums (length - 1) over the cycles in X.
- The signs follow from X = 2-cycle (p_1^2 - p_2), X = empty pair (p_1^2 + p_2) and X = triangle (p_1^3 + 2 p_3).
- At (1, 0, 0, ...) the left side counts listings with no consecutive X-arc, which are the Hamiltonian paths of Xbar.
- Put X = Dbar.
- Checked on all loopless digraphs of order <= 4, all semicomplete and oriented digraphs of order 5, and 110,000
  random digraphs of orders 6 and 7 (`digraph_weights.c`).

**(2).** A cycle using only T-pairs lies in D exactly when its reverse lies in Dbar. Their two permutation cycles
contribute 1 + (-1)^(len-1). All other cycles of D use a B-pair, and all other cycles of Dbar use a Z-pair; neither
kind has a partner.

**(3).** Weights 0 and 2 vanish mod 2 and weights +-1 are 1 mod 2. The defect cycles of D (B-cycles of D and
Z-cycles of Dbar) are the same set as the defect cycles of Dbar. So (3) gives Berge's theorem ham(D) = ham(Dbar)
(mod 2), which it refines.

**(4a).** A Hamiltonian path either ends at n-1, or contains n-1 -> x with x in S. In the second case the path after
n-1 is x followed by any increasing subset of {x+1, ..., n-2}.

**(4b).** The route is Session234's (PROOFS.md DG3, DG4, section 3). Its DG5 already gives the case 7, and its DG3/DG4
give the case 1. The exclusion of 4 is new.
- H is multiplicative over the strong components, which are linearly ordered (Session234 DG3). The B-pair lies in one
  component D_0, and the other components are tournaments with odd H.
- A strong D_0 of order k >= 3 with one B-pair has a strong tournament resolution. This is the b = 1 case of
  Session234 DG4: a strong semicomplete digraph is Hamiltonian, and keeping the cycle's arcs gives a strong spanning
  tournament. A self-contained check: two non-strong resolutions would give cuts (S, S-bar) and (R, R-bar). Their four
  cells force S∩R-bar = {u} and S-bar∩R = {v}, with S∩R dominating and S-bar∩R-bar dominated, so k = 2.
- Hence H(D_0) >= 2 always, and H(D_0) >= f(k) when k >= 3.
- So H(D) is not 1. If H(D) is 4 or 7, the odd cofactor forces H(D_0) in {4, 7}. Then f(k) <= 7 < f(5) = 9 forces
  k <= 4, and the exhaustive census of orders <= 4 has neither value.

**(4c), (5).** Exhaustive enumeration: `digraph_weights.c` (part 2) and `core_weights.py`. The witnesses are in
`out/one_defect_witness.out`.

## Remarks

- **On the owner's conjecture** ("7 or 21 forces a diagonal 1 or a mirrored 00 pair").
  - The diagonal never enters ((2); Session234 DG1), so that clause is vacuous.
  - Session234's semicomplete examples refute the statement.
  - What is true:
    - 7 and 21 force some defect pair (THM-1370);
    - a single 00 pair suffices for both (4c; the 7 case is in Session234's census);
    - a single 11 pair suffices for 21 but never for 7. Both facts are Session234's (DG5 for 7); (4b) adds that it
      never gives 4 either.
- **Structure.** The factor 2 of OCF is the pairing of each odd cycle with its reverse in the complement. Reciprocal
  pairs break that pairing for every cycle through them (weight 1), and missing pairs create unpartnered complement
  cycles (signed weight). This is how 7 = 1 + 6 escapes the tournament exclusion of alpha = (3, 0) (THM-4940 (3),
  THM-343).
