---
id: THM-4840
title: "The Tutte–Camion dichotomy for arithmetic square-class banks: the binary matroid of the square classes of sigma(p) = p+1 over a prime set has an ordinary-graph cycle model iff it avoids the four Tutte minors F7, F7*, M*(K5), M*(K3,3), has ANY integer or rational incidence model (signed, signless, gain) iff it is regular iff it avoids F7 and F7*, and can be signed to a totally unimodular matrix iff regular; the seven-prime Fano bank therefore has no incidence model in characteristic zero at all, and the Kuratowski obstructions of THM-4701 and the Fano obstruction are two halves of one theorem. Hamming [15,11,3] is the set of edge subsets of K6 whose six degrees share one parity (PG(3,2) = M(K6)/1, lines = 20 triangles + 15 perfect matchings, 105 weight-four words = 45 four-cycles + 60 claw-plus-edge), so the fifteen-prime bank's square-total selections are the equal-parity-degree subgraphs of K6 and the sixteenth element of [16,11,4] is the global degree parity. Rank-5 and rank-6 prime banks exist (31 and 63 least primes d k^2 - 1), so every simple binary matroid of rank <= 6 is a square-class matroid of distinct primes. The line hypergraphs of PG(2,2), PG(3,2), PG(4,2) have chromatic numbers 3, 3, 4 (least sum-free partitions of F2^r minus 0 for r = 3, 4, 5), hence every two-colouring of {2,5,7,19,29,59,89} has a monochromatic triple with square sigma. The unit-Cayley tournaments of Z[omega] and Z[i] at the primes of norm 3, 5, 7 are C3, T5, P7."
status: >
  PROVED (klein-2026-10-28, elementary, with CITED classical inputs named at their point of use; not independently
  audited): (A) given Tutte 1958 (binary and representable over a field of characteristic not 2 implies regular; regular
  iff no F7, F7* minor), Tutte 1959 (graphic iff no F7, F7*, M*(K5), M*(K3,3) minor), Camion 1963 (a binary matrix signs
  to a totally unimodular matrix iff its matroid is regular) and the classical fact that F7 is representable only in
  characteristic 2; (B) the duad model of PG(3,2) (S6 = Sp(4,2), CITED) plus a direct proof, with the weight distribution
  and the line and weight-four decompositions FINITE-EXACT by enumeration of the 2^15 edge sets of K6 and of the 2^15
  Hamming syndromes; (C) FINITE-EXACT: the 31 and 63 primes p = d k^2 - 1, each the least prime for its squarefree d over
  the row primes (2,3,5,7,11) and (2,3,5,7,11,13), all distinct, extending session85's rank <= 4 banks; (D) FINITE-EXACT
  by exhaustive colouring search: chromatic numbers 3, 3, 4 of the line hypergraphs of PG(2,2), PG(3,2), PG(4,2), with
  explicit colourings; (E) elementary. INHERITED: the square-class matroid and the seven-prime Fano bank
  (unifying_relations THEORY A, A1; incidence_attractors A3), the fifteen-prime bank and the [16,11,4] lift
  (session85 ARITHMETIC_FREE_BITS sections 3 and 6), the Fano lines as the Paley7 out-neighbourhood design
  (collatz_mod6_20260922_paley_fano_octonion_design S6-S7, arithmetic_braids2_20260917_fano_code) and the 16 Paley
  orientations of a fixed Fano plane (session85 FANO_PALEY_INFORMATION), the pure-blue alphabet K1, C3, T5 (THM-4730),
  THM-4839's complete unit Cayley graphs K4, K5, K7. OPEN: the all-rank realization (HYP-9438 (1)) and the chromatic
  ladder beyond r = 5 (HYP-9438 (3)).
session: klein-2026-10-28 (codex/session-synthesis-20261028)
source: 05-knowledge/results/unifying_binary_dependence_free_bits_7_15_56_20261028/README.md
depends_on:
  - 05-knowledge/results/unifying_relations_20261028/THEORY.md
  - 05-knowledge/results/session85_arithmetic_codes_free_bits_20261028/ARITHMETIC_FREE_BITS.md
  - 01-canon/theorems/THM-4701-kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs.md
  - 01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md
related:
  - 05-knowledge/hypotheses/HYP-9438-every-binary-matroid-is-a-prime-square-class-matroid-arithmetic-r10-and-kuratowski-banks-sum-free-chromatic-numbers-square-sum-descent-and-the-unbounded-dual-word-signature-of-cryptids.md
  - 01-canon/theorems/THM-4839-unit-cayley-graphs-toroidal-maps-petersen-klein-four-belyi-pythagorean-frey.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 05-knowledge/results/session85_arithmetic_codes_free_bits_20261028/FANO_PALEY_INFORMATION.md
scripts:
  - 05-knowledge/results/unifying_binary_dependence_free_bits_7_15_56_20261028/scripts/synthesis_checks.py (+ data/synthesis_checks.out, receipts/)
---

# THM-4840 — binary is the field of arithmetic dependence, six is the number behind fifteen

Notation. For a finite set `P` of distinct primes, the **square-class matroid** `M(P)` is the column matroid over `F_2`
of `A[l, i] = v_l(p_i + 1) mod 2`; a selection `x` of primes has square `σ(∏ p_i)` iff `A x = 0` (unifying_relations
THEORY §A). `F_7` is the Fano matroid, `F_7*` its dual, `M*(G)` the cographic matroid of `G`.

## (A) The characteristic boundary

1. `M(P)` is graphic iff it has no minor among `F_7, F_7*, M*(K_5), M*(K_{3,3})`.
2. `M(P)` is representable over `Q` (equivalently has a signed-incidence, signless-incidence, gain-graph or totally
   unimodular model) iff it is regular iff it has no `F_7` or `F_7*` minor.
3. A binary representation of `M(P)` signs to a totally unimodular `±1` matrix iff `M(P)` is regular; this signing is
   the passage from square-sum (unsigned) to graceful (signed) incidence of unifying_relations README §"Signed
   incidence", and it exists exactly for regular banks.

Proof. 1 and 2 are Tutte's theorems applied to the binary matroid `M(P)`, with the remark that a binary matroid
representable over some field of characteristic `≠ 2` is regular (Tutte 1958); 3 is Camion's theorem. Since `F_7` is
representable only in characteristic 2, the seven-prime bank `{2,5,7,19,29,59,89}`, whose square-class matroid is `F_7`
(THEORY A1), has no incidence model over `Z` or `Q`; the six-prime bank `{2,3,5,7,11,13}` is graphic (a loop at 3, a
parallel pair `{2,11}`, two triangles on the edge `5–7`). ∎

The two cographic Kuratowski minors are the matroids behind THM-4701's `K_5` and `K_{3,3}`; so for arithmetic banks the
planarity obstructions and the Fano obstruction are the four excluded minors of a single theorem.

## (B) Hamming [15,11,3] on six vertices

Map the edge `{a,b}` of `K_6` to `e_a + e_b ∈ F_2^6` modulo the all-ones vector `𝟙`; the even-weight subspace modulo `𝟙`
is `F_2^4`, and the fifteen edges go bijectively to the fifteen nonzero vectors (two distinct edges differ by a vector of
weight 2 or 4, never `0` or `𝟙`). Hence the parity-check matrix of Hamming `[15,11,3]` (all nonzero vectors of `F_2^4`)
is the edge–vertex incidence matrix of `K_6` reduced modulo `𝟙`, and an edge set is a codeword iff its degree vector
is `0` or `𝟙` mod 2, i.e. **all six degrees have the same parity**. The code is the cycle space of `K_6` (dimension 10)
together with the coset of all-odd subgraphs; `PG(3,2) = M(K_6)/𝟙`. Lines (weight 3): the 20 triangles and the 15
perfect matchings; weight 4: the 45 four-cycles and the 60 subgraphs `K_{1,3} ∪ K_2`; the whole weight distribution
`1, 35, 105, 168, 280, 435, 435, 280, 168, 105, 35, 1` at weights `0, 3, 4, …, 12, 15` coincides with the Hamming
code's (enumeration). Adjoining `𝟙` itself as a sixteenth element gives `[16,11,4]`, its five checks being the four
reduced vertex parities and the global parity, dual `RM(1,4)`. For the fifteen-prime bank of session85 (least primes
`p_v = d_v k_v² − 1` over the row primes `(2,3,5,7)`: `7, 2, 5, 19, 89, 59, 29, 251, 13, 83, 41, 139, 5669, 419, 839`),
whose matroid is `PG(3,2)`, a selection of primes has square `σ` iff the corresponding edge set of `K_6` has all degrees
of one parity. The correlation of `W(2)` (the outer automorphism of `S_6`) matches each prime with one of the fifteen
square-total triples that are perfect matchings. ∎

## (C) Rank five and six

For every squarefree `d` dividing `2·3·5·7·11` there is a least prime `p = d k² − 1`; the sixteen with `11 | d` are
`(d,k,p) = (11,2,43), (22,3,197), (33,2,131), (66,2,263), (55,6,1979), (110,1,109), (165,2,659), (330,2,1319),
(77,2,307), (154,9,12473), (231,8,14783), (462,1,461), (385,6,13859), (770,1,769), (1155,6,41579), (2310,1,2309)`;
for `d | 2·3·5·7·11·13` all 63 exist (largest 1,261,259, largest `k = 42`). The squarefree part of `p + 1` is `d`, so
these primes realize `PG(4,2)` and `PG(5,2)` as square-class matroids, and every simple binary matroid of rank `≤ 6` is a
restriction of one of them (every simple binary matroid of rank `r` is a restriction of `PG(r−1,2)`). Loops: `p + 1`
is a square only for `p = 3`, so a prime bank has at most one loop. ∎

## (D) Chromatic numbers and the arithmetic Ramsey statement

A colouring of the points of `PG(r−1,2)` with no monochromatic line is a partition of `F_2^r ∖ 0` into sum-free sets.
Exhaustive search: `PG(2,2)` needs 3 colours (not 2), `PG(3,2)` needs 3 (not 2; a 3-colouring: `{1,2,4,8,15}, {3,5,7,10,11},
{6,9,12,13,14}` in the integer labels of `F_2^4`), `PG(4,2)` needs 4 (not 3). Through (A1 of THEORY) and (B): every
two-colouring of `{2,5,7,19,29,59,89}` has a monochromatic triple whose product has square `σ`; every two-colouring of the
fifteen-prime bank has a monochromatic triangle or perfect matching of `K_6` (this case is also `R(3,3) = 6`); three
colours suffice in both; four are needed for the rank-5 bank. ∎

## (E) Unit-Cayley tournaments

Let `K ∈ {Q(ω), Q(i)}` and `𝔭` a prime of `O_K` of norm `N ∈ {3, 5, 7}`. The unit group `O_K^×` has order 6 or 4 and
contains `−1`; a *half* is a set `S` of units with `S ∩ −S = ∅`, `S ∪ −S = O_K^×`. The Cayley digraph of `O_K/𝔭` with
connection set the image of a half is a tournament iff the image avoids `−1`-pairs, which holds exactly when `−1 ≠ 1` in
the residue field, i.e. `N` odd. Norm 3 (`𝔭 = (2+ω)`, all units `↦ ±1`, the half `{1,ω,ω²} ↦ {1}`): `C_3`. Norm 5
(`𝔭 = (2+i)`, `i ↦ 2`, half `{1, i} ↦ {1, 2}`): the regular 5-tournament `T_5`. Norm 7 (`𝔭 = (3+ω)`, `ω ↦ 2`, half
`{1, ω, ω²} ↦ {1, 2, 4}` = the quadratic residues = the `(7,3,1)` difference set): `P_7`. At norm 4 (`Z[ω]/(2) ≅ F_4`)
`−1 = 1` and no half exists, matching the absence of a regular 4-tournament. THM-4839's complete unit Cayley graphs
`K_4, K_5, K_7` are the symmetrizations; THM-4730's pure-blue alphabet `K_1, C_3, T_5` is the norm `1, 3, 5` part of this
list, and `P_7` (the Fano out-neighbourhood design, inherited) is where the alphabet stops being pure-blue. ∎

Scope. (A) is a decision procedure, not a classification of which banks occur (HYP-9438 (1), (2)). (B) and (C) are
finite identities. (D) leaves the ladder `3, 3, 4, …` open beyond `r = 5`. Nothing here bears on the Hamiltonian-path
counts 7 and 21 of THM-1370, on the Borwein indices 7 and 56 of THM-4837, or on `A000568(6) = 56`: those are typed
COINCIDENCES in the packet.
