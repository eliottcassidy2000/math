---
id: THM-4840
title: "Arithmetic square-class matroids: graphic and regular boundaries, standard-form TU signing, Hamming15 as equal-parity-degree K6 subgraphs, the corrected even-weight Hamming16 extension, rank-five/six prime banks and specified unit-Cayley tournaments"
status: >
  CORRECTED PROVED (2026-10-28; see correction note and hostile receipts). PROVED (klein-2026-10-28, elementary, with CITED classical inputs named at their point of use; not independently
  audited): (A) given Tutte 1958 (binary and representable over a field of characteristic not 2 implies regular; regular
  iff no F7, F7* minor), Tutte 1959 (graphic iff no F7, F7*, M*(K5), M*(K3,3) minor), the standard-form TU signing theorem (for a prescribed matrix A, regularity is tested on [I | A]) and the classical fact that F7 is representable only in
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

> **CORRECTION 2026-10-28.** Rational representability, signed-graph realizability, and signing a prescribed matrix are different assertions. The added Hamming coordinate records edge-count parity, not common degree parity. Unit halves must be checked after reduction. Exact counterexamples, repaired maps, and a nine-prime M*(K3,3) bank are in [the repair](../../05-knowledge/results/free_bits_seven_fiftysix_20261028/INCOMING_REPAIR.md). The fifteen-column identity and finite rank-five/six banks survive.

Notation. For a finite set `P` of distinct primes, the **square-class matroid** `M(P)` is the column matroid over `F_2`
of `A[l, i] = v_l(p_i + 1) mod 2`; a selection `x` of primes has square `σ(∏ p_i)` iff `A x = 0` (unifying_relations
THEORY §A). `F_7` is the Fano matroid, `F_7*` its dual, `M*(G)` the cographic matroid of `G`.

## (A) The characteristic boundary

1. `M(P)` is graphic iff it has no minor among `F_7, F_7*, M*(K_5), M*(K_{3,3})`.
2. `M(P)` is representable over `Q` iff it is regular iff it has no `F_7` or `F_7*` minor. This is equivalent to the existence of a totally unimodular representation, not to signed-graph incidence representability.
3. A prescribed binary r-by-n matrix A can be signed totally unimodular iff the column matroid of `[I_r | A]` is regular. In standard representation `[I_r | D]` of a regular matroid, D therefore has a TU signing. Regularity of the column matroid of A alone is insufficient.

Proof. The graphic and regular assertions are Tutte's classical excluded-minor and representability theorems. For3, the standard-form signing theorem applies to `[I_r | A]`; conversely adjoining I to a TU signing stays TU. The matrix with columns(111,110,101,011) is a hostile: it represents a triangle plus a coloop, hence a regular matroid, but adjoining I gives Fano, so its original support has no TU signing. All512 sign choices are checked in the repair. Signed-graphic matroids form a proper subclass of regular matroids when restricted to regular examples: [Qin–Slilaty–Zhou](https://doi.org/10.1017/S0963548309990241) give regular excluded minors R15,R16 among others.

Since F7 is representable only in characteristic2, the seven-prime bank has no rational representation preserving the same column dependencies. Integer valuation lifts do exist; reducing their coefficients modulo2 changes their dependencies through torsion, as the [relation-gluing packet](../../05-knowledge/results/relation_gluing_7_15_20261028/THEOREMS.md#t4-changing-coefficients-can-create-dependencies-from-torsion) explains. Changing signed/unsigned equations is not generally Camion signing. ∎

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
code's (enumeration). The extended `[16,11,4]` is obtained by adjoining the bit z equal to the parity of the number of selected edges. Its check columns are `(1,h_e)` for the fifteen quotient columns h_e, and `(1,0)` for the new bit; its dual is `RM(1,4)`. Merely adjoining the all-ones six-vector to the original incidence matrix gives `[16,11,3]`: a triangle survives with new bit0. The common degree bit and the correct even-weight bit are different functions. For the fifteen-prime bank of session85 (least primes
`p_v = d_v k_v² − 1` over the row primes `(2,3,5,7)`: `7, 2, 5, 19, 89, 59, 29, 251, 13, 83, 41, 139, 5669, 419, 839`),
whose matroid is `PG(3,2)`, a selection of primes has square `σ` iff the corresponding edge set of `K_6` has all degrees
of one parity. The correlation of `W(2)` (the outer automorphism of `S_6`) matches each prime with one of the fifteen
square-total triples that are perfect matchings. ∎

## (C) Rank five and six

For every squarefree `d > 1` dividing `2·3·5·7·11` there is a least prime `p = d k² − 1`; the sixteen with `11 | d` are
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
contains `−1`; a *half* is a set `S` of units with `S ∩ −S = ∅`, `S ∪ −S = O_K^×`. The Cayley digraph of `O_K/𝔭` is a tournament precisely when the reduced connection set R satisfies R disjoint(-R) and R union(-R)=all nonzero residues. An abstract half of the units need not have this property: at norm3, `{1,-ω,ω²}` is a half but its image contains both1 and-1.

The following specified choices work. Norm3, `𝔭=(2+ω)`: ω maps to1, so `{1,ω,ω²}` maps to{1} and gives C3. Norm5, `𝔭=(2+i)`: i maps to3, so `{1,i}` maps to{1,3} and gives the regular T5. Norm7, `𝔭=(3+ω)`: ω maps to4, so `{1,ω,ω²}` maps to{1,2,4}, giving P7. At norm4, the original unit group still has halves, but no reduced set can satisfy the tournament condition since -1=1. The symmetrizations are K3,K5,K7 respectively; THM-4839 separately includes the characteristic-two K4. The named pure-blue/Fano comparison survives with these explicit quotient maps. ∎

Scope. (A) is a decision procedure, not a classification of which banks occur (HYP-9438 (1), (2)). (B) and (C) are
finite identities. (D) leaves the ladder `3, 3, 4, …` open beyond `r = 5`. Nothing here bears on the Hamiltonian-path
counts 7 and 21 of THM-1370, on the Borwein indices 7 and 56 of THM-4837, or on `A000568(6) = 56`: those are typed
COINCIDENCES in the packet.

**Update 2026-10-28 (same session, after the owner's go-ahead; FINITE-EXACT, `scripts/bank_search.py`, receipts identical).** The least prime banks (minimal largest prime, exhaustive inside the stated class ranges, isomorphism verified) realizing Tutte's four excluded minors and `R_{10}`: `F_7` at 89 (`7, 2, 19, 5, 89, 59, 29` — the owner's seven primes are the least Fano bank), `F_7*` at 89 (`19, 89, 13, 59, 83, 29, 41`), `M*(K_{3,3})` at 139 (`7, 2, 19, 89, 13, 59, 83, 139, 41`; a second witness with the same maximum is `7, 19, 5, 89, 13, 83, 29, 139, 41`), `R_{10}` at 307 (`19, 5, 43, 13, 59, 131, 139, 41, 307, 109`), `M*(K_5)` at 311 (`7, 19, 43, 103, 131, 139, 41, 311, 109, 181`). So the first bank with no characteristic-zero model appears by `p = 89`, and the smallest regular non-graphic bank has nine primes (`M*(K_{3,3})`), largest 139. HYP-9438 (2) is resolved: the codex repair of the same day (nine-prime bank at 839, nine sharp by Tutte's minors; the MISTAKES.md entry '2026-10-28 — THM-4840: the wrong extra bit and an omitted representation quantifier' (codex)) settled the size, and this search settles the largest-prime minimization it left open.
