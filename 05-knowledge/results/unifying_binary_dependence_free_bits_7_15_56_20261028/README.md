# Binary dependence, free bits and the number six: the theorem behind "graph cycles cannot represent arithmetic dependence", the K₆ reading of the fifteen-prime Hamming bank, rank-5 and rank-6 prime banks, the information principle behind guarded observations, and the instances of 7, 56, 5 and 15, typed

> **CORRECTED 2026-10-28:** the prescribed-matrix signing criterion, the sixteenth parity bit, and the explicit unit quotient maps are repaired below. The unrestricted distinct-prime realization has the one-loop boundary. HYP-9438's ten-element regular-nongraphic minimum is refuted by a nine-prime bank. See [exact repairs and receipts](../free_bits_seven_fiftysix_20261028/INCOMING_REPAIR.md); the fifteen-column identity and rank-five/six computations survive.

**Session:** klein-2026-10-28 (worktree `codex/session-synthesis-20261028`). Canon:
[THM-4840](../../../01-canon/theorems/THM-4840-tutte-camion-dichotomy-for-arithmetic-square-class-banks-hamming-15-is-the-equal-parity-degree-space-of-k6-rank-6-prime-banks-unit-cayley-tournaments-c3-t5-p7-and-fano-chromatic-numbers.md);
open problems: [HYP-9438](../../hypotheses/HYP-9438-every-binary-matroid-is-a-prime-square-class-matroid-arithmetic-r10-and-kuratowski-banks-sum-free-chromatic-numbers-square-sum-descent-and-the-unbounded-dual-word-signature-of-cryptids.md);
reflection: [binary-is-the-right-field-and-six-is-the-right-number](../../../07-reflections/binary-is-the-right-field-and-six-is-the-right-number-klein-20261028.md).

**Status.** PROVED (elementary, with CITED classical inputs named at their point of use: Tutte 1958/1959, Camion
1963, the duad model of `PG(3,2)` through `S_6 ≅ Sp(4,2)`, the `(7,3,1)` difference set, Comtet 1967 / Boas–Wrench
1971): §1, §2, §6.1. FINITE-EXACT: the `K_6` identity by enumeration, the rank-5 and rank-6 prime banks, the chromatic
numbers of `PG(2,2)`, `PG(3,2)`, `PG(4,2)` (`scripts/synthesis_checks.py`, receipts byte-identical). FRAMING: §4 (the
information principle), typed instance by instance. Everything else is typed IDENTITY / ANALOGY / COINCIDENCE /
NUMEROLOGY / OPEN. Not independently audited. No literature-priority claim: the matroid and coding facts are classical;
new here are their application to the arithmetic banks, the `K_6` reading, the rank-5/6 realizations, the chromatic
numbers, the unit-Cayley-tournament identity and the typed ledger.

**Sources read.** [session85_arithmetic_codes_free_bits_20261028](../session85_arithmetic_codes_free_bits_20261028/ARITHMETIC_FREE_BITS.md)
(and its CONSTRAINT_CERTIFICATES, LINEAR_GLUING_INFORMATION), [unifying_relations_20261028](../unifying_relations_20261028/THEORY.md)
(README, THEORY, RECOVERY), [incidence_attractors_20261028](../incidence_attractors_20261028/README.md),
[aliquot_66_five_fifteen_20261028](../aliquot_66_five_fifteen_20261028/README.md), the Borwein ladder
[THM-4837](../../../01-canon/theorems/THM-4837-the-owners-borwein-ladder-first-nonvanishing-of-the-odd-sinc-cosine-transform-at-2-3-and-2pi-closed-form-first-failure.md)
with [HYP-9436](../../hypotheses/HYP-9436-crossing-index-of-the-odd-harmonic-sum-is-floor-exp-2n-minus-gamma-over-4-comtet-type.md),
[THM-4838](../../../01-canon/theorems/THM-4838-von-staudt-clausen-as-reciprocal-sum-thresholds-orbifold-euler-characteristics-triangle-groups-2-3-p-and-the-twelve-behind-the-owners-sets.md),
[THM-4839](../../../01-canon/theorems/THM-4839-unit-cayley-graphs-toroidal-maps-petersen-klein-four-belyi-pythagorean-frey.md),
session85's [FANO_PALEY_INFORMATION.md](../session85_arithmetic_codes_free_bits_20261028/FANO_PALEY_INFORMATION.md) (read after it landed, 2026-10-28),
[THM-1370](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md), [THM-4730](../../../01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md),
[THM-4701](../../../01-canon/theorems/THM-4701-kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs.md),
[THM-4622](../../../01-canon/theorems/THM-4622-logical-status-dictionary-absoluteness-choice-and-automorphisms-for-jc2-collatz-hadwiger-nelson.md),
the classification ledger [PROBLEM_CLASSIFICATION.md](../analogy_reduction_atlas_20261025/PROBLEM_CLASSIFICATION.md).

**Owner's prompt (verbatim).** "take idas and see the bigger pictures and unfiying theorems." followed by the pasted
directive: "four erased selection bits leave two parity-compatible possibilities, but the known coprimality constraint
rules one out. In another context, the same erased positions remain ambiguous between 70 and 364, so the recovery rule
has a sharp boundary. The seven-prime code also extends to explicit 15- and 16-prime constructions. These give
structural links between 7, 15 and 5 through Hamming codes and affine geometry. keep those same ideas from the previos
sessions and merge in the square sum conjecture and collatz and planar jacobian and hadwiger nelson, and also propose
a few of your own open problems on which we may be able to make new progress or extend prior work. also spend this
session synthesizing and exploring and thinking in terms of information looking for tricks to find ways to get a bit
for free. since ordinary graph cycle structures cannot represent ever arithmetic dependence, then work to come up with
new theories and out of the box ideas and unifying theorems that do, and see how they apply to all areas of math we
have studied. The square-class calculation points to binary matroids and error-correcting codes. There is an explicit
seven-prime example whose square-total subsets form the classical Hamming code. That also gives a sharp limit to the
graph analogy: its Fano matroid cannot be represented by an ordinary graph's cycle structure. There is also a
signed-incidence connection between graceful labels and square-sum labels. think {2,5,7,19,29,59,89} and consider
past instances of 7 we have studied in the repo and explore them and their connections extensively looking for any
possible abstract analogy you can notice. notice the connection between those instances of 7, and other past instances
of 5 and 15. think about concepts around prior work on lowest common ancestor and consider ways to go beyond the
paradigm of "guarded observations": what later operations can test determines what information a representation must
retain. We now have precise applications connecting: F = U + S balance, its divisibility boundaries, and a complete
abundance classification. Berggren's tree, {2,3,11}, and square-sum aliquot edges through explicit guarded maps.
Compact arithmetic formulas versus exact memory: Berggren's missing sign needs just one bit, while richer arithmetic
operations can force full reconstruction. also Put 2cos(x) in front of the Borwein integrals and the pattern gets
dramatically longer. a single cosine stretches the run from 7 integrals to 56, and the size of each miss is computed
exactly, as the volume of a corner cut off a many-dimensional box, from nothing but whole-number fractions.
Baillie–Borwein–Borwein (2008), Example 1(b), shows the sinc-sum = sinc-integral identity holds to N = 40248 and fails
at 80499 (Poisson summation at 2π), with Crandall's error estimate. That paper was on my own citation list. the
single-corner first failure (Borwein–Borwein 2001, Cor. 1), the 1/113 value (Almkvist–Gustavsson 2014) and the
threshold ladder (Bäsel–Baillie 2016). think about prior work on 7 and 56 each deeply, trying to find as many examples
as you can for each."

Vocabulary map: "square-total subset" = a set of primes whose product has square `σ` (the repository's "square-sigma"
or "zero-syndrome" selection); "the recovery rule has a sharp boundary" = session85's conditional rule
`x_0 + x_1 = x_4 + x_6` (CONSTRAINT_CERTIFICATES §3); "guarded observations" = session84/85's observers that report the
output and the legality of each future operation; "LCA" = session56's strong-module tree and session85's rank-ownership
formula.

## 1. One theorem behind "graph cycles cannot represent arithmetic dependence" (PROVED with CITED inputs)

Fix distinct primes `P` and the binary matrix `A[l, i] = v_l(p_i + 1) mod 2` (THEORY §A): for squarefree `n_x =
∏_{x_i = 1} p_i`, `σ(n_x) = ∏ (p_i + 1)` is a square iff `A x = 0`. Call the column matroid `M(P)` the **square-class
matroid** of the bank. It is binary by construction, because exponent parity lives in `F_2`.

**Theorem 1 (the characteristic boundary).**

1. `M(P)` is the cycle matroid of an ordinary graph iff it has no minor isomorphic to `F_7`, `F_7*`, `M*(K_5)` or
   `M*(K_{3,3})` (Tutte 1959, CITED).
2. `M(P)` is representable over Q iff regular iff it has no F7 or F7* minor (Tutte, CITED). This gives some TU representation, not necessarily a signed-graph incidence model.
3. A specified matrix A has a TU signing iff the binary matroid of `[I | A]` is regular. For a standard representation `[I | D]` of a regular matroid, D has such a signing. The arbitrary-matrix version formerly stated here is false; columns(111,110,101,011) provide the exact counterexample.

**Consequences.** The seven-prime bank `{2,5,7,19,29,59,89}` is `F_7` itself (THEORY §A1), representable over no field
of characteristic `≠ 2`; it has no characteristic-zero representation preserving the same column dependencies. Integer valuation lifts with different rational dependencies remain valid. The Kuratowski obstructions of THM-4701 (`K_5`, `K_{3,3}`, whose cographic matroids are two
of Tutte's four excluded minors) and the Fano obstruction of the seven-prime bank are the two halves of one theorem:
for arithmetic banks, Kuratowski and Fano are the same kind of object. The six-prime bank `{2,3,5,7,11,13}` (columns
`010, 000, 110, 100, 010, 101` on `(2,3,7)`) is graphic: a loop at 3, a parallel pair `{2,11}`, two triangles sharing
the edge `5–7`; the "positive multigraph model" of unifying_relations is this graph, and Theorem 1 says exactly where
such a model stops. The right category for "relations plus composition laws" is binary matroids (equivalently binary
codes with their duals); Theorem 1 is the decision procedure for when the graph picture is recoverable.

## 2. The number six: Hamming [15,11,3] is the equal-parity-degree space of K₆ (PROVED, enumerated)

The fifteen-prime bank of session85 (least primes `p_v = d_v k_v² − 1` for the fifteen squarefree `d_v` over the row
primes `(2,3,5,7)`: `7, 2, 5, 19, 89, 59, 29, 251, 13, 83, 41, 139, 5669, 419, 839`) has square-class matroid
`PG(3,2)`, and its square-total selections form Hamming `[15,11,3]`. The classical duad model (`S_6 ≅ Sp(4,2)`, CITED)
makes this a statement about six vertices: map the duad `{a,b}` of a six-set to `e_a + e_b` in the even-weight subspace
of `F_2^6` modulo the all-ones vector; this is a bijection from the fifteen edges of `K_6` to the fifteen nonzero
vectors of `F_2^4`. Hence

**Theorem 2.** A set of edges of `K_6` is a Hamming `[15,11,3]` codeword iff all six vertex degrees have the same parity.
The code is the cycle space of `K_6` (all degrees even, dimension 10) together with its all-odd coset. The 35 lines of
`PG(3,2)` are the 20 triangles and the 15 perfect matchings; the 105 weight-four words are the 45 four-cycles and the 60
claw-plus-disjoint-edge subgraphs; the full weight distribution `1, 35, 105, 168, 280, 435, 435, 280, 168, 105, 35, 1`
(weights `0, 3, …, 12, 15`) agrees with the Hamming code's by enumeration (`data/synthesis_checks.out`). Equivalently
`PG(3,2) = M(K_6)/𝟙`: the projective geometry is the graphic matroid of `K_6` with the all-ones vector adjoined and
contracted. Adjoining the edge-count parity bit gives the extended code `[16,11,4]`: use check columns(1,h_e) for the fifteen quotient columns and(1,0) for the added bit. Adjoining the all-ones six-vector to the original incidence matrix instead records the common degree bit and gives minimum distance3, since triangles survive. These are different one-bit extensions.


**What this buys.** (i) *The 7/15/5 link in one picture.* `7 = PG(2,2)` is the Fano plane, which is `K_4` with its four
vertices and three perfect matchings as points; `15 = PG(3,2)` is the duad geometry of six points; `5 = dim RM(1,4) = 1
+ 4` is the affine lift. The hidden number is six. (ii) *The outer automorphism of `S_6` is an arithmetic self-duality.*
The correlation of the generalized quadrangle `W(2)` swaps duads and synthemes, so it matches each of the fifteen primes
with one of the fifteen square-total triples that are perfect matchings (the odd weight-three codewords); it is not an
automorphism of the bank, it is a bijection between its elements and half of its circuits. (iii) *Ramsey.* `R(3,3) = 6`
says every two-colouring of the edges of `K_6` has a monochromatic triangle, hence `PG(3,2)` is not two-colourable; the
Fano plane is the smallest 3-uniform hypergraph without Property B. In arithmetic words: **every partition of the seven
primes `{2,5,7,19,29,59,89}` into two classes has a class containing three primes whose product has square `σ`**, and
likewise for the fifteen-prime bank with triples that are triangles or matchings of `K_6`. (iv) *Tournaments at order
six.* The 15 arc positions of a 6-tournament are the points of `PG(3,2)` (IDENTITY), the 56 isomorphism classes
(A000568(6)) and the strong floor `f(6) = 15` (THM-1370) are counts on the same six vertices but with no mechanism
linking them to the code (COINCIDENCE, recorded so they are not re-promoted).

## 3. Rank-5 and rank-6 prime banks (FINITE-EXACT, new) and the realization conjecture

For row primes `(2,3,5,7,11)` the least prime `p = d k² − 1` exists for all 31 squarefree `d`; the new sixteen (those with
`11 | d`) are `(d,k,p) = (11,2,43), (22,3,197), (33,2,131), (66,2,263), (55,6,1979), (110,1,109), (165,2,659),
(330,2,1319), (77,2,307), (154,9,12473), (231,8,14783), (462,1,461), (385,6,13859), (770,1,769), (1155,6,41579),
(2310,1,2309)`; largest prime 41,579, largest `k = 9`. For row primes `(2,3,5,7,11,13)` all 63 exist, largest prime
1,261,259, largest `k = 42` (`data/synthesis_checks.out`). Since the squarefree part of `p + 1` is exactly `d`, the
31 and 63 primes realize `PG(4,2)` and `PG(5,2)` as square-class matroids, so **every simple binary matroid of rank at
most six is a restriction of an explicit prime bank**, extending session85's rank `≤ 4`. The all-rank statement needs,
for each squarefree `d`, one prime of the form `d k² − 1`: this is a Bunyakovsky/Schinzel-type statement for the
quadratic `d k² − 1`, OPEN in general and conditional on Schinzel's hypothesis H (HYP-9438 (1)); parallel elements
(two primes with the same `d`) need infinitely many such primes for one `d`. The prime 3 is the unique loop
(`p + 1` square forces `p = 3`), so loops are the one structure a prime bank cannot repeat.

### 3b. The least prime banks for Tutte's four excluded minors and for R₁₀ (FINITE-EXACT, added after the owner's go-ahead)

`scripts/bank_search.py` takes every prime below `10^6` (78,498 primes; 60,956 distinct square classes `d = sqf(p+1)`;
the loop class `d = 1` only at `p = 3`), keeps for each class its least prime `L(d)`, and searches exhaustively inside
`D(N) = {d : L(d) < N}` for configurations isomorphic to each target, minimising the largest prime (any bank can be
replaced class by class with least primes, so this is the true minimum). Fano: all seven nonzero points of a 3-space.
`F_7* = AG(3,2)` minus a point: seven points of a coset of a 3-space avoiding `0`. `M*(K_{3,3})`: the `3×3` grid
`{a_1, a_2, a_1a_2, b_1, b_2, b_1b_2, a_1b_1, a_2b_2, a_1a_2b_1b_2}` on four independent classes (the hyperbolic quadric
of `PG(3,2)`), every row-to-row matching tried. `R_{10}` (affine, rank 5): grown from four-circuits through hyperplanes
with at least six points to five-spaces with at least ten, then every affine-hyperplane section. `M*(K_5)` (affine, rank
6): its only four-circuits are the five vertex stars of `K_5`, pairwise meeting in one element, so the search walks the
four-circuit hypergraph. Every candidate is checked by rank, circuit profile and an explicit backtracking isomorphism
to a reference representation; receipts byte-identical.

| target | elements, rank | least max prime | square classes `d` | primes |
|---|---|---:|---|---|
| `F_7` (Fano) | 7, 3 | **89** | 2, 3, 5, 6, 10, 15, 30 | 7, 2, 19, 5, 89, 59, 29 |
| `F_7*` (`AG(3,2)` minus a point) | 7, 4 | **89** | 5, 10, 14, 15, 21, 30, 42 | 19, 89, 13, 59, 83, 29, 41 |
| `M*(K_{3,3})` | 9, 4 | **139** | 2, 3, 5, 10, 14, 15, 21, 35, 42 | 7, 2, 19, 89, 13, 59, 83, 139, 41 |
| `R_{10}` | 10, 5 | **307** | 5, 6, 11, 14, 15, 33, 35, 42, 77, 110 | 19, 5, 43, 13, 59, 131, 139, 41, 307, 109 |
| `M*(K_5)` | 10, 6 | **311** | 2, 5, 11, 26, 33, 35, 42, 78, 110, 182 | 7, 19, 43, 103, 131, 139, 41, 311, 109, 181 |

Minima are unique as values, not as witnesses: `M*(K_{3,3})` at 139 is also realized by `7, 19, 5, 89, 13, 83, 29, 139, 41`
(classes `2, 5, 6, 10, 14, 21, 30, 35, 42`), found by an earlier traversal order. Exhaustive ranges: `D(200)`, `D(400)`, `D(6000)`, `D(15000)`, `D(5400)` respectively (a construction gave an `M*(K_5)`
bank at 5,393, so the minimum lies in `D(5394)`; the `R_{10}` search found 611 candidate ten-sets, the `M*(K_5)` search
9,503). Readings. (i) **The owner's seven primes are the least Fano bank**: no seven primes with largest prime below 89
have square classes forming `F_7`. (ii) Both Fano minors live below 89, so **the first prime bank with no
characteristic-zero model appears by `p = 89`**; the `F_7*` bank is seven of the eight classes `{5, 10, 14, 15, 21, 30,
35, 42}` (an `AG(3,2)` coset of `2·{...}`-type), sharing `5, 10, 15, 30` with the Fano bank. (iii) **The smallest regular
but non-graphic prime bank has nine primes, not ten** (HYP-9438 (2) said ten; `M*(K_{3,3})` has nine elements; the codex repair of the same day already refuted the ten with a nine-prime bank of largest prime 839 and proved nine sharp by Tutte's minors), with
largest prime 139; `R_{10}` and `M*(K_5)` follow at 307 and 311, four orders of magnitude below the `10^6` the
hypothesis allowed. (iv) All five least banks reuse the small classes `5, 10, 14, 15, 30, 35, 42` (primes 19, 89, 13,
59, 29, 139, 41): the square classes of `p + 1` for small primes are dense in the low-rank geometries, which is why the
excluded minors appear so early. Typed: FINITE-EXACT within the stated ranges; the all-rank realization question
(HYP-9438 (1)) is unchanged.

## 4. The information principle, stated once (FRAMING with exact instances)

"What later operations can test determines what information a representation must retain" has a precise form that the
repository reached from three sides this week: THEORY §E1 (the coarsest lawful quotient, Myhill–Nerode with legality as
an observable), session85's linear criterion (a target functional is recoverable iff it lies in the row span of
constraints plus observations; the exact lost-bit budget `d = dim W − dim W_local`, and its junction-tree form
`d = Σ_bags − Σ_separators − dim C`, a conditional mutual information), and its nonlinear form (recoverability iff
`x_k + g ∈ I` in the Boolean coordinate ring). One sentence covers all three:

> **A representation must retain the annihilator of the future test algebra; a bit is free exactly when its test is
> dependent on constraints already retained; the price of a test is its new rank, and nonlinear legality is a
> Boolean-ideal membership, not a rank.**

"Beyond guarded observations" is then a change of object: stop listing guards and compute the dual space (the
annihilator) of the tests the later operations can perform; the LCA tree of session56/85 is one factorization of that
dual space along separators (a guard is charged only for new rank at the LCA of its support). The repository's one-bit
phenomena sort cleanly:

| phenomenon | the test | status of the bit | mechanism |
|---|---|---|---|
| Hamming-7 overall parity from any Fano line (session85) | dual weight-four word | free (dependent) | the dual code |
| 70 versus 364 (session85) | four erased coordinates, one conflict edge | conditionally free: free iff `x_4 + x_6 = 0` | nonlinear legality, Boolean ideal `x_0 = x_0(x_0+x_1) + x_0 x_1` |
| Berggren's sign `χ = (−1)^{#A}` (THM-4808) | residues at every modulus | not free from residues; a two-state observer | a word character, not a residue function |
| Rédei: `H(T)` is odd | parity of an involution's non-fixed points | zero bits (a tautology) | fixed-point parity of an involution |
| Terras: parity vectors of Collatz | cylinders `mod 2^k` | zero free bits (full shift) | bijective coding; no dependent test |
| Antihydra halting `Σ(2v_i − 3) = −6` | a linear form of unbounded support | never recoverable by a bounded scope | the dual word's support exceeds every bag |
| the all-ones projection `K_6 → PG(3,2)` (§2) | global degree parity | the one bit the graph model must forget | elementary quotient by `𝟙` |
| the affine lift `[15,11,3] → [16,11,4]` | the sixteenth element | one element buys one unit of distance | adjoin the parity of the selected edge count |
| the bipartite sign gauge (THEORY §D) | edge sums versus edge differences | one free root parameter on bipartite graphs; fixed, with `Z/2` torsion, otherwise | Camion signing exists iff regular |

The last two rows are the same move seen from both sides: adjoin or contract the all-ones vector. The Antihydra row is
the bridge to the halting atlas: a `Π⁰₁` cryptid is a halting predicate whose certificate is a dual word of unbounded
support, so no finite-scope observer family (no regular or bounded-Parikh decider) can see it — the information-theoretic
form of the three-level ladder in
[eth-false-busy-beaver-sizes-and-modular-loops](../../../07-reflections/eth-false-busy-beaver-sizes-and-modular-loops-the-machine-collatz-dictionary-klein-20261027.md).

## 5. The Borwein ladder in the same frame (CITED and typed)

THM-4837: `7, 56, 418` are the crossing indices `m(ξ)` of the odd harmonic sum `1 + 1/3 + … + 1/(2m+1)` at `ξ = 2, 3, 4`,
and HYP-9436 (certified for `2 ≤ n ≤ 3000`) gives `m(n) = ⌊e^{2n−γ}/4⌋`, the odd analogue of Comtet's `⌊e^{n−γ}⌋`. The
factor `2cos x` doubles the first sinc's budget (`2 cos x · sinc x = 2 sinc 2x`), so the run stretches from the threshold
`2` to `3`; each unit of threshold multiplies the run by `e² ≈ 7.39` (`7 → 56 → 418 → 3091`). The first miss is the
corner of the box `∏[−1/(2j+1), 1/(2j+1)]` cut off by the hyperplane `Σ y_j = ξ`: `π ε^m (2m+1)!!/(2m)!!` with
`ε = O_m − ξ` (Borwein–Borwein 2001, Cor. 1). In the vocabulary of §4 the ladder is a *budget*, not a code: there is no
dual word, and the 7 and 56 here are COINCIDENCES with the Fano 7 and the tournament 56 (§6).

## 6. The instances of 7, 56, 5 and 15, typed

**6.1 One identity ties the sevens of this week together (IDENTITY, elementary; the Fano/Paley facts are INHERITED).**
That the Fano lines are the translates of the quadratic-residue difference set `{1,2,4}` and hence the out-neighbourhoods
of the Paley tournament `P_7` is already in canon-adjacent notes
([collatz_mod6_20260922_paley_fano_octonion_design.md](../collatz_mod6_20260922_paley_fano_octonion_design.md) S6–S7,
[arithmetic_braids2_20260917_fano_code.md](../arithmetic_braids2_20260917_fano_code.md)), and session85's
[FANO_PALEY_INFORMATION.md](../session85_arithmetic_codes_free_bits_20261028/FANO_PALEY_INFORMATION.md) (pushed while this
note was being written) owns the finer statement: a fixed Fano plane admits exactly 16 compatible Paley orientations (three
partner bits and one chirality bit), and flipping the three arcs of a line when its prime is selected transports square
`σ` exactly to "the result is Paley". What is added here is the ring reading. Let `π = 3 + ω` in `Z[ω]`, of norm 7,
so `Z[ω]/(π) ≅ F_7` with `ω ↦ 4`. The cube roots of unity `{1, ω, ω²}` map to `{1, 2, 4}`, the quadratic residues mod 7,
which is the `(7,3,1)` difference set (re-checked: every pair of points on exactly one line). So **`P_7` is the Cayley tournament
of `Z[ω]/(3+ω)` with connection set the cube roots of unity, the Fano plane is its out-neighbourhood design, the Heawood
graph (THM-4839's toroidal map at norm 7) is the Fano incidence graph, and `K_7 = {3,6}_{(2,1)}` (THM-4839 (a), (b)) is
the symmetrization of `P_7` and the dual of the Heawood map.** The same recipe at the other norms gives the pure-blue
alphabet of THM-4730: at norm 3, `Z[ω]/(2+ω) ≅ F_3` with all six units mapping onto `{±1}` and the half `{1, ω, ω²}`
onto `{1}`, the Cayley tournament is `C_3`; at norm 5, `Z[i]/(2+i) ≅ F_5` with `i ↦ 3`, the half `{1, i}` maps to
`{1, 3}` and the Cayley tournament is the regular 5-tournament `T_5`. **The unit-Cayley tournaments of the Eisenstein and
Gaussian integers at the primes of norm 3, 5, 7 are `C_3`, `T_5`, `P_7`**: the first two are exactly the strong pure-blue
atoms (THM-4730), and the third is the Fano plane. At norm 4 (`Z[ω]/(2) ≅ F_4`) there is no reduced tournament half because `−1 = 1`, which
is why `K_4` has no regular tournament and the alphabet stops at `T_5`. The group behind the seven-prime bank
(`GL(3,2)`, the Fano collineation group, order 168) is the group behind the triangle group `(2,3,7)` of THM-4838 (`X(7)`
is Klein's quartic, automorphism group `PSL(2,7) ≅ GL(3,2)`): the finite-geometry and modular families that HYP-9437's
lexical census found anti-associated are structurally one group at `7`. Inside it, `F_{21} = Z_7 ⋊ Z_3 = Aut(P_7)`
(THM-4743) is the normalizer of the Singer cycle and acts regularly on the 21 flags of the Fano plane: the
`7`/`21` of the Fano world are points and flags, `Z_7` and `F_{21}`.

**6.2 Sevens.** Fano plane, Hamming `[7,4,3]`, the seven-prime bank, `P_7`, `F_{21}`, `GL(3,2)`, `(2,3,7)`, Klein's quartic,
Heawood, `K_7` on the torus (all one object, 6.1); the Borwein crossing `m(2) = 7` (first failure at `1/15`;
COINCIDENCE); `χ(R²) ≤ 7` by the hexagonal 7-colouring and the Moser spindle's seven vertices (THM-4622's bounds
`5 ≤ χ_alg ≤ χ(R²) ≤ 7`; COINCIDENCE with Fano's `7`/`5` access structure); THM-1370's forbidden Hamiltonian-path count
`7` (no mechanism; COINCIDENCE, as the seven_twentyone notes already typed); `vSC(6) = {2,3,7}` and `χ(2,3,7) = −1/42 =
−B_6` (THM-4838; the `7` is the Hurwitz triangle group, IDENTITY with 6.1 through `PSL(2,7)`); the divisor-height
five-cycle `4 → 2 → 5 → 3 → 7 → 4` (Session43; `3^7 ∥ R`; NUMEROLOGY); `7 = 111_2`, `2^3 − 1`, Catalan's `2^3 + 1 = 3^2`
(the Mersenne-line notes; separate mechanism).

**6.3 Fifty-sixes.** `56 = C(8,3)` is the number of triples of the eight points of `AG(3,2)`, each in a unique affine plane
(the 14 weight-four words of the extended Hamming code `[8,4,4] = RM(1,3)`, Steiner system `S(3,4,8)`): this is the
`7 → 8` affine lift, the exact analogue of session85's `15 → 16` (IDENTITY within the Hamming family). `56 = A000568(6)`,
the tournaments on six vertices (COINCIDENCE with the lift; it is a count on the same six vertices as §2 but by a
different mechanism). `56 = m(3)` of the Borwein ladder (COINCIDENCE). `56 = C(2·4, 3)` in THM-3010's bronze sequence
`1, 4, 15, 56, 210 = C(2n, n−1)` (an identity of binomials with the previous row, `15 = C(6,2)` likewise; no mechanism
connecting Newton-circuit alternation to codes: COINCIDENCE). `56` as the class number `h(−4N_6)` on the square-sigma ray
(aliquot_66 P4; non-monotone, NUMEROLOGY). Outside the repository: the 56 minuscule weights of `E_7`, the 56 lines on a
degree-two del Pezzo surface (`2 × 28` bitangents), `|PSL(2,7)|/3`; the last is the only one with a thread to 6.1 and it
is not pursued here.

**6.4 Fives and fifteens.** `5 = dim RM(1,4)`, the five checks of `[16,11,4]`, and the Hamming-7 access structure (overall
parity from every five-coordinate set, never from two, from a four-set exactly when it contains a line: session85 §1) — IDENTITIES
of §2. `15 = C(6,2)`, the points of `PG(3,2)`, the perfect matchings of `K_6`, the edges of the Petersen graph (`L(K_5)`'s
crossing graph needs five hierarchies: LCA note), the odd weight-three codewords — one object (§2). `f(6) = 15`
(THM-1370's strong floor) and `A000568(6) = 56` — same six vertices, no mechanism (COINCIDENCE). `5 ≤ χ(R²)` (de Grey)
and the five Fermat primes, `R_5 = T_5` (6.1: the norm-5 unit-Cayley tournament), the five-cycle of divisor heights,
`|σ^{−1}(144)| = 5`, Ogg's fifteen supersingular primes (aliquot_66 ATLAS §5) — distinct mechanisms; typed COINCIDENCE
except where 6.1 says otherwise. The repository's one earlier structural `15/25` pair, `f(6) = 15, f(7) = 25`
(THM-1370), stays as it was.

**6.5 What the repository already held (inventory by a reader over canon titles, the results and hypotheses indices,
TANGENTS and the two atlases; file pointers as given there).** Sevens and twenty-ones: `H ≠ 7, 21` at every order with the
multiplicative mechanism `21 = 3·7` forcing a strong component with `h ∈ {7, 21}` (THM-115, THM-517, THM-1370); the
presentation fibres `tc = H/|Aut|` never `7` or `21` (THM-4812); `|Aut(P_7)| = 21`, frozen at `F_{21}` along the Mersenne
doubling tower from order 7 to 255 (THM-448, THM-4557); `AGL(1,7)` on `QR_7` typed "the genuine 7/21 group link"
(results INDEX); `C(7,2) = 21` forcing the `x^5` coefficient of every 7-tournament skew polynomial (THM-1455); the Paley
spectrum `x(x²+7)³` with Gauss sum `±i√7` (THM-1450) — all consistent with 6.1 (the Gauss sum is the `S² = J − 7I`
side of THM-4814 (E)); `P_7`'s cyclic triples form two Fano planes (TANGENTS T106), the two `GL(3,2)`-orbits of
session85's sixteen orientations. Fifty-sixes: the 56 classes of 6-tournaments are the 56 metagraph vertices at order
six (METAGRAPH-ATLAS; 12 self-converse + 44 chiral, THM-1420; switching quotient 6, THM-1415); `P_7` has 63 minimum
blocking 6-arc sets, 56 Hall obstructions and 7 exotic stars (THM-4553); `P_7`'s invariant triple `(H, disc, s·disc) =
(189, 8, 56)` (THM-1950); LRC(14)'s best near-blocker misses exactly 56 cells, typed "numerology with structure" there
(T372); the ballot column `1, 4, 15, 56, 210` (THM-3010). Fifteens: `2^15` labelled 6-tournaments (LEM-003); the maximum
`H` at order 5 is 15 (THM-2010); `f(6) = 15` (THM-1860, THM-4741); 15 prime tournaments at order 6 (THM-1960);
`K(A_5) = 15` matching the fifteen involution axes (THM-870); the square-sum chain first exists at 15 (THM-4505,
THM-4510); the lacunary threshold "sharpened from 15 to 7" (THM-933); the `q = 15` sunflower as a "contracted Fano
stalk" (THM-1242). Fives: `A_5` cannot be a tournament symmetry (THM-868); `K_5` Kuratowski threshold (THM-519);
the untouchable `5` with its `{7, 21}` link typed "no map" (results INDEX); the quintic threshold `C_5` (THM-436). The
reader's two inferences not stated in the repository ("56 = 35 strong + 21 non-strong" from HYP-9429's counts; the
`15 & 56` adjacency in THM-3010 having no mechanism) are recorded here as COINCIDENCE pending a mechanism. Against this
ledger, THM-4840 adds exactly four structural statements: the characteristic boundary, `PG(3,2) = M(K_6)/𝟙`, the rank-5/6
banks, and the unit-Cayley reading `C_3, T_5, P_7`; everything else in this section was already typed.

## 7. The requested merges, typed

- **Square-sum conjecture (Hamiltonian paths in the square-sum graph on `{1..n}`; onset 15, 23, 25 inherited from
  aliquot_66 ATLAS §6 / THM-4510).** THEORY §D gives the exact reconstruction: a path is bipartite, so its square edge
  sums determine the vertex labels up to one free root parameter, `x_v = ε_v t + c_v`. Hence a square-sum Hamiltonian path
  on `[n]` is the same thing as a sequence of `n − 1` squares whose alternating partial sums, shifted by one integer `t`,
  enumerate `[n]`: the information content is the square sequence plus one integer, and the "bit for free" is that the
  parity of `t` is forced by the parity of `n`'s labels. What is OPEN is a descent in the sense of unifying_relations'
  third theory (a unique inverse with decreasing height, as Berggren's) on these sequences; HYP-9438 (4).
- **Collatz.** Terras's parity vector is a full shift (zero free bits); the only exact recovery equations are
  Archimedean linear forms in the Parikh vector (Steiner–Eliahou), outside the `F_2` test algebra of cylinders; Antihydra's
  criterion `Σ(2v_i − 3) = −6` has unbounded support (§4). Nothing here changes a Collatz status; it names the information
  reason for the difficulty.
- **Planar Jacobian.** By Bass–Connell–Wright the inverse of a Keller map of degree `d` in `n` variables has degree at most
  `d^{n−1}`, so in the plane the certificate "`∃ G` of degree `≤ d` with `G ∘ F = id`" is a bounded-scope test, and `JC(2)`
  at each degree is a finite linear-algebra statement (true through degree 100, Moh 1983, CITED). In the vocabulary of
  §4, `JC(2)` asserts that the local guard (`det J = 1`) is a dependent test of the bounded-scope certificate; the
  dimension-three refutation in canon (the owner's Keller map) is a case where the local guard is independent of global
  bijectivity. ANALOGY only; the matroid content is nil.
- **Hadwiger–Nelson.** `5 ≤ χ(R²) ≤ 7` (THM-4622) shares its integers with Fano's `7` points and `5`-set access structure
  by COINCIDENCE. The structural analogue is Property B: the Fano plane is the smallest 3-uniform hypergraph that is not
  two-colourable, as the Moser spindle is the smallest unit-distance graph that is not three-colourable; the chromatic
  numbers of the line hypergraphs of `PG(2,2), PG(3,2), PG(4,2)` are `3, 3, 4` (computed; equivalently the least number
  of sum-free sets partitioning `F_2^r ∖ 0` is `2, 3, 3, 4` for `r = 2..5`), so the arithmetic banks have Ramsey numbers
  of their own (HYP-9438 (3)). Unit-distance dependence itself is graphic (closed unit-step walks are cycles), so
  Theorem 1 says the code machinery does not transfer to the plane.

## 8. Open problems proposed (HYP-9438)

(1) every binary matroid is a prime square-class matroid (conditional on Schinzel H; rank `≤ 6` done here); (2) the
arithmetic `R_{10}` and the arithmetic cographic Kuratowski banks — RESOLVED in §3b (least maximal primes 307, 311, 139;
nine elements, not ten); (3) the sum-free chromatic ladder
`χ(PG(r−1,2)) = 2, 3, 3, 4, …` and its arithmetic Ramsey reading; (4) a Berggren-type descent for square-sum paths;
(5) the exact Nerode quotient of guarded square-sigma multiplication as a function of the matroid and the conflict
graph (session85's 86/87/96); (6) the unbounded-dual-word signature as a formal criterion separating cryptids from
decidable machines; (7) whether the duad/syntheme self-duality of §2 is induced by an arithmetic involution on the
fifteen primes.

## 9. Reproduction and scope

`scripts/synthesis_checks.py` (standard library; ~90 s, dominated by the `PG(4,2)` three-colourability refutation) prints
the `K_6` weight distributions, the Fano/`P_7` check, the three prime-bank tables and explicit colourings;
`receipts/` holds the `python3` and `python3 -O` outputs (byte-identical). Nothing in §4 is a theorem beyond the three
exact instances it cites; nothing here closes a conjecture; the owner's "seven-prime code extends to 15- and 16-prime
constructions" is session85's and is cited, not re-derived.
