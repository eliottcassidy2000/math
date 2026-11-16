# Collatz on the honeycomb, the three tilings and the five solids, and the pattern "3 into 5" taken apart

Session klein-2026-11-16-hex. A synthesis of prior sessions (a read-only survey of 47 items; the files are named
where used) with the row picture of path_splice_rows_continuation_20261116/COLLATZ_ROWS.md and the tree tournament of
COLLATZ_TOURNAMENT.md. Every statement below is labelled PROVED (cited or shown here), FINITE-EXACT, ANALOGY or
NUMEROLOGY.

## 1. Collatz as a hexagonal tiling of edges — what the repository already has, and what the rows add

- **The honeycomb cluster of 2026-11-02** (`collatz_hexagonal_growth_20261102`, `session140_collatz_hexagonal_geometry_20261102`,
  `honeycomb_gap_aliquot_20261102`, `clocks_carriers_genus_20261102/COLLATZ_GEOMETRY.md`, `joint_incidence_clocked_frames_20261102/COLLATZ_LATTICE.md`,
  `eleven_wedge_honeycomb_20261102`): the ordinary Collatz graph has a faithful directed subdivision inside the honeycomb
  lattice, but its routing lengths must be unbounded (PROVED, independently reviewed); the capacity obstruction is
  exact — 704 integers reach 5 in at most 21 steps while a honeycomb ball of that radius has 694 vertices (the deficit
  is ten, MISTAKES.md); the cycle `1 → 4 → 2 → 1` becomes a literal hexagon after two-edge subdivision.
- **The brick-wall models** (`session178_prime_flow_honeycomb_20261104/honeycomb/RESULTS.md`): with the chart
  `Φ(n) = (a, q − 1)` for `n = 2^a q`, every even arrow `n → n/2` is a nearest-neighbour honeycomb edge and the odd
  arrows are long jumps (constructive, FINITE-EXACT).
- **The hexagon of residues** (THM-4521, `collatz_coherence_20260930_trichotomy_hexagon_motifs.md`): the level-2 unit
  cycle is the hexagon `(Z/9)^× = {1, 2, 4, 8, 7, 5}` (PROVED).
- **The two self-similarities** (`hn_interiors_atoms_tilings_20261113` §4): in `Z[ω]` the prime 3 ramifies with index 3
  and 2 is inert with index 4, so `3x + 1` and `x/2` are the two primitive self-similarities of the hexagonal tiling
  (CITED dictionary).

What the row picture adds (PROVED, elementary): the brick-wall chart's rows `{2^a q : a ≥ 0}` are the doubling chains
of `q`, i.e. in count coordinates the rows `L_c = {1 + c 2^s}` of COLLATZ_ROWS.md are the chart of `N − 1`. The even
arrows (nearest-neighbour edges) walk along a row; the odd arrows (long jumps) are the Syracuse map on row indices,
`q ↦ oddpart(3q + 1)`. So "Collatz as a hexagonal tiling of edges" and "odd numbers in rows of chains growing by
doubling" are one drawing: rows are the inert direction (×2 stays in the row), jumps are the ramified direction
(×3 + 1 changes the row), and the 704/694 obstruction is the statement that the jumps cannot be made short. In that
drawing the forbidden pair `{7, 21}` is one jump: the head of `L_3` to the second rung of `L_5`.

## 2. Three tilings, five solids: the pattern taken apart

The regular tilings `{3,6}, {4,4}, {6,3}` and the Platonic solids `{3,3}, {3,4}, {4,3}, {3,5}, {5,3}` are the integer
solutions of one Schläfli inequality, `(p − 2)(q − 2) = 4` and `< 4` (HYP-3772; `hn_interiors_atoms_tilings_20261113` §3
tabulates them). The pattern "3 splits into 5" is: **one self-dual object stays put (`{4,4} ↦ {3,3}`) and one dual pair
splits into two dual pairs (`{3,6}/{6,3} ↦ {3,4}/{4,3}` and `{3,5}/{5,3}`, the 6 dropping to 4 or 5)**. Counting:
`3 = 1 + 2·1`, `5 = 1 + 2·2`. The reflection `the-angle-defect-spine-…` records the scope correction that the five
solids and the five planar lattice types are different classifications (the five Bravais lattices are not this
pattern). Deconstructed, the observation has three ingredients: (i) a classification by a boundary inequality,
(ii) a duality with fixed points, (iii) a count of the form `fixed + 2 · pairs`. The same three ingredients recur in
the programme:

| where | the inequality / classification | the duality | fixed + pairs |
|---|---|---|---|
| order-4 tournaments (THM-5124's diamond `1, 3, 3, 5`; `procgen_kuratowski_20260925`) | counts of the four classes | converse (arc reversal), `H(T^op) = H(T)` | `TT_4` (self-converse, 1), the 3-cycle-with-source / with-sink pair (3, 3), the strong one (self-converse, 5): `1 + 2 + 2` — THM-5124's "only additive identity" `5 = 1 + 2 + 2` (PROVED) |
| the five odd numbers around the Collatz image `3q` of a count `q` | the digits `σ ∈ {−2, −1, 0, 1, 2}` of `3q + 2σ` | `x ↦ −x`, which conjugates `3x + 1` to `3x − 1` | centre `3q` (Collatz, digit 0); pairs `{3q ± 4}` (the `3x − 1` image `3q − 4` and the digit-2 child `3q + 4 = 3(q + 1) + 1`) and `{3q ± 2}` (`3(q − 1) + 1` and `3(q + 1) − 1`): one member of each pair is a `+1`-type image, the other a `−1`-type image (PROVED identities, COLLATZ_ROWS.md) |
| the rows `L_3`, `L_5` (`hn_interiors` §3 "3\|5 of Chan–Pak") | the two subordinate lines `c = 3 = F_0`, `c = 5 = F_1`; each loses exactly one point, `7 = 1 + 3·2` and `21 = 1 + 5·4` | the Syracuse step `3 ↦ 5` | "3 into 5" is literally the only Collatz step inside the exceptional set: `7 ↦ 21`, the head of `L_3` to the second rung of `L_5` (PROVED, COLLATZ_ROWS.md) |
| the fixed-root continuation | `D = E − O ≥ −g` (FINITE-EXACT to order 8) | converse swaps root-start and root-end sums `(E, O) ↔ (E_e, O_e)` | the self-dual quantity `X = E + O` survives a sink; the pair `(D, D_e)` is what a sink kills (CONTINUATION.md) |
| the Busch floors (THM-4744/4930) | `f(3q) = 3·5^(q−1)`, `f(3q + 1) = 5^q`, `f(3q + 2) = 9·5^(q−1)` | — | three residues of the order, one factor 5 per three vertices: "every 3 vertices split off a 5" (PROVED, CITED) |

The mindset, stated as a method: when a small integer pattern appears (3 and 5, 8 and 9, 24), look for the inequality
whose boundary it sits on and the involution that pairs its solutions; the count `fixed + 2·pairs` then predicts the
next member. Applied to tournaments of order 5: the twelve classes split under converse into self-converse classes and
pairs; their counts `H` are converse-invariant, so each count value is attained by a `fixed + 2·pairs` number of
classes — the Hamiltonian spectrum's multiplicities carry the same signature (FINITE-EXACT reading of THM-1370's
census; not pursued further here).

## 3. Platonic solids, lattices and tournaments: typed verdicts

- `THM-4839`: the honeycomb modulo `2 + ω` (norm 3) is `K_(3,3) = {6,3}_(1,1)`, modulo `3 + ω` (norm 7) the Heawood graph;
  Petersen, `K_5`, `K_(3,3)` are the smallest non-planar skeletons of types `{5,3}`, `{4,4}`, `{6,3}` (PROVED). This is
  the honest link between the three tilings and the small non-planar graphs, and through THM-4839's unit Cayley
  tournaments (`C3/T5/P7`) the link to tournaments: the Paley tournament on 7 vertices is a quotient of the `{6,3}`
  honeycomb story at norm 7. ANALOGY beyond that sentence.
- `collatz_pentagon_operator_coherence_20260930`: the tetrahedron and the icosahedron are the fixed points of the induced
  cycle operators on the Platonic graphs — the "fixed + pairs" shape again (FINITE-EXACT).
- `gilbreath_fermat_platonic_20260926`, five Fermat primes against five solids: NUMEROLOGY (no map).
- The tree tournament of COLLATZ_TOURNAMENT.md is the lattice-free version of the honeycomb picture: it keeps the tree
  and the global order and forgets the embedding; its invariant (`H = 1` iff depth-first) is what survives when the
  honeycomb's bounded routing is given up.

## 4. More patterns generated by the method (CONJECTURES and QUESTIONS, labelled)

- **Q1.** Does the order-`n` converse pairing explain the multiplicities of the small counts in the Hamiltonian
  spectrum — is the number of classes of count `H` congruent to the number of self-converse classes of count `H` modulo 2?
  (True for `n = 4`: counts 1 and 5 once each, 3 twice.)
- **Q2.** In the five-number window `{3q − 4, …, 3q + 4}`, the digit programme realises the right half `3q, 3q + 2,
  3q + 4` with one continuation; is the left half `3q − 4, 3q − 2, 3q` realised by the converse construction (reverse
  the splice through `x ↦ −x`)? If so the ternary deposit has a mirror image and the five-window is one object.
- **Q3.** The tree tournament's count sequence for the breadth-first Collatz order: is `H(T_N)` asymptotically
  `C·ρ^N` with `ρ` determined by the branching rate of the inverse tree (two children at `n ≡ 4 mod 6`), i.e. by the
  honeycomb's growth, and does the even/odd preference change only `C`?
