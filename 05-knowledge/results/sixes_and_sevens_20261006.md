# Sixes and sevens: stripping one square is the knight's Hall law on every torus tested; the 6x6 and 7x7 knight tori are the two-place product and the Paley torus; one point-stripping ladder (PSL(2,7) > Borel > torus) carries Paley, THM-4524, the Collatz Frobenius and the octonionic J; the S6 manuscript and the Lyapunov n = 6 case

**Session:** mac-mini-2026-10-06-sixseven (worktree `math-wt-chessboard-20261006`), 2026-10-06. This continues the chessboard session ([chessboard_weave_20261006.md](chessboard_weave_20261006.md)).

**Owner's prompt (verbatim core):** "think how 'The knight analogue: on the 6×6 torus, the fewest moves whose removal kills every closed tour is 7, and the only way is to strip one square.' relates to our previous work on 6 and 7 and https://alpo.ge/s6.pdf". The pasted context (the 11-square Lean formalization, its degree-8 polynomial and Langlands, and arXiv:2610.06783 against "18 and 19") was answered in parallel by opus S16 in [square11_octic_langlands_collatz_20261006.md](square11_octic_langlands_collatz_20261006.md). Section 7 below only adds what that note does not have.

**Canon:** [THM-4552](../../01-canon/theorems/THM-4552-exceptional-knight-tori-six-is-a-two-place-product-seven-is-a-paley-torus.md) and [THM-4553](../../01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md). **Hypotheses:** [HYP-9211](../hypotheses/HYP-9211-knight-torus-blocking-number-is-seven-by-stars-only.md), [HYP-9212](../hypotheses/HYP-9212-lyapunov-symmetric-maximizer-holds-at-order-six.md).

**Status.**

* **PROVED (elementary):**
  * the CRT isomorphism `G_6 ≅ C_4 × Paley(9)`;
  * the `F_49` description of `G_7`;
  * the transversality lemma;
  * the point-stripping ladder;
  * the Fano/octonion pairing;
  * the 3SUM-against-cycle-search comparison (section 7.3).
* **FINITE-EXACT:**
  * knight-torus blocking numbers and star-only minima, `n = 5, 6, 7, 8`;
  * automorphism groups, `n = 5..14`;
  * linear stabilisers, `n = 5..26`;
  * the Paley `P_7` blocking census;
  * the octic class group (PARI `bnfcertify`).
* **LEAN-CHECKED:** the `P_7` facts, in [sixseven_20261006_paley_seven_certificate.lean](../../04-computation/lean/standalone/sixseven_20261006_paley_seven_certificate.lean). It uses core Lean 4.30 and `native_decide`; each theorem rests on `propext` and its own native axiom, with no `sorry`.
* **FINITE-NUMERICAL:** the Lyapunov `n = 6` searches, with positive controls at `n = 7, 8, 9` (HYP-9212).
* **ANALOGY / NUMEROLOGY:** typed where they occur.
* **Audit:** the independent audit is in section 10.

**Scripts:**

* [structure](../../04-computation/experiments/sixseven_20261006_structure.py) (`.out` beside it, ending `ALL CHECKS PASSED`; needs nauty's `dreadnaut`);
* [N x N blocking census in C](../../04-computation/experiments/sixseven_20261006_knight_torusblock_n.c), with the CP-SAT lazy-cut cross-check [knight_torus_block.py](../../04-computation/experiments/sixseven_20261006_knight_torus_block.py) and the census output [knight_torusblock_n.out](../../04-computation/experiments/sixseven_20261006_knight_torusblock_n.out);
* Lyapunov: [L-BFGS, KV branch, stripping](../../04-computation/experiments/sixseven_20261006_lyapunov_six.py) and [Adam searches with controls](../../04-computation/experiments/sixseven_20261006_lyapunov_adam.py), output [sixseven_20261006_lyapunov.out](../../04-computation/experiments/sixseven_20261006_lyapunov.out);
* [octic class group, PARI](../../04-computation/experiments/sixseven_20261006_octic_class_group.gp).

## 0. Answers in brief

1. **The sentence is not about 6.** "Seven deletions, and only by stripping a square" holds on every knight torus tested:
   * `n = 5, 6, 7, 8` exhaustively, stars only. On `8×8`, for example, all `8,637,487,551` six-sets through a fixed move leave a closed tour, and of the `359,895,314,625` seven-sets exactly 14 block, all stars.
   * The 7 is `8 − 1`, the knight's degree minus one. It is the knight's Hall law, the same one HYP-9168 states for tournaments.
   * Conjecture: this holds for all `n ≥ 5` (HYP-9211).
2. **But 6 and 7 are exactly the exceptional knight tori, in opposite ways.**
   * **`G_6` is a two-place product.** The 6×6 knight graph is, by CRT, the tensor product `C_4 × Paley(9)` of a mod-2 rook step and a mod-3 bishop step. Squares `x` and `x + (3,3)` are twins with identical neighbourhoods, and `|Aut| = 2^18 · 144`. This is the shape of note 18's "two places" (`n < 6^L` is fixed by its residues mod `2^L` and `3^L`) at the first level, typed ANALOGY in 2.1.
   * **`G_7` is a Paley torus.** On `7×7 = F_49`, the eight knight moves are the elements of norm 5, a coset of the norm-one torus `μ_8`, and all of them are non-squares. The nightrider graph is the complement of Paley(49), and the queen graph is Paley(49). The queen and nightrider slopes split `P^1(F_7)` four and four, by `QR_7` against `NQR_7`. Among `n >= 6`, only here does the knight torus gain an order-8 "45-degree rotation" `ω = 5 + 5i` (stabiliser of order 16). It is the same quadratic character as note 19's `QR_7 | NQR_7` split of the trivial cycle's codes, an EXPLAINED COINCIDENCE (2.2).
   * 7 is also the least `n > 1` with `gcd(n, 30) = 1`, so the four knight direction classes are pairwise transversal (pairwise determinants 3, 4, 5).
   * **Langlands-adjacent: `G_7` is a Kloosterman graph.** It is the finite Euclidean graph "squared distance 5" over `F_7`. Its nontrivial eigenvalues are exactly the Kloosterman sums `-Kl_7(c)`, so Weil's bound `2 sqrt 7` makes it Ramanujan. The Ramanujan knight tori are exactly `n = 5, 6, 7, 8, 10` (section 2.3).
3. **The repo's 6/7 phenomena are one ladder of point-strippings** (THM-4553). The steps:
   * `PSL(2,7)` acts 2-transitively on the 8 points of `P^1(F_7)`, so it preserves no tournament.
   * Strip `∞`. The stabiliser is the Borel subgroup (order 21) `= Aut(P_7)`, whose only invariant tournaments are `P_7` and its reverse.
   * Strip `0`. The stabiliser is the split torus `{x ↦ 2^k x}` (order 3) `= Aut(P_7 − 0)`. Here `P_7 − 0` is THM-4524's unique all-odd 6-tournament.

   The same order-3 group is note 19's Frobenius, the Collatz map on its trivial cycle `{1, 4, 2}`. On the octonion side, the three Fano lines through a point pair its out-neighbours `QR_7` with its in-neighbours `NQR_7` by `q ↦ 3q`. That pairing is the almost-complex structure `J` of `S^6` on the tangent space at a basis point. That the knight's degree 8 equals `|P^1(F_7)|` is NUMEROLOGY: the knight's "strip one square" is a degree count, not this ladder.
4. **Tournaments are less rigid.** In the Paley tournament `P_7` the Hamiltonian-path blocking number is 6, its Hall bound `N − 1`. Of its 63 minimum blocking sets:
   * 56 are Hall obstructions: 21 "two sources", 21 "two sinks", and 14 "three vertices with one common out- or in-neighbour";
   * 7 are stars. A star blocks by isolating a vertex but leaves a 1-path-cycle factor (the vertex plus two cyclic triangles), so it is an "exotic minimum" in the chessboard note's sense.

   On the knight torus the star *is* the Hall (2-factor) obstruction, and nothing else blocks. The count 63 and the star are Lean-checked; the type split is checked in Python.
5. **The S6 manuscript.** The relation is a typed ANALOGY: the manuscript puts all of `χ(S^6) = 2` in one singular fibre (repo ledger: `χ(W) = 2`; THM-3991's one-Euler-fibre grammar), and the knight torus puts its entire 2-factor deficiency at one square. Both are local–global index counts with one nonzero local term. Nothing transfers in either direction. The manuscript's status in the repo, MANUSCRIPT CLAIM / UNDER AUDIT, is unchanged.
6. **Lyapunov, 6 against 7 (HYP-9212).** The search uses Kressner–Vandereycken's own discovery recipe (Adam on the gap, which escapes the "equality ridge" that L-BFGS cannot leave). It finds order-7 counterexamples in 37 of 1712 runs, order-8 in 25 of 100, and order-9 in 4 of 9, but **0 of 6000 at order 6**.

   Along the KV branch, the 7x7 margin found by continuation decays to 0 as the matrix is pushed toward `A_6 ⊕ 0` (log-slope 1.2–1.5, optimizer-dependent), and every stripped 6x6 block falls short. The strongest order-7 counterexamples share KV's shape (five large singular values and two small), but weaker ones do not, so there is no `5 + 2` mechanism.

   This is evidence, not proof: `n = 6` stays OPEN, now recorded as HYP-9212 ("the conjecture holds at 6").

## 1. The knight sentence holds on every torus tested (FINITE-EXACT; HYP-9211)

`G_n` is the knight graph on `Z_n × Z_n`. For `n ≥ 5` its 8 moves are distinct mod `n`, so `G_n` is 8-regular. It is edge-transitive under translations and the point group `D4`, which acts regularly on the 8 moves. So every deletion set may be assumed to contain a fixed move `e0`. Let `β(n)` be the fewest deletions leaving no Hamiltonian cycle. The upper bound `β(n) ≤ 7` holds because deleting 7 of the 8 moves at a square leaves it with degree 1.

| `n` | 6-sets through `e0` (all leave a tour) | 7-sets through `e0` that block | of which strip a square | time |
|---|---|---|---|---|
| 5 | 71,523,144 | 14 (of 1,120,529,256) | 14 | 1.3 min |
| 6 | 464,306,843 | 14 (of 10,679,057,389) | 14 | (chessboard session) |
| 7 | 2,231,243,664 | 14 (of 70,656,049,360) | 14 | 1.3 min + 30 min |
| 8 | 8,637,487,551 | 14 (of 359,895,314,625) | 14 | 4 min + 2 h (split runs) |

The method is a pool of about 1000 random Hamiltonian cycles. Nested loops keep, at each depth, the sub-pool avoiding the prefix, and the last edge is tested against the AND of the sub-pool masks. Only uncertified sets go to an exact DFS. Every DFS cycle is re-verified and added to the pool.

The cross-check is a CP-SAT lazy-cut hitting-set model ([knight_torus_block.py](../../04-computation/experiments/sixseven_20261006_knight_torus_block.py)): the master problem is a minimum hitting set of the pool's `Z_n^2 ⋊ D4` orbits, and the oracle is `AddCircuit`. It independently gives `β(5) = 7` with stars only.

Each star leaves 14 blocking 7-sets through `e0`: 2 endpoints of `e0`, times 7 choices of which move survives. In total there are `8 n^2` minimum sets. On `6×6` the Hall (2-factor) bound is also 7, by Ore's bipartite `f`-factor criterion (chessboard note, section 7). For even `n`, the same count gives `d ≥ 6(|S| − |T|) + 1 ≥ 7`. Equality needs `|N(T)| = |T| + 1` with exactly 8 edges leaving `N(T)`; `T = ∅` and `T = B − b` both give the star. Whether any intermediate `T` is tight is left to the census, which finds only stars.

**Reading.** The sentence's 7 is `8 − 1` and its "strip one square" is the min-degree Hall obstruction. Neither depends on the board being 6×6. The sentence is the knight case of the repo's Hall genus (HYP-9168: the cheapest way to kill every Hamiltonian object is to starve a set).

## 2. But 6 and 7 are the exceptional knight tori (PROVED + FINITE-EXACT; THM-4552)

### 2.1 `G_6` is a two-place product

Reduce mod 2 and mod 3. The 8 knight moves map onto `S_2 × S_3` exactly:

* `S_2 = {(1,0), (0,1)} ⊂ Z_2^2`, the rook step mod 2;
* `S_3 = {±(1,1), ±(1,−1)} ⊂ Z_3^2`, the bishop step mod 3;
* `8 = 2 · 4`.

Hence

`G_6 ≅ Cay(Z_2^2, S_2) × Cay(Z_3^2, S_3) = C_4 × Paley(9)` (tensor product),

with `Paley(9) = K_3 □ K_3`. Since `C_4 = K_{2,2}` has twins, `x` and `x + (3,3)` have the same 8 neighbours in `G_6`. There are 18 twin pairs and `|Aut(G_6)| = 2^18 · 144`.

This is the "two places" of note 18 ([collatz_connectivity_from_rigidity_20261001.md](collatz_connectivity_from_rigidity_20261001.md), section 7: "a positive `n < 6^L` is determined by its residues modulo `2^L` and `3^L`") at level `L = 1`.

**Typed.** The decomposition is an exact MAP at `n = 6`. As a bridge to Collatz it is an ANALOGY: the knight's 2-part and 3-part are a rook and a bishop, not halving and tripling, and the product does not persist to `6^L`.

### 2.2 `G_7` is a Paley torus

Identify `Z_7^2` with `F_49 = F_7[i]` (−1 is a non-residue mod 7), with norm `N(a, b) = a^2 + b^2`.

* The 8 knight moves are exactly the elements of norm 5, the coset `(1 + 2i) μ_8` of the norm-one torus `μ_8 = ker(N)`. Since `5 ∈ NQR_7`, all are non-squares of `F_49^*`.
* The nightrider steps (all multiples of knight moves) are `{N ∈ NQR_7}`, and the queen's steps are `{N ∈ QR_7}`, the squares. The torus graphs are Paley(49) (queen) and its complement (nightrider), both `srg(49, 24, 11, 12)`. Together they are `K_49`: the 8 slopes of `P^1(F_7)` split as `{0, ∞, ±1}` (queen, square norms) and `{2, 3, 4, 5}` (knight, non-square norms).
* `ω = 5 + 5i` has order 8 and maps knight moves to knight moves. It is a "45-degree rotation", and the vertex stabiliser of `G_7` is `μ_8 ⋊ Gal = D_8`, of order 16.

**Linear stabiliser census** (`n = 5..26`, FINITE-EXACT). The knight set has linear symmetries beyond `D4` only for `n ∈ {5, 6, 7, 8, 10}`. Those of order 8 occur only for `n = 5` (degenerate: `G_5 = K_5 □ K_5`, since `2^2 = −1` mod 5 puts the moves on the two isotropic lines) and `n = 7`. The extras at 6 and 10 are CRT swaps, and at 8 they are shears by 4.

Full automorphism groups (nauty), with generic value `8 n^2`:

| `n` | 5 | 6 | 7 | 8 | 9 | 10 | 11–14 |
|---|---|---|---|---|---|---|---|
| \|Aut\|/n^2 | 1152 | 2^20 | 16 | 512 | 8 | 16 | 8 |

**Transversality lemma (PROVED).** The pairwise determinants of the four knight directions `(1,2), (2,1), (1,−2), (2,−1)` are `±3, ±4, ±5`. So every two direction classes of `n`-cycles form a grid (each cycle of one class meets each cycle of the other once) iff `gcd(n, 30) = 1`. The least such `n > 1` is 7.

Compare the toroidal `n`-queens problem, solvable iff `gcd(n, 6) = 1` (Pólya, CITED), and Euler's 36 officers: there are no two orthogonal Latin squares of order 6, so no affine plane of order 6 (Tarry, CITED). Six fails every one of these field-like structures, and seven has all of them.

**Typed bridge to note 19.** Note 19 ([collatz_paley_bridge_20261001.md](collatz_paley_bridge_20261001.md)) reads the trivial Collatz cycle as the code `QR_7` (real place) and `NQR_7` (2-adic place). The 7×7 torus reads the queen as `QR_7` and the knight as `NQR_7`. Both splits are exact. The coincidence of labels is an EXPLAINED COINCIDENCE: both are the quadratic character of `F_7`. There is no map between Collatz cycles and leaps.

### 2.3 Seven is a Kloosterman graph, hence Ramanujan (PROVED; CITED Weil)

`G_7 = Cay(F_49, {N = 5})` is the finite Euclidean graph `E_7(5)` of Medrano–Myers–Stark–Terras: the squares at squared distance `5 = 1^2 + 2^2`. The prime 7 is the only one where the knight set is the whole sphere (`p + 1 = 8` points).

* **Spectrum.** Its eigenvalue at the character `(a_1, a_2)` is `-Kl_7(3 N(a))`, with `Kl_7(c) = Σ_x e((x + c/x)/7)`. So apart from 8, the spectrum is the six Kloosterman sums over `F_7` (values `-2.049, 2.357, 1.604, 2.692, -4.494, -1.110` for `c = 1..6`), each with multiplicity 8.
* **Ramanujan.** Weil's bound `|Kl_q| <= 2 sqrt q` makes `G_7` a Ramanujan graph, because the Ramanujan bound `2 sqrt(d - 1)` for `d = 8` is Weil's `2 sqrt q` for `q = 7`.
* **Census.** The Ramanujan knight tori are exactly `n ∈ {5, 6, 7, 8, 10}`. These are the same five `n` that carry extra linear symmetries (2.2), a finite coincidence with no mechanism claimed. For `n >= 12` the character `(1,0)` already exceeds `2 sqrt 7`.

**This is the Langlands-adjacent face of the knight question (typed MAP, exact).** `Kl_7(c)` is the trace of Frobenius on Deligne's rank-2 Kloosterman sheaf, and the Weil bound is its Riemann hypothesis. Kloosterman sheaves are geometric-Langlands eigensheaves (Heinloth–Ngô–Yun). So the expansion of the 7x7 knight torus is literally a Riemann hypothesis over `F_7`.

Compare the repo's Langlands ladder in S16's note:

* GL(1): the Collatz arithmetic and the octic's sign character;
* GL(2): `X_0(11)`;
* GL(7): the octic's open Artin representation.

The 7x7 knight torus adds a function-field rung: a rank-2 local system on `G_m/F_7` (geometric monodromy `SL_2`, Katz). Over function fields, Langlands for GL(2) is a theorem (Drinfeld), and Heinloth–Ngô–Yun construct the Kloosterman case explicitly. The one statement used here, the Weil bound, is proved. Nothing follows for tours or for Collatz.

## 3. One ladder of point-strippings (PROVED + FINITE-EXACT + LEAN-CHECKED; THM-4553)

| step | object | symmetry | repo item |
|---|---|---|---|
| 8 points | `P^1(F_7)` | `PSL(2,7)`, order 168, 2-transitive, so no invariant tournament | the 8 slopes of the 7x7 torus (2.2); the queen/knight 4+4 split is **not** `PSL(2,7)`-invariant (stabiliser of order 4; 8 in `PGL(2,7)`) |
| strip `∞` (7 points) | Paley `P_7` on `F_7` | Borel `{x ↦ ax + b : a ∈ QR_7}`, order 21 `= Aut(P_7)`; its 2 orbitals give exactly `P_7` and its reverse | flip-rank apex, LRC Paley heptagon |
| strip `0` (6 points) | `P_7 − 0` | split torus `{x ↦ x, 2x, 4x}`, order 3 `= Aut(P_7 − 0)` | THM-4524's unique all-odd 6-tournament |
| the torus acting on `QR_7` | `{1, 2, 4}` | `x ↦ 2x` | note 19 Prop. 5: Collatz `T` on its trivial cycle |

So the Frobenius `<x2>`, which by note 19 (Proposition 5) is the part of `Aut(P_7)` that Collatz realises on its trivial cycle, is exactly the full automorphism group of THM-4524's 6-tournament. The 6 non-trivial translations that the carry kills make up, with the identity, the unipotent radical of the Borel. Note 19's "maximal symmetry against no symmetry" is, in this language, Borel against torus.

**Modular curves (CITED, standard; typed DICTIONARY).** The Borel is `Γ_0(7)` mod 7, and the torus is `Γ_0(7) ∩ Γ^0(7)` mod 7, conjugate to `Γ_0(49)`. The ladder `8 -> 7 -> 6` is therefore the tower `X(1) <- X_0(7) <- X_0(49)`, ending at the genus-1 CM curve 49a1 (CM by `Q(sqrt(-7))`). `QR_7`, the code of the trivial Collatz cycle, is the set of residues of the primes split in `Q(sqrt(-7))`. In this reading, the symmetry the trivial cycle keeps is the one uniformising a CM elliptic curve, a GL(1)-over-`Q(sqrt(-7))` Langlands object. This is consistent with S16's finding that all exact Collatz arithmetic is abelian. No Collatz consequence follows.

**Octonions (PROVED).** Take the octonion table `e_x e_{x+1} = e_{x+3}` (indices mod 7; the alternative law was checked exactly on random elements). The lines are translates of `{0, 1, 3}`, and `{1, 2, 4}` is one of them. The three lines through 0 are `{0,1,3}, {0,2,6}, {0,4,5}`. They pair each out-neighbour `q ∈ QR_7` of 0 in `P_7` with the in-neighbour `3q ∈ NQR_7`. Left multiplication `J = L_{e_0}` on `T_{e_0}S^6 = span(e_1, ..., e_6)` sends `e_q ↦ e_{3q}` and `e_{3q} ↦ −e_q`. So the almost-complex structure of `S^6` at a basis point is the Fano matching of the six stripped points, `QR → NQR`. The torus `x ↦ 2x` permutes the three complex lines cyclically.

Within `P_7 − 0`:

* `QR_7` and `NQR_7` are cyclic triangles;
* the Fano matching arcs run `QR → NQR`;
* the only `NQR → QR` arcs are the antipodal arcs `u → −u`, the orbit proved odd in THM-4524 (recorded in HYP-9167).

**Typed.** All identities are exact. "S^6 is the octonionic shadow of THM-4524's tournament" is an ANALOGY: no parity statement about Hamiltonian paths follows from `J`, and `J` is not integrable (the repo's Nijenhuis check, 2026-08-24).

## 4. Blocking in tournaments versus the knight torus (FINITE-EXACT + LEAN-CHECKED)

For the Paley tournament `P_7`, Hamiltonian-path blocking:

* every 5-arc deletion leaves a Hamiltonian path (20,349 sets);
* the star at a vertex blocks;
* exactly 63 six-arc sets block.

So `β(P_7) = 6 = hall = N − 1`, as the chessboard session proved for every regular tournament (`hall = N − 1`). The 63 minimum sets split as follows. Hall obstruction means: kills every 1-path-cycle factor (HYP-9168).

| type | count | Hall obstruction? |
|---|---|---|
| two sources | 21 | yes: `|S| = 2`, empty in-neighbourhood |
| two sinks | 21 | yes: `|S| = 2`, empty out-neighbourhood |
| three vertices with one common out- (or in-)neighbour | 14 | yes: `|S| = 3`, neighbourhood of size 1 |
| strip a vertex (star) | 7 | **no**: exotic. The isolated vertex is a 1-vertex path and `P_7 − v` has two cyclic triangles, so a 1-path-cycle factor survives; Hamiltonian paths die simply because the isolated vertex has no arcs |

On the 8-regular bipartite knight torus, Ore's count `d ≥ 6(|S| − |T|) + 1` reaches 7 only when `|S| = |T| + 1`, `N(T) = S`, and the vertex set `S ∪ T` has exactly 8 boundary edges.

* The two trivial cases (`T = ∅`, or `T` = a colour class minus a square) both produce the star.
* Any other case needs an 8-edge cut with at least 2 vertices on each side. But the restricted edge connectivity of `G_n` is 14 (every cut separating two disjoint moves has `>= 14` edges) for `n = 5, 6, 7, 8, 10, 12`, by our max-flow computation, which reproduces the audit's.

So for even `n ∈ {6, 8, 10, 12}`, the 7-sets killing every 2-factor are exactly the stars (PROVED, given the cut computation). The census shows the same for Hamiltonian cycles at `n <= 8`. On a 3-regular tournament, by contrast, Hall sets of sizes 1, 2 and 3 all cost exactly `2k = 6`, and stars are exotic.

The knight's "only one way" is a degree-ratio phenomenon (8 against 2), not a 6- or 7-phenomenon. The "only stars" pattern has one more repo precedent: Sumner `n = 5` fails on 7-vertex tournaments only at the regular ones and only for the pure stars (degree 3 against the star's 4; [sumner_t5_petersen_collatz_20261005.md](sumner_t5_petersen_collatz_20261005.md): "exactly the 3 regular 7-tournaments, exactly the two pure stars").

## 5. The S6 manuscript (typed ANALOGY; status unchanged)

The manuscript [*The (3,4,∞) modular family of 2-tori, completed at its three special points, is a complex structure on S6*](https://alpo.ge/s6.pdf) is the file hashed and audited in [CORE-PAPERS-HOPF-S6-2026-08-24.md](../reference/CORE-PAPERS-HOPF-S6-2026-08-24.md). It remains a **MANUSCRIPT CLAIM / UNDER AUDIT**. The repo does not record the Hopf problem as solved.

| knight torus | S6 manuscript | type |
|---|---|---|
| global question: is there a closed tour (2-factor plus connectivity)? | is the compactified torus family `X` homeomorphic to `S^6`? | — |
| index: Ore/Tutte deficiency, a sum of local terms | Euler characteristic: `χ(X) = Σ χ(singular fibres)` because the 2-torus fibres have `χ = 0` | ANALOGY |
| the whole deficiency at one square (the star) | the whole `χ = 2` in the cusp fibre `W` (ledger: `χ(W) = 2` given the manuscript's stated quotient, homology ranks `(1,2,4,2,1)`; bielliptic fibres contribute 0) | ANALOGY |
| "only stars" (no other minimum) | THM-3991: one-Euler-fibre periodic toric cusps have `χ(W) = d·n!`. For an irreducible fibre (`d = 1`), `χ = 2` forces torus dimension `n = 2`; `(n, d) = (1, 2)` is excluded by irreducibility | ANALOGY ("uniqueness of the place where the index lives") |
| 8 moves, blocking `8 − 1` | `S^6 ⊂ Im 𝕆 = R^7 ⊂ 𝕆 = R^8`; `T_x S^6 = C^3` from the 3 Fano lines through `x` (section 3) | exact on the octonion side only |

**Numerology rejected.** The knight's pairwise determinants 3, 4, 5 set beside the manuscript's triangle group `Δ(3,4,∞)` and `p = 12ℓ_0 − 4ℓ_1 − 3ℓ_2` is NUMEROLOGY. Its orbifold orders come from the monodromy matrices `T_1^3 = T_2^4 = 1`, not from any leap.

**What would make it more than analogy.** A combinatorial model of the cusp fibre `W` (a hexagon with opposite sides glued, an `A_2` triangulation of the torus) in which an Euler-type count and a Hall-type count coincide. None is known here.

## 6. Lyapunov, 6 against 7 (FINITE-NUMERICAL; `n = 6` OPEN)

The setting is the symmetric-maximizer conjecture for `L_A(X) = AX + XA^T` ([octonions reflection, section 8](../../07-reflections/octonions-at-the-center-s6-bott-lyapunov-and-class-rank-frontiers-root-20260824.md); Kressner–Vandereycken, arXiv:2608.20875): it is known for `n ≤ 5`, false for `n ≥ 7`, and open at `n = 6`. We maximise `h(A) = log(σ_skew/σ_sym)` with exact gradients.

* **The KV basin.** KV's integer matrix has `h = 1.80e−4`. L-BFGS-B ascent from it stops at `h = 7.31e−3`, which is **not** a local maximum: the audit's ascent reaches `h = 7.759e−3` (ratio `1.00779`, about 43 times KV's margin). The structure persists: two near-zero rows and a lower block pattern 2 → 3 → 2.
* **Stripping a dimension.** Put `μ(A) = λ_min(A^T A + A A^T)/|A|_F^2`. This is zero iff `A` is orthogonally a padded 6×6 matrix `A_6 ⊕ 0`, which is a counterexample iff `A_6` is.
  * The forward direction is the reflection's padding identity (section 8.1).
  * Converse: `||L_A|Sym||^2 >= 2||A||_2^2` (take `X = v v^T` for a top right singular vector `v`), so the cross block, of norm `||A||_2`, never decides.

  Along the branch, maximise `h` subject to `μ ≤ m` (26 continuation steps, `0.15 ≥ m ≥ 0.0032`). Our continuation gives `h_7 ≈ 0.146 · m^1.47`. These are optimizer-dependent lower bounds on the frontier: the audit's feasible points are larger (`3.94e−3` at `m = 0.06`, `8.09e−4` at `0.015`, `1.09e−4` at `0.0032`; log-slope 1.2–1.3), and the ratio `h_6/h_7` drifts from `−1.0` to `−6.8`. What is robust:
  * `h_7 → 0` as `m → 0`;
  * every stripped 6×6 block has `h < 0`;
  * 7 local re-optimisations of each stripped block never exceed `−1.3e−6`.

  So, numerically, this basin's counterexamples need the seventh dimension: by the fit, their margin tends to 0 as it is removed.
* **L-BFGS / structured random searches.**
  * At `n = 6`: 63,190 starts (lower block-triangular layer patterns, sparse diagonals, random relabelings), best `h = -1.2e-12`.
  * **Hostile control:** at `n = 7` the same method also finds nothing (12,845 starts, best `-2.1e-11`), although counterexamples exist. KV report the same failure of L-BFGS: it stays on the equality ridge `g = 0`.
* **Adam, KV's discovery recipe (section 4 of their paper), with positive controls.**
  * Maximise `g = σ_skew - σ_sym` on the unit sphere with tangent-projected gradients.
  * Adam `(0.9, 0.999)`, step 0.02 for 800 iterations, then 0.006 decaying quadratically to 15% (800 more); variants with steps 0.01 and 0.05 and twice the length.

  | order | runs with `g > 0` | best `σ_skew/σ_sym − 1` |
  |---|---|---|
  | 9 | 4 / 9 | `1.9e-6` |
  | 8 | 25 / 100 | `7.6e-3` |
  | 7 | 37 / 1712 | `7.7e-3` |
  | **6** | **0 / 6000** | runs end on the ridge, `~2e-15` |

  A success probability of 0.1% per run at order 6 would give zero hits in 6000 with probability 0.25%.
* **Shape of the counterexamples found.**
  * In 1000 more order-7 runs, 20 counterexamples were collected.
  * The recurrent top basin (gap `~6.2e-3`, reached 7 times) has normalised singular values `(1, .85, .85, .74, .69, .14, .05)`, five large and two small like KV's rank-5 matrix. Weaker counterexamples do not: one has smallest singular value `0.185`.
  * The `5 + 2` profile is a feature of the top basin, not a necessary condition.
  * A rank-constrained Adam (`A = U V^T`) failed its own `n = 7` control (0/200), so its order-6 runs (0/1200) are not counted.

**Verdict.** No 6x6 counterexample. With working positive controls at orders 7, 8 and 9, the order-6 failure rate is real evidence, recorded as **HYP-9212 (the conjecture holds at `n = 6`)**. `n = 6` is still OPEN; a rarer basin cannot be excluded. The decay of the KV-branch margin under stripping and the top-basin singular profile are the new structural data. In this note's motif, removing the seventh coordinate kills the 7x7 phenomenon here, just as removing one square kills every tour. That shared motif is an ANALOGY.

**Next experiment.** Run Adam at order 6 with structured starts near the stripped top basin (`A_7` with `μ` small) and with step-size schedules tuned on the order-7 hit rate. Seek a classification of the order-7 basins (the top basin recurs; how many others?).

## 7. Additions to opus S16's answer to the pasted prompt

### 7.1 The octic field: S16's open thread 4 closed (FINITE-EXACT, PARI)

* `polredabs` gives the small model `x^8 − 2x^7 + 4x^6 − 6x^5 + 3x^4 − 2x^2 − 5` for `K = Q(T) = Q(u)`.
* `Res_u(U(u), s(1 + 2u − u^2) − (6u + 4)) = −1792 P(s)`, agreeing with S16.
* `K` has no proper subfields, so the Galois action is primitive, consistent with S16's `Gal = S8`.
* **Class number 1**, certified unconditionally (`bnfcertify`). The regulator is `73.89902...`, the unit rank is 4 (signature `(2, 3)`), and the torsion is `±1`.
* Splitting:
  * `2 = 𝔭^4` (`f = 2`);
  * 3 has Frobenius type `(6,2)`, 7 has `(7,1)`, 11 has `(4,3,1)`, 19 has `(4,2,1,1)`;
  * 5, 31 and 104053 each ramify with `e = 2` at one prime.

So the optimal 11-square side lives in an S8 octic with trivial class group, while its quadratic resolvent `Q(√−16128215)` has class number 5444 (S16).

### 7.2 "18 and 19" read as the S15 series' eighteenth and nineteenth Collatz notes

S16 read 18 and 19 as the procgen sessions S18 (Artin, `20/19`) and S19 (Mazur). The opus S15 series also has an eighteenth note ([connectivity from rigidity](collatz_connectivity_from_rigidity_20261001.md)) and a nineteenth ([Paley bridge](collatz_paley_bridge_20261001.md)). Against arXiv:2610.06783 (Alman–Vassilevska Williams):

| paper | note 18 / note 19 | type |
|---|---|---|
| every wanted output string has one **private leaf** (order 0, shared with nothing) | rigidity: the infinite backward tree of `n` determines `n` (THM-4523) | ANALOGY |
| leaves of order `d` are shared by `C(L − m + d, d)` outputs, and few in total (`β_d` geometric) | every acyclic finite view of every integer recurs in `Γ_1` (note 18, Corollary 1) | ANALOGY; **lost**: Collatz backward leaves at a fixed depth are disjoint (a functional graph shares nothing backward), and finite-view counts grow (`2^(R+1) 3^R` classes) instead of decaying |
| Lemma 11 splits at order `m/9`: charge low orders per output, count high orders globally | the quantifier exchange (note 18, Prop. 2): below depth `log_3 n` witnesses are shared, beyond it the only witness is `n` itself | ANALOGY ("sharing horizon"); no Collatz bound follows |
| strings over a 10-letter alphabet, Kronecker/Yates encodings | cycle words over `{0,1}`, codes as `×2`-orbits mod `2^L − 1` (note 19) | no correspondence: the recursion is not shift-invariant and the codes are |

**Where the paper's 18 and 19 come from (PROVED from the paper's displays).** Leaves of order `d` number `β_d = C(L, m−d) 9^(L−m+d)`, with ratio `β_d/β_(d−1) = 9(m−d+1)/(L−m+d) ≤ 9m/(L−m+1)`. This is below `1/2` for every `1 <= d <= m`, uniformly in `m`, iff `L − m + 1 > 18m`. So `L = 19m` is the least integer ratio, `19 = 2·9 + 1`, where 9 is the number of outer-product terms `P_ij` of Schönhage's identity. The exponent in `N ≥ D^18` then makes a tile fit, `N ≥ K N_0 = C(19m, m) 3^(18m)`, because `C(19m, m) 3^(18m) ≤ (19e · 3^18)^m ≤ 4^(18m) = D^18`, and makes the encodings affordable (the paper's (6)). Both are presentation constants (the paper says so); they carry no Collatz content.

**Collatz 18/19 ledger (NUMEROLOGY, recorded so it is not re-derived):**

* the `−17` cycle has non-shortcut length 18 (= 11 + 7, `139 = 3^7 − 2^11`);
* its Mersenne clock is `2^18 − 1 = 3^3 · 7 · 19 · 73`, divisible by the Paley clock 7 and by 19, because `ord_19(2) = 18`;
* `19/12` is a convergent of `log_2 3` (the Pythagorean comma).

No mechanism links any of these to the paper.

### 7.3 The one exact map, and why the new 3SUM bound buys nothing (PROVED, elementary)

A Collatz cycle with shortcut word `w` (length `L`, `p` odd steps) exists iff `c_w ≡ 0 (mod 2^L − 3^p)`. Splitting `w = w_1 w_2 w_3` gives

`c_w = 3^(p_2 + p_3) c_(w_1) + 2^(L_1) 3^(p_3) c_(w_2) + 2^(L_1 + L_2) c_(w_3)`.

So exhaustive cycle search at fixed `(L, p)` is literally modular 3SUM on lists of size `n_3 ≈ C(L/3, p/3) ≈ 2^(H L/3)`, with `H = H(p/L)`. Two halves give modular 2SUM on lists of size `n_2 ≈ 2^(H L/2)`, solvable by hashing in `O~(n_2)`. Any 3SUM algorithm of exponent `2 − ε` costs `2^((2 − ε) H L/3)`, which exceeds `2^(H L/2)` whenever `ε < 1/2`. The new exponent (`ε = 0.0008`) therefore gives no gain over meet-in-the-middle.

Real cycle exclusion (Eliahou, Simons–de Weger, Hercher) does not enumerate words at all; it uses Diophantine bounds. **Typed:** an exact MAP whose consequence is negative.

## 8. Verdicts

| claim | status |
|---|---|
| knight torus: `β(n) = 7`, minimum sets exactly the stars, `n = 5, 6, 7, 8` | FINITE-EXACT (C census; CP-SAT cross-check at `n = 5`) |
| same for all `n ≥ 5` | OPEN (HYP-9211) |
| `G_6 ≅ C_4 × Paley(9)`, twins, `|Aut| = 2^18·144` | PROVED + FINITE-EXACT |
| `G_7`: knight = norm-5 coset of `μ_8` (non-squares), nightrider/queen = Paley(49) complement/Paley(49), order-8 rotation `ω` | PROVED + FINITE-EXACT |
| extra linear symmetries only at `n ∈ {5,6,7,8,10}`, order 8 only at 5 and 7 (`n ≤ 26`) | FINITE-EXACT |
| transversality iff `gcd(n, 30) = 1` | PROVED |
| `G_7 = E_7(5)`: spectrum `= {8} ∪ {-Kl_7(c)}` (each ×8), Ramanujan by Weil; Ramanujan knight tori exactly `n ∈ {5,6,7,8,10}` | PROVED (CITED Weil, Kloosterman norm identity) |
| ladder `PSL(2,7) > Borel = Aut(P_7) > torus = Aut(P_7 − 0) = note 19's Frobenius` | PROVED + LEAN-CHECKED (orders 21, 3) |
| octonionic `J` at `e_0` = Fano matching `q ↦ 3q`, `QR → NQR` | PROVED |
| `β(P_7) = 6`, 63 minima: 56 Hall obstructions (three types) and 7 exotic stars | FINITE-EXACT (count and star LEAN-CHECKED; types in Python) |
| knight sentence vs S6 manuscript | ANALOGY (index localisation); manuscript status unchanged |
| Lyapunov `n = 6` | OPEN (HYP-9212: holds); FINITE-NUMERICAL: Adam 0/6000 at `n = 6` against 37/1712, 25/100, 4/9 at `n = 7, 8, 9`; KV-branch margin `-> 0` under stripping (optimizer-dependent rate) |
| octic: polredabs model, `h(K) = 1` certified, unit rank 4 | FINITE-EXACT |
| 3SUM paper vs S15 notes 18/19 | ANALOGY (private leaf / sharing horizon); 18, 19 = presentation constants; cycle-search map exact, gain none |

## 9. Directions

* **D-a (HYP-9211).** Prove `β(n) = 7` for all `n ≥ 5`. Exclusion of non-star 7-sets needs a robustness lemma: an 8-regular, edge-transitive, 8-connected graph stays Hamiltonian after any 6 deletions. Possible routes are Hamilton-connectedness of `G_n − v` or a Pósa-rotation argument. For odd `n`, a Tutte `f`-factor count would replace Ore's.
* **D-b.** Does the twin structure of `G_6` (pairs `x, x + (3,3)`) force a parity law for 6×6 torus tours, in the way THM-4524's antipodal orbit forces oddness?
* **D-c.** Is there a ladder statement one level up? The knight lives on `μ_8`, the norm-one subgroup of the non-split torus `F_49^*` (order 48). The Collatz trivial cycle lives on `<x2>`, the image in `PSL_2(F_7)` of `SL_2`'s split torus, inside the split torus `F_7^* × F_7^*` of `GL_2` (order 36). The cuspidal and principal-series representations of `GL_2(F_7)` are indexed by characters of these two maximal tori (Deligne–Lusztig, CITED). This is a dictionary only; is there any finite statement it predicts?
* **D-d (Lyapunov, HYP-9212).** Classify the order-7 counterexample basins (one recurrent top basin with KV's rank-5 shape, plus weaker ones). Then run seeded order-6 searches from each basin's stripped limit.

## 10. Audit record

**Audit 1** (2026-10-06, blind subagent, own code in the session scratchpad `audit_sixseven/`; none of the author's scripts run).

**Reproduced:**

* THM-4552 (i)–(iv): the CRT isomorphism (nauty canonical forms), 18 twin pairs, `|Aut(G_6)| = 2^18·144`; the `F_49` facts, including `z` square iff `N(z)^3 = 1` on all 48 units; transversality for `n = 5..60`; the `|Aut|/n^2` table; the linear stabilisers to `n = 26`.
* The blocking census, by an independent method (branching on an unhit Hamiltonian cycle, validated against brute force and CP-SAT on 1,350 instances): `n = 5, 6, 7` stars only, `n = 8` all six-sets. There were also about 800k structured and random spot checks at `n = 7, 8`. The Ore count was checked.
* THM-4553 (i)–(iii): the `PSL(2,7)` ladder, the associator alternating on all 512 basis triples, and `J`.
* `H(P_7) = 189`, the count 63, and the Lean file (no `sorryAx`).
* The octic data, including `h(Q(√−16128215)) = 5444`.
* The 18/19 derivation (`L = 19m` is the least `L` for every `m <= 199`), the numerology ledger, and the cycle-composition formula (at `(L, p) = (11, 7)` the search returns exactly the `−17` cycle).
* The KV certificate.

**Corrections applied (MISTAKE-567):**

* The 7 stars of `P_7` are exotic minima, not Hall obstructions. The split is 56 + 7.
* The KV-basin point was not a local maximum: ascent reaches `7.759e−3`.
* The stripping "law" was retyped as optimizer-dependent lower bounds.
* "8 → 7" is NUMEROLOGY.
* The 4+4 slope split is not `PSL(2,7)`-invariant.
* The Deligne–Lusztig tori were renamed correctly.
* `χ = 2` forces `n = 2` only for irreducible fibres.
* Credit for the antipodal orbit goes to THM-4524.
* Smaller wording fixes: affine plane of order 6; 6 non-trivial translations; torus symbol `H`; "every `d <= m`"; 104053 added to the ramified primes.

**Items the audit flagged that were already fixed before it finished:**

* the `n = 7` row;
* the `n = 8` census, completed with 14 stars among `C(255,6)` seven-sets and recorded in the `.out`;
* "order 8 only at 7" changed to "among `n >= 6`".

The audit also supplied the restricted-edge-connectivity value 14, which this session re-derived after fixing a multi-edge bug in its own first attempt (13 by mistake), and the padding converse.

**Audit 2** (the late additions: section 2.3 Kloosterman/Ramanujan, the modular-curve dictionary, the Ore/λ' paragraph, the Adam table and HYP-9212): *(in progress)*

## 11. Reproduction

```bash
python3 04-computation/experiments/sixseven_20261006_structure.py          # 10 s, needs dreadnaut
cc -O2 -DN=7 -o tb7 04-computation/experiments/sixseven_20261006_knight_torusblock_n.c && ./tb7 6 1000 2 && ./tb7 7 1000 3
python3 04-computation/experiments/sixseven_20261006_knight_torus_block.py 5   # CP-SAT cross-check, 90 s
python3 04-computation/experiments/sixseven_20261006_lyapunov_six.py kv; python3 04-computation/experiments/sixseven_20261006_lyapunov_six.py strip
gp -q 04-computation/experiments/sixseven_20261006_octic_class_group.gp < /dev/null
lean 04-computation/lean/standalone/sixseven_20261006_paley_seven_certificate.lean   # 60 s
```
