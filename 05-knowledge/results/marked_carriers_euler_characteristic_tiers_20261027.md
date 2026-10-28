# Tournaments and graphs with one fewer edge than vertices: the Euler-characteristic trijection, its tiers, and procedural counterexample families

kind-pasteur-2026-10-27-S8. This builds on the owner's trijection (kind-pasteur S3: a tournament with a marked line,
cycle+solo or star is a tournament on one vertex fewer modulo 1, `Z_{n-1}`, `S_{n-1}`) and on the symmetry calculus
of session80 (composition, cycle indices, explicit collisions). It was written concurrently with session81, which read
the same prompt as an encoding-capacity question (credited in §6).

- **Code and data:** `04-computation/tournament_carriers/`.
- **Pre-registration:** `PREREGISTRATION.md` (P1–P3), pushed in efdf0b0a5f before the confirmatory runs.
  - P1 PASSED.
  - P2 FAILED and is recorded.
  - P3 PASSED.
  - The census that suggested the main criterion was exploratory.
- **Canon:** THM-4835, THM-4836, HYP-9435.

## 1. Inheritance, board, niche, wildcard

- **Closest proved mechanism.** S3 Theorem 1 (three marked bijections) and LEM-003 (automorphisms act freely on
  Hamiltonian paths).
- **Hostile.** Session80's collision: equal counts, different groups. Here, the bigon `C_4 + K_1` has the same
  marked count as `Tour(4)/Z_2`, yet admits no pair-level trijection.
- **Corrected near miss.** S3's numeral `8 = 4 + 4` (lines = cycles + stars) holds only at n = 4 (MISTAKE-529 style).
- **Least-used sidecar.** The Burnside ring / table of marks (tom Dieck's equivariant Euler characteristic).
- **Niche.** Non-functional orientations of disconnected excess -1 graphs.
- **Wildcard.** Maxwell rigidity counts `jn - C(j+1,2)`.
- **Board.** Marked carrier and Aut(D); signed pair sets; marks and `χ(G^K)`; linear vs weak vs count tiers; coherent
  carrier sets and session80's composition calculus; the Maxwell ladder.
- **META-PATTERNS cards used.**
  - "The same representation is not the same carrier": a permutation character does not decide a signed A-set.
  - "Controlled forgetting and unlabeled quotients require a sidecar": the tiers are exactly the sidecars lost by
    counting.

## 2. The answer in one paragraph

The owner's "graphs with one fewer edge than vertices" are the Euler-characteristic-1 graphs: trees, and disconnected
graphs whose extra cycles are paid for by extra components. For every such graph the complement has exactly
`C(n-1,2)` pairs, the size of `K_{n-1}`. So a tournament with a marked oriented copy D is, as a *set*, a tournament on
`n-1` vertices.

The trijection is *equivariant and linear*, with the reduced vertex set equal to the edge set of the carrier, exactly
when every symmetry of D fixes a sub-carrier of Euler characteristic 1. Equivalently, `V ≅ E ⊔ {pt}` as Aut(D)-sets
(THM-4835).
- Trees (any orientation) and functional orientations always qualify, with explicit maps. The S3 trijection and the
  vertex-deleted deck are three instances.
- Symmetric bigons, alternating even cycles, dense cores with isolates and odd thetas do not.

The obstruction has two tiers:
- involutions (the 2-part) break only the linear structure and are invisible to counts;
- odd-order symmetries break counts as well.

Read as "each tournament maps to a pair or a triple" of carrier kinds, the statement is exact precisely for `n <= 5`
(THM-4836). The first singleton is `C_3 ⇒ C_3`.

## 3. Results

### 3.1 The criterion and its proof (THM-4835 (2), (4))

- **The mark argument.** View the ordered non-edges `Ω` as an `(A x Z2)`-set (Z2 = swap). A linear trijection is an
  `(A x Z2)`-isomorphism with the ordered pairs of an `(n-1)`-set X. Burnside's marks reduce this to two families of
  equations over subgroups K:
  - `C(|V^K|,2) - |E^K| = C(|X^K|,2)`;
  - `|V^{K0}| - |V^K| = |X^{K0}| - |X^K|` for `K0` of index 2 in K.
- **Sufficiency.** With `X = V \ {r}`, the equations are exactly `|E^K| = |V^K| - 1`.
- **Necessity.**
  - At 2-subgroups: the second family makes `|V^K| - |X^K|` constant (= 1) along index-2 chains.
  - At odd prime order (full A): moved components must have `χ = 0` (else a transposition creates an involution
    defect). Relative homology gives `δ >= 0`. The mark arithmetic `δ = pa(2v-3+pa)/2` then forces `a = 0`.

### 3.2 Explicit trijections (THM-4835 (3))

- **Trees.** `Tour(V) ≅ Or(G) x Tour(E(G))` for every spanning tree G, `Aut(G)`-equivariantly. The map sends a
  non-edge `{u,w}` to the end edges of the u-w path, with `u -> w` sent to `e_u -> e_w`.
  - For the path this is S3's staircase.
  - For the star it is `T -> T - v`.
- **Functional orientations.** The root-swap `{v, f(v)} -> {r, v}`, with `v -> f(v)` corresponding to `r -> v`, is an
  equivariant linear trijection onto `Tour(V \ {r})`.
  - For the directed cycle plus solo it is S3's map.
- **Positive controls.** 600 random trees and functional orientations at n = 8 and n = 9: all bijective, all
  equivariant (P1(d)).

### 3.3 The tiers (THM-4835 (5)–(8))

| tier | sees | broken by | witnesses |
|---|---|---|---|
| counts `\|{T ⊇ D}/A\|` | odd-order subgroups (odd-order rule) | odd-order defects | `θ_p^k + (p-1)K_1` (all odd primes p, `k >= 2`); `K_4(C_3+apex) + 3K_1`; `K_{3,2} + 2K_1` |
| A-set `{T ⊇ D}` (weak) | odd-order subgroups | same | same; equal to the count tier for n <= 7 |
| signed A-set (linear) | all subgroups | involution defects (always even) | bigons `B_k`, alternating `C_{2k} + K_1`, dense cores `+ (j+1)K_1` |

Census, n = 3..7, 1250 carriers: linear 1035 (including the two n = 7 stars), weak only 210, neither 5.
- **Every non-linear carrier has an involution defect.** Conjectured not to persist (HYP-9435(3)).
- **Group dependence.** `K_4 + 3K_1` is linear for `Z_3` but not for the full `S_3`.

### 3.4 Pairs and triples (THM-4836)

- **The three classes.** `Σ(T)` records L (always, Rédei), C (strong, or strong `(n-1)` plus source/sink; Camion–Moon)
  and S (source/sink).
- **Class sizes.**

  | `Σ(T)` | types |
  |---|---|
  | `{L,C}` | `s_n` |
  | `{L,C,S}` | `2s_{n-1}` |
  | `{L,S}` | `2t_{n-1} - t_{n-2} - 2s_{n-1}` |
  | `{L}` | `[x^n] S_{>=3}^2 T` |

- **Singleton counts.** 1, 3, 16 at n = 6, 7, 8.
- **The pair-or-triple law** holds iff `n <= 5`.
- **Failures of the law** are exactly the two-sided strong paddings of session80's ordered-sum calculus.

### 3.5 The Maxwell ladder (HYP-9435(2), pre-registered P2 FAILED)

- **Setup.** Marking `jn - C(j+1,2)` edges reduces by j vertices ("pairs" j = 2, "triples" j = 3).
- **Proved.**
  - The count `|E^K| = j|V^K| - C(j+1,2)` with `|V^A| >= j` is sufficient.
  - Its 2-subgroup instances are necessary.
- **The pre-registered equivalence failed** (j = 2, n = 6: 6/3817; j = 3: 14/2790), always at order-3 subgroups.
  So j = 1, the trijection, is the rigid rung.

### 3.6 A free identity (THM-4835 corollary)

`Σ_T Σ_r ∏_{v≠r} s_v(T) = F_n 2^{C(n-1,2)}`, where `F_n` = 9, 72, 815, 11700 counts rooted functional digraphs
without loops or 2-cycles. This is checked by brute force for n = 3..6.

## 4. Procedural generators of counterexamples

Each generator comes with a proof, not only a search.

1. **Involution gadgets (break linearity only).** Insert a symmetric bigon, an alternating even cycle, or two
   swappable isolates next to a dense core. Any involution whose fixed sub-carrier has `χ = 1 ± 2` kills the linear
   trijection. Counts and A-sets survive whenever odd-order defects vanish (the THM-4835 (5) lemma).
2. **Odd gadgets (break everything).** p identical directed u-w paths plus `p-1` isolates. The `Z_p` mark forces a
   fixed count `p+1 ≢ 0 (mod p)`, so the A-set, and here even the count, has no `Tour(X)` model.
3. **Chirality (predicted).** Remove all involutions by chiral decorations so that only odd defects remain. This
   would make non-linearity without an involution witness (HYP-9435(3)).

## 5. Unifying proposal (typed)

**PROVED.**
- Marking an oriented excess -1 carrier is the same thing as passing to the edge set. The price is one vertex, paid
  equivariantly iff the carrier's Burnside Euler characteristic is a point.
- This one statement covers:
  - S3's trijection;
  - the deck (T - v);
  - Hamiltonian-path staircases;
  - functional digraphs (choice functions, `∏ s_v`);
  - arborescences (spanning in-trees).

**PROVED (tier law).**
- The 2-part of a symmetry group is seen only by the signed (linear) structure.
- The odd part is seen by counts.
- This sharpens the classical odd-order rule (a permutation fixes a tournament iff it has odd order) into an
  obstruction theory.

**ANALOGY (typed, no implication).**
- *The Maxwell ladder:* marking rigidity-count graphs removes j vertices.
- *Session80's information hierarchy:* cycle index, groups, subgroup marks.
- *Session81's "the representation must be judged against its future operations":* our tiers are an instance, with
  the operations being the subgroup restrictions.

## 6. Concurrent and prior work

- **Session81** (`session81_counterexample_families_20261028`, first-pushed today) reads the same prompt as a
  capacity problem.
  - Independent pairs/triples of `(n-1)`-edge graph types cannot encode tournament types (failures at n = 3, 9, 14).
  - One padded sign-even carrier can.
  - It also finds `K_4 + 3K_1` as the first non-forest `E = V-1` odd-cycle conflict graph.

  That is a different role from our sporadic non-weak carrier `K_4(C_3+apex) + 3K_1`; the coincidence of the graph is
  noted, not identified. Our results are about equivariant marked trijections, which session81 lists as needing the
  mark's symmetry and orientation as "indispensable sidecars"; THM-4835 says exactly which sidecar is needed.
- **Royle–Praeger–Glasby–Freedman–Devillers** (arXiv:2204.01947): tournaments and sign-even graphs are equinumerous.
  Their "sign" (parity of edges reversed by an automorphism) is the global form of our involution reversals. No
  implication is claimed.
- **Classical inputs.** Rédei; Camion; Moon (pancyclicity); Burnside's table of marks; fixed subtrees of finite group
  actions on trees.

## 7. Reproduction

From `04-computation/tournament_carriers/`:

| command | covers |
|---|---|
| `python carrier_census.py 7` | census, about 20 min |
| `python families_and_controls.py` | P1 |
| `python maxwell_ladder.py` | P2 |
| `python coherent_carrier_sets.py` | P3 |
| `python tiers_and_sporadics.py` | tiers, `θ_p`, group dependence, identity |
