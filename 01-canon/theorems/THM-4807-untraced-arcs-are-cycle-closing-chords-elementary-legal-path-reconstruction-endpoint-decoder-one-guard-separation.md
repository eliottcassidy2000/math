---
id: THM-4807
title: "Untraced arcs are cycle-closing chords. Let T be a strong tournament and u -> v an arc lying on no Hamilton path (untraced). (a) T - u is not strong and v lies in the terminal strong component of T - u; dually T - v is not strong and u lies in its initial component. (b) u -> v is the chord of a Hamilton path from u to v, all of whose arcs are traced, and no Hamilton path runs from v to u; reversing u -> v, together with any further untraced arcs, closes that path into a Hamiltonian cycle. (c) Exact tail-deleted characterisation: u -> v is untraced iff T - u has ordered components C_1, ..., C_k with k >= 2, v in C_k, u dominates C_(k-1), and C_k has no partition into a path ending at an in-neighbour of u and a path starting at v. (d) One-guard separation and reconstruction: for labelled tournaments T != T' there is an explicit vertex word that is a Hamilton path of exactly one of them, so the set of Hamilton paths determines every finite labelled tournament, by an elementary proof that needs only Camion's theorem. (e) Endpoint decoder: inside a strong block, an untraced pair {x, y} is oriented x -> y iff some legal word restricted to the block starts at x and ends at y."
status: >
  PROVED, self-contained and elementary; the only external inputs are Redei's and Camion's theorems. The reconstruction
  consequence (d) was proved concurrently, the same day, relative to the CITED Busch-Jacobson-Reid structure theorem
  (T - z has at least four components with the untraced arc from the first to the last) by
  guard_reconstruction_20261026 (G1-G3, two-step-count decoder), by Session70 (session70-guards G1-G5,
  deleted-vertex reachability decoder) and by growing_sieve_path_20261026 (paths P1-P4, deleted-vertex reachability
  on the used-arc union; merged in sieve_path_merge_20261026). This theorem removes that dependency and adds (b),
  (c), (e) and the explicit separating word.
  Literature: (a)'s first half is Busch's Proposition 4.1 (thesis, 2005): untraced arcs join two cut vertices.
  (c) is the tail-deleted analogue of Busch's Lemma 4.13, which is stated with a separating vertex z rather than the
  tail u. Busch's Corollary 4.8 (an arc x -> y with d+(y) > d+(x) is traced) gives d+(u) >= d+(v) for every untraced
  arc. No priority is claimed for the elementary statements. No reconstruction statement was found in the thesis.
  FINITE-EXACT: every tournament class for n <= 9, labelled n <= 6, ordered labelled pairs n <= 5, all three decoders
  compared n <= 8; stress tests with planted cut vertices to n = 18.
  The same-day S61 THM-4806 (R) proves the single-flip case (i) from Camion and the strong-component reduction (ii),
  and leaves all orders conditional on a "strong form" (HYP-9417 (1b)): every untraced arc u -> v of a strong C
  has a Hamilton path of C^(uv) through v -> u avoiding every other untraced arc. (b) below proves that strong form:
  the path P'u uses only traced arcs and v -> u.
  AUDITED (independent, 2026-10-26): (a)-(e) CONFIRMED at every order; one major error (the traced core was said not
  to suffice, which is false) and twelve minors were corrected before the first push.
  Report: 04-computation/experiments/hamilton_path_recovery_20261026/audit_thm4807/REPORT.md. The audit extends the
  census to every class at n = 9 and every strong class at n = 10.
session: mac-mini-2026-10-25-pyth (2026-10-26)
scripts:
  - 04-computation/experiments/hamilton_path_recovery_20261026/verify_reconstruction.py (+ .out, and _n9.py / _n9.out): (a), (b), (c) and the decoder (e) on every class n <= 9 (191,536 classes at 9; 200 untraced arcs in strong classes; (c) exact on 6,412,788 arcs), labelled n <= 6
  - 04-computation/experiments/hamilton_path_recovery_20261026/separating_word.py (+ .out): the explicit word of (d) on all 1,047,552 ordered labelled pairs at n = 5 and on biased random pairs n = 6..9 (cases alpha / beta / gamma all exercised; 0 failures)
  - 04-computation/experiments/hamilton_path_recovery_20261026/three_decoders.py (+ .out): endpoint (this), two-step count (G1) and deleted-reachability (Session70 G4) decoders agree and are exact on every class n <= 8
  - 04-computation/experiments/hamilton_path_recovery_20261026/stress_cut_vertex.py (+ .out): bitmask-DP decoder on 1000 relabelled tournaments with planted (and nested) cut vertices, n <= 18: 3363 untraced arcs, 945 random flip sets, 0 failures
  - 04-computation/experiments/hamilton_path_recovery_20261026/ham_recovery.py, missed_arcs.py, injectivity_flips.py (+ .out): the first census (cut positions, skipping separators, flip sets n <= 8)
related:
  - 05-knowledge/results/guard_reconstruction_20261026/README.md                       # G1-G3, concurrent, relative to BJR
  - 05-knowledge/results/session70_attachment_guard_alignment_20261026/session70-guards/REPORT.md   # G1-G5 (G6-G7 on elementary relocations), concurrent, relative to BJR
  - 05-knowledge/results/growing_sieve_path_20261026/paths/RESULT.md                   # P1-P4, concurrent, relative to Busch Thm 4.12
  - 05-knowledge/results/sieve_path_merge_20261026/README.md                           # integrates the three concurrent proofs
  - 05-knowledge/results/strong_atom_continuation_vitali_20261026/README.md            # C5: the question
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md                          # step (1): the same ordered-join factorisation
  - 01-canon/theorems/THM-4806-named-legality-language-single-pair-lemma-goldbach-free-shapes-and-row-companions-of-even-perfect-numbers.md   # (R): single flip, reduction, strong form (proved here by (b))
  - 05-knowledge/hypotheses/HYP-9418-legal-path-sub-observers-berggren-sieve-levels-and-the-branch-character.md
---

# THM-4807 — Untraced arcs are cycle-closing chords

## Conventions

- Tournaments are on a named finite vertex set.
- `HP(T)` is the set of Hamilton paths, written as vertex words.
- An arc is **traced** if it joins consecutive letters of some word in `HP(T)`, and **untraced** otherwise.
- `T − x` is the subtournament on the other vertices.
- The strong components of any tournament are linearly ordered, `C_1 ⇒ C_2 ⇒ ⋯ ⇒ C_k`, with every arc between different
  components pointing forward.
- Two facts are used throughout:
  - Rédei: every tournament has a Hamilton path.
  - Camion: a strong tournament on at least 3 vertices has a Hamiltonian cycle.
- Consequently every vertex of a strong component `C` starts some Hamilton path of `T[C]`, and ends some Hamilton path of
  `T[C]`. For a singleton this is trivial. Otherwise `|C| ≥ 3`, since no 2-vertex tournament is strong, so `T[C]` has
  a Hamiltonian cycle; cut it.
- **The ordered-join factorisation:** every Hamilton path of `T` is a concatenation `Q_1 Q_2 ⋯ Q_k` of Hamilton paths
  `Q_i` of the components, and every such concatenation is a Hamilton path. A path never returns to a component it
  has left, since every cross arc points forward, so to cover every vertex it runs through `C_1, …, C_k` contiguously
  and in order. Conversely, the last vertex of `Q_i` beats the first vertex of `Q_(i+1)`. Hence:
  - an arc inside a component is traced in `T` iff it is traced in that component;
  - every arc between consecutive components is traced;
  - no arc between non-consecutive components is traced.

## Statement and proofs

### (a) Where an untraced arc can sit

Let `T` be strong and let `u → v` be untraced. Then `T − u` is not strong, and `v` lies in the terminal component of
`T − u`. Applying this to the converse tournament: `T − v` is not strong, and `u` lies in its initial component.

*Proof.*
- **`T − u` strong (so `n ≥ 4`).** `T − u` has a Hamiltonian cycle. Cut it so that it starts at `v`, and prepend `u`.
- **`T − u` not strong.** Let its components be `C_1, …, C_k` with `v ∈ C_j`. Strongness of `T` gives two facts:
  - nothing can enter `C_1` except from `u`, so `u` has an out-neighbour in `C_1`;
  - `C_k` can leave only to `u`, so some `a ∈ C_k` has `a → u`.
- **`j = 1`.** Then `u Q_1 Q_2 ⋯ Q_k`, with `Q_1` starting at `v`, is a Hamilton path through `u → v`.
- **`2 ≤ j ≤ k − 1`.** Then

      Q_1 ⋯ Q_(j−1) Q_(j+1) ⋯ Q_k  u  Q_j,      Q_k ending at a,  Q_j starting at v,

  is a Hamilton path through `u → v`. The jump `Q_(j−1) → Q_(j+1)` points forward, and so does `a → u`.

So `j = k`. ∎

The first half, that `u` and `v` are cut vertices, is Busch's Proposition 4.1. The middle-component detour is what
excludes `2 ≤ j ≤ k − 1`.

### (b) The chord, and the cycle it closes

Keep the notation of (a), with `v ∈ C_k`.
- Choose `c ∈ C_1` with `u → c`. Let `Q_1` start at `c`, and `Q_k` end at `v`.
- Then `P' = Q_1 ⋯ Q_k` is a Hamilton path of `T − u` from `c` to `v`.
- So `uP'` is a Hamilton path of `T` from `u` to `v`. All its arcs are traced, since it is itself a word of `HP(T)`.
- `|C_k| ≥ 3`, because `C_k = {v}` would force `v → u`.
- No Hamilton path runs from `v` to `u`. Such a path, with its first letter `v` removed, would be a Hamilton path of
  `T − v` ending at `u`. But it visits the components of `T − v` in order (at least two of them, by the dual of (a)),
  so it ends outside the initial one, which contains `u`.

If `D` is any set of untraced arcs containing `u → v`, reverse `D` to get `T_D`.
- The arcs of `uP'` are traced, so they are untouched, and `v → u ∈ T_D`.
- Hence `uP'u` is a Hamiltonian cycle of `T_D`.
- Cutting it just after `u` gives the word `P'u`. That word is a Hamilton path of `T_D` and not of `T`. ∎

**Mechanism.** The arc hidden from every legal word is the shortcut of a legal word from its tail to its head. Reversing
`u → v` alone turns every Hamilton path of `T` from `u` to `v` into a new Hamiltonian cycle through `v → u`. These are
all the new Hamiltonian cycles, and there can be several: 1, 2, 3 or 5 of them for `n ≤ 8`, and up to 45 at `n = 10`.
Shorter new cycles appear too, for instance a 3-cycle `v → u → w → v` for every `w` with `u → w → v` (2 to 5 of them
for `n ≤ 8`).

### (c) Exact tail-deleted characterisation

Let `T` be strong and `u → v` an arc. Then `u → v` is untraced iff all of the following hold:
- `T − u` has components `C_1, …, C_k` with `k ≥ 2` and `v ∈ C_k`;
- `u` dominates every vertex of `C_(k−1)`;
- `C_k` admits no partition into a nonempty path `α` ending at an in-neighbour of `u` and a path `β` starting at `v`.

*Proof.*

**Necessity.** The first condition is (a). Suppose `T − u` is not strong with `v ∈ C_k`, and `x → u` for some
`x ∈ C_(k−1)`. Then `Q_1 ⋯ Q_(k−1) u Q_k`, with `Q_(k−1)` ending at `x` and `Q_k` starting at `v`, is a Hamilton path
through `u → v`. A partition `(α, β)` of `C_k` gives the path `Q_1 ⋯ Q_(k−1) α u β` through `u → v`.

**Sufficiency.** Let `A u v B` be a Hamilton path. The path `vB` lies in `T − u` and starts in `C_k`, so it stays in
`C_k`. The path `A` therefore covers `C_1 ∪ ⋯ ∪ C_(k−1)`, which forces `A` to be nonempty. Its last vertex beats `u`.
Because `A` is a path of `T − u` and every arc between components points forward, `A` meets the components in
nondecreasing order. So that last vertex lies in the last component `A` meets, which is `C_(k−1)` (dominated by `u`)
or `C_k`. In the second case the part of `A` inside `C_k` is a final segment of `A`: a nonempty path ending at an
in-neighbour of `u`. With `vB` it forms a forbidden partition. ∎

This characterisation is exact on every arc of every strong class with `n ≤ 9`: 6,412,788 arcs at `n = 9` alone
(178,133 classes, 36 arcs each), and 6,589,019 arcs for `n ≤ 9` in total. The independent audit re-checked it on all
421,017,705 arcs of the strong classes at `n = 10`. None of the three conditions is redundant: at `n = 8`, 35 arcs
fail only the second and 6,308 fail only the third. Busch's Lemma 4.13 states the same type of condition for the
separating vertex `z` of his Theorem 4.12, for which `T − z` has four or more strong components. Here the deleted
vertex is the tail itself.

### (d) One explicit separating guard; reconstruction

Let `T ≠ T'` be tournaments on the same vertex set. Exactly one of the following cases holds, and each gives a word in
`HP(T) △ HP(T')`.

- **(α)** Some traced arc of `T` is reversed in `T'`. Then a word of `HP(T)` through it is legal in `T` and illegal in
  `T'`.
- **(β)** Otherwise, suppose some untraced arc `u → v` inside a strong component `B` of `T` is reversed.
  - In `T[B]`, the arc `u → v` is untraced, so (b) supplies `P'_B` with `P'_B u` legal in `T'[B]`. The arcs of
    `P'_B` are traced in `T` and hence unchanged.
  - Concatenate Hamilton paths of `T'[C]` for the other components `C` (Rédei), in component order, with `P'_B u` at
    `B`'s place.
  - The arcs between consecutive components are traced and hence unchanged.
  - The resulting word is legal in `T'` and illegal in `T`.
- **(γ)** Otherwise, only arcs between non-consecutive components differ.
  - Let `i` be the least component index of a tail of a differing arc. For that `i`, let `j` be the largest head index.
    Take a differing arc `u → v` with `u ∈ C_i` and `v ∈ C_j`, so `j ≥ i + 2`.
  - The word is

        Q_1 ⋯ Q_(i−1)  Q_(i+1) ⋯ Q_(j−1)  [Q_j ending at v]  [Q_i starting at u]  Q_(j+1) ⋯ Q_k,

    with `Q_t` Hamilton paths of the unchanged `T[C_t]`.
  - It uses `v → u` together with at most two possibly-skipping arcs (the first is absent when `i = 1`, the second
    when `j = k`):
    - `C_(i−1) → C_(i+1)`: its tail index is below `i`, so it is not a differing arc;
    - `C_i → C_(j+1)`: its head index is above `j`, so it is not a differing arc.
  - So the word is legal in `T'` and illegal in `T`.
  - Both extremal choices are needed. Taking the smallest head index instead fails for 3,120 of the 8,160 case-(γ)
    ordered pairs at `n = 5`, and taking the largest tail index fails for 3,360.

Hence `T ↦ HP(T)` is injective on labelled tournaments of every finite order. In the language of vitali-response C5, the
complete named legal re-marking guard signature separates all of `X_n`. The minimum operation-stable quotient for that
protocol has `2^((n−1)(n−2)/2)` states, at every `n`. These are valid presentation states, as in C5; an absorbing
failure state, if counted, adds one (guard_reconstruction G2). ∎

The separating word is constructive and has one letter per vertex. It was checked on all 1,047,552 ordered labelled
pairs at `n = 5`.

### (e) The endpoint decoder

From `𝓗 = HP(T)` alone:
1. Keep every traced arc.
2. Find the **cut positions**: the positions `m` at which all words of `𝓗` share the same prefix set. They are exactly
   the boundaries between strong components, because every nonsingleton component rotates its Hamiltonian cycle freely
   inside the factorisation.
3. Orient every pair between different blocks forward.
4. Inside a block `B`, orient an untraced pair `{x, y}` as `x → y` iff some word of `𝓗`, restricted to `B`, starts at
   `x` and ends at `y`.

By (b), applied to `T[B]` (whose Hamilton paths are exactly the restrictions), exactly one direction qualifies.

The decoder reads three things, all from `𝓗`: the traced arcs, the cut positions, and the block endpoint pairs. The
two halves differ:
- the traced core `D` (the set of traced arcs) does not display the orientation of an untraced pair, but it loses no
  information: `D ⊆ T` and every word of `𝓗` uses only traced arcs, so `HP(D) = HP(T)` (Session70 G1). By (d), the
  traced core alone therefore determines `T`. Reading the untraced orientations off `D` needs a non-local rule: this
  decoder applied to `HP(D)`, the two-step count of G1, or the deleted-vertex reachability of Session70 G4;
- the endpoint relation `{(first letter, last letter)}` alone does not determine `T`: it collides from `n = 4`, and
  never for `n ≤ 3`. There are 44 distinct endpoint relations among 64 labelled 4-tournaments, 614 among 1,024 at
  `n = 5`, and 15,343 among 32,768 at `n = 6`.

On every class `n ≤ 8`, three decoders return `T` and agree:
- this one;
- the two-step count of G1, which orients `u → v` when it has at least two traced two-step routes;
- the one-way deleted-vertex reachability of Session70 G4.

At untraced pairs, the two-step counts `(max, min)` take only the values `(2,1), (3,0), (3,1), (4,0), (4,1), (5,1)`
through `n = 8`, and the maximum is always the true direction. The minimum is never above 1, at every order, but
this needs more than (a)–(e). Two routes `v → w → u` would be two internally disjoint paths from `v` to `u`, and
Busch's Theorem 4.7 would then put `u → v` on a Hamilton path. Equivalently, BJR's separator `z` lies on every path
from `v` to `u`; this is G1's argument. The audit found the minimum at most 1 on every untraced pair of every class
with `n ≤ 9`, and of every strong class with `n = 10`.

## What is and is not new

- **Concurrent proofs.** The reconstruction (d) itself is proved concurrently, relative to the Busch–Jacobson–Reid
  four-block structure theorem (Busch's thesis Theorem 4.12), in three places: guard_reconstruction_20261026 G1,
  Session70's guard report G4/G5, and growing_sieve_path_20261026 paths P1–P4, integrated with the other two in
  sieve_path_merge_20261026. All three recover `T` from the traced core alone. Each already notes that one named
  guard separates two distinct states, which is immediate from injectivity. What is new here is the explicit
  construction of that guard.
- **Contribution here:**
  - an elementary proof that does not need that theorem: (a) and (b) need only Camion, and the four-block structure is
    never used;
  - the chord and cycle-closing mechanism (b), which is exactly the "strong form" left open in S61's HYP-9417 (1b),
    so S61's single-flip route also closes at every order;
  - the exact tail form (c);
  - the explicit single separating word (d);
  - the endpoint decoder (e).
- **Cut vertices.** The vertex cut named in the owner's prompt is the tail `u` itself. In the BJR form it is a separate
  vertex `z`. Both are visible in `𝓗`: `u` and `v` are cut vertices, and `𝓗` contains a word from `u` to `v` and none
  from `v` to `u`.
- **Connection to 7 and 21.** The same ordered-join factorisation gives the multiplicativity `H(T) = ∏ H(C_i)`, which
  is step (1) of THM-1370 (H never equals 7 or 21). That step reduces 7 and 21 to strong components; the exclusion
  itself comes from the strong floor (Moon's theorem and the exhaustive small spectra). Reconstruction reads the
  factorisation off the cut positions. Untraced arcs occur only where a strong block is "one vertex away" from such
  a factorisation: by (a), both ends are cut vertices. The converse fails. The strong 4-tournament has cut vertices
  but no untraced arc, and 3,600 of the 3,633 strong 8-classes with a cut vertex have no untraced arc. HYP-9418 (1)
  records the quantitative version that remains open.

## Open

See HYP-9418 (1):
- the smallest sub-observers of `𝓗` that still reconstruct;
- the number of untraced arcs (`1, 2, 11, 39, 200` in strong classes, `n = 5..9`; `1,207` at `n = 10`, audit). The
  maximum per strong `n`-tournament is already known: at most `(n² − 4n + 3)/8`, attained at odd `n` by Busch's
  `T_max` (Busch's thesis, Theorem 4.26; Busch–Jacobson, JCMCC 63 (2007)). The census maxima are 1, 1, 3, 3, 6, 6
  for `n = 5..10`;
- the cut-vertex form of the H-floor.
