# Lean formalization of THM-4524, THM-4525, THM-4526 and S15's Collatz drop results (procgen lean lane, 2026-10-01)

**Lane:** procgen lean lane (resumed instance; the first instance's files were lost in a reboot), session
`collatz-procgen-20260922`, 2026-10-01, machine mac-mini (Apple M2, 8 GB).
**Package:** `04-computation/lean/ProcgenSelfieEdim/` (Lean 4 core `v4.30.0` with the bundled `Std`; no Mathlib, no downloads).
**Verifier:** `python3 verify.py` in the package: PASS, **503 theorems** in 45 modules, every one within
`{propext, Quot.sound}` (152 use no axiom, 54 only `propext`, 1 only `Quot.sound`, 296 `propext, Quot.sound`).
**No theorem uses `Classical.choice`.** No `sorry`, no `native_decide` / `decide +native`, no custom axiom.
**Status:** FORMALIZED (Lean-checked) for the statements in sections 3-6e. Independent audit: OWED (orchestrator).

## 0. Status table

| Target (orchestrator's list) | Lean module | Status |
|---|---|---|
| 1(a) THM-4525: `Q_d`, edges, edge-vertex distance, histograms, edge-multiset resolving; the paper's 15-set resolves `Q_6` | `Hypercube`, `EdimQ6` | FORMALIZED (`paperSet_resolving`) |
| 1(b) L2 antipodal reversal, all `d` | `Hypercube` | FORMALIZED (`hist_antipode`; plus the parity corollaries `hist_not_palindrome`, `no_antipodal_collision`) |
| 1(c) L4: no set of size `<= 6` resolves `Q_6` (counting + pigeonhole) | `CountingBound` | FORMALIZED (`not_resolving_of_length_le_six`) |
| 1(d) L1: trivial stabilizer | `Stabilizer` | FORMALIZED (`trivial_stabilizer`, all `d >= 2`) |
| 2(a) THM-4524 A1, gauge theorem, all `N` | `Tournament` | FORMALIZED (`path_arc`, `selfie_fibre`, `loops_ne_compl`, `selfie_onto`) |
| 2(b) `QR_7 - v` all-odd (H = 45; 13 and 23); `QR_7` every arc 54; no all-odd 5-tournament | `SelfieFinite` | FORMALIZED (`qr7del_counts` for all 7 vertices, `qr7_counts`, `no_allOdd_three_to_five` for N = 3, 4, 5) |
| 2(c) HP infrastructure, `sum_e c(e) = (N-1) H`, parity corollary | `HamPath`, `ArcParity`, `Redei` | FORMALIZED (`mem_hamPathsR`, `nodup_hamPathsR`, `sum_arcCount`, `allOdd_mod_four`), and **Rédei's theorem is formalized too** (`redei`), so the corollary is unconditional (`allOdd_mod_four_of_tournament`) |
| 3 Collatz drop multiplicity (S15 Thm 2) and label map (Thm 1) | `CollatzDrop` | FORMALIZED |
| 4 Shaved: `H_4` in every 4-tournament; exactly the `C3[1,C3,1]` copies avoid `H_5`; `#copies(H_n) = H - n hc` | `Shaved`, `CycleCount` | FORMALIZED (`every_four_contains_H4`, `avoids_H5_iff`, `copiesH_add`) |
| extra: THM-4525 explicit upper bounds | `KeyBlocks`, `EdimQ7*`, `EdimQ8*`, `EdimQ9*` | FORMALIZED: `edim_m(Q_7) <= 19`, `edim_m(Q_8) <= 26`, `edim_m(Q_9) <= 38` (`Q_10..Q_12` not: kernel cost) |
| extra: THM-4524 C1 (Cayley all-even) for cyclic groups | `Circulant` | FORMALIZED (`circulant_arcCount_even`: every circulant digraph on `Z_n`, `n` odd) |
| extra: THM-4524 A2 (loop-set sum) | `SwitchSum` | FORMALIZED (`two_loop_sets`, `switch_sum`: `Σ_L H(switch_L T) = 2 N!`; `hpCount_complete`: `N!` orderings) |
| extra: THM-4524 §3, no HP-covered tournament has an arc on every HP; HP existence | `HPExist` | FORMALIZED (`no_covered_universal_arc`, `exists_hamPath`) |
| extra: **Rédei's theorem** and corollaries (THM-4524 (C) unconditional; Rédei shaving; THM-4526 Theorem A for even `n`) | `Redei` | FORMALIZED (`redei`, `allOdd_mod_four_of_tournament`, `allEven_odd_of_tournament`, `oddArcs_parity_of_tournament`, `exists_odd_arc_of_even`, `shave_keeps_odd_iff`, `copiesH_odd_of_even`, `containsH_of_even`) |
| extra: THM-4524 A3(e), constant-`H` switching classes need `N = 2^k` | `ConstantH` | FORMALIZED (`constant_switching_class`, via `switch_sum`, `redei` and Legendre's formula `legendre_two`; attained at `N = 4`, `constant_class_example`) |
| extra: THM-4525 L3 (alternating sum) | `AltSum` | FORMALIZED (`edgeDist_projection`, `alt_sum`, `alt_sum_hist`, `alt_sum_two_values`) |
| extra: THM-4524 §3, parity-break number `ρ = 2` for all-even tournaments | `ParityBreak` | FORMALIZED (`parity_break_two`, with `delete_two`, `out_arcs_sum`, `end_sum`, `arc_split`) |
| out of scope / not done | | `edim_m(Q_6) >= 15` (1.4e10 leaves); L5; THM-4524 A3 (a)-(d), B1, C2, E1, the dead-arc formula, the census, the Hall bound; THM-4526 Theorem A for odd `n > 5`, Theorems B, C; THM-4527 Props 3-5 |

## 1. Reproduction and resources

```text
cd 04-computation/lean/ProcgenSelfieEdim && python3 verify.py
```

`verify.py` scans for proof escapes, checks that every `theorem` of the 45 library modules is listed in
`AxiomAudit.lean`, runs `lake clean`, builds the modules **one at a time** in dependency order (at most one Lean
process), runs `lake env lean AxiomAudit.lean` (which imports only the library root), rejects any axiom outside
`{propext, Quot.sound}`, and writes `verification.json` (axioms per theorem, sha256 of all sources, build times,
peak child RSS) and `verification.log` (the transcript). The lakefile passes `--threads=1 --memory=1400` to Lean,
so a kernel check that would exceed about 1.4 GB aborts with "excessive memory consumption" instead of swapping.

Measured (final run): clean sequential build **129 s** in total over 45 modules (largest: `SelfieFinite` 13.3 s,
`Shaved` 6.5 s, the twelve `EdimQ9P*` parts 4-6.4 s each), axiom audit 0.4 s, **peak RSS 1028 MB** (one Lean
process at a time; about 400 MB of it is the baseline for loading `Init`/`Std`). Every module stays below the
1.4 GB cap set in the lakefile. The build log has no warnings.

## 2. Toolchain findings: which tactics introduce `Classical.choice` (Lean 4.30.0)

The previous instance had started this investigation; it is now complete. Each row was checked with
`#print axioms` on small test theorems (scratch), and the rules were then enforced on the whole package.

| Tactic / lemma | Axioms |
|---|---|
| `decide`, `decide +kernel`, `rfl`, `exact` with term proofs | none |
| `by_cases h : p` with a `Decidable p` instance (e.g. `Nat` equalities, `Bool`) | none |
| `by_cases h : p` without an instance (classical fallback) | `propext, Classical.choice, Quot.sound` |
| `omega`, goal atomic, `False`, or a disjunction of atoms | `propext, Quot.sound` |
| `omega`, goal containing a conjunction (`A ∧ B`, `(A ∧ B) ∨ C`) | **`Classical.choice`** |
| `omega` on a non-arithmetic goal (e.g. an existential), closing it from contradictory hypotheses | **`Classical.choice`** (it falls back to classical `by_contra`; `exfalso` first avoids it) |
| `omega` with a hypothesis `(2^1 - 3 : Int) ∣ x` in context | produced a kernel "application type mismatch" (an `omega` bug with a power literal inside `∣`); clearing the hypothesis avoids it |
| `simp` (typical rewriting) | `propext` |
| `simp [h]` proving `(a == b) = true` from `h : a = b`; `simp` closing `0 = k + 1` or `(x :: l).length = 0` | **`Classical.choice`** |
| `beq_iff_eq.2 h` for the same goal | `propext, Quot.sound` |
| `beq_self_eq_true` | **`Classical.choice`** |
| `grind` | **`Classical.choice`** |
| `split` | none |
| `funext` | `Quot.sound` |
| core `List.nodup_range`, `List.length_erase_of_mem`, `List.mem_erase_of_ne` | **`Classical.choice`** |
| core `List.mem_flatMap`, `List.nodup_append`, `List.mem_range`, `List.filter_congr`, `List.append_inj`, `Int.mul_ediv_cancel'` | `propext` and/or `Quot.sound` |

Adopted rules: no `grind`; `omega` only on goals without `∧` (split conjunctions first) and after `exfalso` when the
goal is not arithmetic; `beq_iff_eq` instead of `simp` for `==` goals; own `Classical`-free proofs of
`nodup_range`, pigeonhole (core has no `List.Subperm`), and filter-length lemmas. The audit tooling
(`verify.py`) treats any `Classical.choice` as a failure.

## 3. THM-4525: the edge multiset dimension of `Q_6`

### Definitions (module `Hypercube`, verbatim)

```lean
def hamming : Nat → Nat → Nat → Nat
  | 0, _, _ => 0
  | d + 1, u, v => (if u % 2 = v % 2 then 0 else 1) + hamming d (u / 2) (v / 2)
def IsEdge (d u v : Nat) : Prop := u < 2 ^ d ∧ v < 2 ^ d ∧ hamming d u v = 1
def edgeDist (d u v s : Nat) : Nat := min (hamming d u s) (hamming d v s)
def hist (d : Nat) (S : List Nat) (u v r : Nat) : Nat :=
  (S.filter (fun s => edgeDist d u v s == r)).length
def SameEdge (u v u' v' : Nat) : Prop := (u = u' ∧ v = v') ∨ (u = v' ∧ v = u')
def Resolving (d : Nat) (S : List Nat) : Prop :=
  ∀ u v u' v', IsEdge d u v → IsEdge d u' v' →
    (∀ r, hist d S u v r = hist d S u' v' r) → SameEdge u v u' v'
def antipode (d u : Nat) : Nat := 2 ^ d - 1 - u
```

Vertices of `Q_d` are `0, …, 2^d - 1` (bit `i` = coordinate `i`); `hamming d` is the graph distance.
`Resolving` quantifies over all edges and all levels `r` (equal histograms = equal distance multisets).

### Theorems

* **(a)** `paperSet_resolving : Resolving 6 paperSet` with
  `paperSet = [1, 3, 13, 15, 17, 22, 29, 31, 33, 37, 44, 45, 51, 53, 57]`;
  `paperSet_is_vertex_set` checks length 15, no repetition, all `< 64`, and equality with the vertices of the
  paper's mask `0x02283022a042a00a`. Proof: kernel check `check_paperSet : checkResolving 6 paperSet = true`
  (no axioms) plus the soundness theorem `resolving_of_check`. The checker compares one numeric key per edge,
  `Σ_{s ∈ S} 16^{d(e,s)}`; soundness needs only that the key is a function of the histogram
  (`weighted_sum_eq : Σ_{s∈S} g(d(e,s)) = Σ_{r<d} g(r) H_e(r)`), not injectivity. `edgeList_six`: `Q_6` has 192
  edges.
* **(b) L2** `hist_antipode (d) (S) (u v) (he : IsEdge d u v) (r) (hr : r < d) :
  hist d S (antipode d u) (antipode d v) r = hist d S u v (d - 1 - r)`, for every `d`. Supporting:
  `hamming_antipode`, `edgeDist_antipode` (`d(ē, s) = d - 1 - d(e, s)`), `isEdge_antipode`.
* **(c) L4** `not_resolving_of_length_le_six (S : List Nat) (hS : ∀ s ∈ S, s < 2 ^ 6) (hm : S.length ≤ 6) :
  ¬ Resolving 6 S`, i.e. `edim_m(Q_6) >= 7`. It holds even for lists with repetitions. The proof follows the
  lane note: an edge with `H(0) >= 1` or `H(5) >= 1` contains a landmark or a landmark's antipode
  (`touched_of_bad`); each vertex and its antipode touch 12 edges (`touch_count`, a 64-vertex kernel check);
  union bound `length_filter_any_le`; the good edges have histograms `[0, a, b, c, d, 0]`, `a + b + c + d = m`,
  in the explicit list `comps m` (`histVec6_mem_comps`, using `rsum_hist : Σ_r H_e(r) = |S|`); resolvability
  makes them distinct; `pigeonhole`; `comps_small : |comps m| + 12 m < 192` for `m <= 6`.
* **(d) L1** `trivial_stabilizer (d) (hd : 2 ≤ d) (S) (hS : ∀ s ∈ S, s < 2 ^ d) (hnd : S.Nodup)
  (hres : Resolving d S) (σ) (hσ : IsAutomorphism d σ) (hfix : ∀ s ∈ S, σ s ∈ S) : ∀ u, u < 2 ^ d → σ u = u`.
  `IsAutomorphism d σ` = a bijection of `{0, …, 2^d-1}` (with an explicit inverse) preserving adjacency both ways.
  `automorphism_isometry` proves such maps preserve Hamming distance (`geodesic_step` + triangle inequality);
  `hist_invariant` gives `H_{σe} = H_e`; then every edge is fixed as a pair, and a vertex on the two distinct
  edges `flipBit u 0`, `flipBit u 1` is fixed. Corollary `paperSet_trivial_stabilizer`.

Not formalized: `edim_m(Q_6) >= 15` (out of scope), L5, the 229-orbit census, the `Q_10..Q_12` sets
(the `Q_7`, `Q_8`, `Q_9` sets are, section 6b; L3 is in section 6e).

## 4. THM-4524: selfie tournaments and arc parity

### Definitions (module `Tournament`, `HamPath`)

```lean
def IsTournament (N : Nat) (T : Nat → Nat → Bool) : Prop :=
  ∀ a b, a < N → b < N → a ≠ b → T b a = !T a b
def SameOn (N : Nat) (T T' : Nat → Nat → Bool) : Prop := ∀ a b, a < N → b < N → a ≠ b → T a b = T' a b
def switch (L : Nat → Bool) (T : Nat → Nat → Bool) (a b : Nat) : Bool := if L a = L b then T a b else T b a
def tiling (t : Nat → Nat → Bool) (a b : Nat) : Bool :=
  if a = b + 1 then true else if b = a + 1 then false else if a < b then t a b else !t b a
def SameTiles (N : Nat) (t t' : Nat → Nat → Bool) : Prop := ∀ a b, a + 2 ≤ b → b < N → t a b = t' a b
def SameLoops (N : Nat) (L L' : Nat → Bool) : Prop := ∀ a, a < N → L a = L' a
def compl (L : Nat → Bool) (a : Nat) : Bool := !L a
def IsPath (T : Nat → Nat → Bool) : List Nat → Prop
  | a :: b :: rest => T a b = true ∧ IsPath T (b :: rest)
  | _ => True
def IsHamPath (n : Nat) (T : Nat → Nat → Bool) (p : List Nat) : Prop :=
  p.length = n ∧ p.Nodup ∧ (∀ x ∈ p, x < n) ∧ IsPath T p
def usesArc (u v : Nat) : List Nat → Bool
  | a :: b :: rest => (Nat.beq a u && Nat.beq b v) || usesArc u v (b :: rest)
  | _ => false
noncomputable def hpCount (n : Nat) (T : Nat → Nat → Bool) : Nat := (hamPathsR n T).length
noncomputable def arcCount (n : Nat) (T : Nat → Nat → Bool) (u v : Nat) : Nat :=
  ((hamPathsR n T).filter (usesArc u v)).length
```

Change of convention: the base path is `N-1 → ⋯ → 0` (the note's `N → ⋯ → 1` shifted by one); tile bit
`t a b` (`a + 2 <= b`) says `a → b`. `usesArc_iff` proves `usesArc u v p = true ↔ ∃ l1 l2, p = l1 ++ u :: v :: l2`.

### (a) Gauge theorem A1 (all `N`)

* `path_arc : switch L (tiling t) (i + 1) i = true ↔ L i = L (i + 1)` (the base-path arc `i+1 → i` is reversed
  iff exactly one of `i, i+1` is looped).
* `selfie_fibre : SameOn N (switch L (tiling t)) (switch L' (tiling t')) ↔ SameTiles N t t' ∧ (SameLoops N L L' ∨ SameLoops N L (compl L'))`.
* `loops_ne_compl (hN : 1 ≤ N) : ¬ SameLoops N L (compl L)` (the two fibre points are distinct).
* `selfie_onto (hT : IsTournament N T) : ∃ t L, SameOn N (switch L (tiling t)) T` (via the explicit gauge
  `gauge`, `gauge_path`).

Together: `(t, L) ↦ switch_L(T_t)` is exactly 2-to-1 onto labeled tournaments with fibres `{(t, L), (t, Lᶜ)}`.
Also `reversedCount_parity` (#reversed path arcs `≡ [0 ∈ L] + [N-1 ∈ L] mod 2`), `selfie_parameter_count`
(`pairs (N-1) + N = pairs N + 1`, `pairs n = C(n,2)` by `pairs_eq : 2 * pairs n = n * (n - 1)`).

### (c) Infrastructure: verified enumerator, double counting, parity

* `mem_hamPathsR : q ∈ hamPathsR n T ↔ IsHamPath n T q` and `nodup_hamPathsR : (hamPathsR n T).Nodup`, for every
  relation `T` and every `n`. The enumerator is a DFS that grows paths backwards along arcs of `T`.
  `hpCount_unique` / `arcCount_unique`: any duplicate-free list of exactly the HPs (through `u → v`) has length
  `hpCount` (`arcCount`).
* `sum_arcCount : rsum n (fun u => rsum n (fun v => arcCount n T u v)) = (n - 1) * hpCount n T` (any relation).
* `arcTotal (hT : IsTournament n T) : rsum n (fun u => rsum n (fun v => ind (isArc T u v))) = pairs n`.
* With **`hH : hpCount n T % 2 = 1` as a hypothesis** (discharged by `redei`, section 6c, which gives the unconditional
  versions `oddArcs_parity_of_tournament`, `allOdd_mod_four_of_tournament`, `allEven_odd_of_tournament`):
  * `oddArcs_parity : rsum n (fun u => rsum n (fun v => arcCount n T u v % 2)) % 2 = (n - 1) % 2`;
  * `allOdd_mod_four (hn : 1 ≤ n) (hT : IsTournament n T) (hodd : ∀ u v, u < n → v < n → u ≠ v → T u v = true → arcCount n T u v % 2 = 1) (hH) : n % 4 = 1 ∨ n % 4 = 2`;
  * `allEven_odd (hn : 1 ≤ n) (heven : ∀ u v, u < n → v < n → arcCount n T u v % 2 = 0) (hH) : n % 2 = 1`.
* `shave_arc : hpCount n (deleteArc T u v) + arcCount n T u v = hpCount n T` (THM-4524 §3, single-arc shaving).

### (b) Finite facts (kernel checks, module `SelfieFinite`)

* `qr7_counts : hpCount 7 qr7 = 189 ∧ ∀ u v, u < 7 → v < 7 → qr7 u v = true → arcCount 7 qr7 u v = 54`,
  `qr7 a b` iff `b - a ∈ {1, 2, 4} (mod 7)`.
* `qr7del_counts (z) (hz : z < 7) : hpCount 6 (qr7del z) = 45 ∧ ∀ u v, u < 6 → v < 6 → qr7del z u v = true →
  (antipodalTo z u v = true ∧ arcCount 6 (qr7del z) u v = 23) ∨ (antipodalTo z u v = false ∧ arcCount 6 (qr7del z) u v = 13)`,
  for **every** deleted vertex `z`; `antipodalTo z u v` iff `u + v ≡ 2z (mod 7)` in `QR_7` labels;
  `qr7del_antipodal_arcs`: exactly 3 such arcs (so 12 arcs with 13). `qr7del_allOdd`: all 15 arcs odd.
* `no_allOdd_three_to_five (hn : 3 ≤ n ∧ n ≤ 5) (hT : IsTournament n T) : ∃ u v, u < n ∧ v < n ∧ u ≠ v ∧ T u v = true ∧ arcCount n T u v % 2 = 0`
  (exhaustive over all `2^C(n,2)` labeled tournaments: 8, 64, 1024; the 1024 in 8 chunks of 128). The reduction
  `tourC_code` shows every tournament on `n` vertices agrees with the coded tournament `tourC c` for some
  `c < 2 ^ pairs n`; `hpCount_congr` / `arcCount_congr` transfer counts along `SameOn`.

Not formalized: A3 (a)-(d) (the loop-Walsh expansion and degree bound), B1, C1 for non-cyclic groups, C2, the
dead-arc formula, the `N <= 10` census, the HP-blocking number and the Hall bound (HYP-9168), E1 (OPEN-Q-060).
Formalized beyond this section: A2 and C1 for cyclic groups (section 6b), Rédei's theorem (6c), A3(e) (6d), and the
parity-break number `ρ = 2` (6e).

## 5. THM-4526 (S15): shaved tournaments

* `every_four_contains_H4 (hT : IsTournament 4 T) : ∃ A B C D, A < 4 ∧ B < 4 ∧ C < 4 ∧ D < 4 ∧ [A, B, C, D].Nodup ∧ T A B = true ∧ T B C = true ∧ T C D = true ∧ T A D = true`.
* `avoids_H5_iff (hT : IsTournament 5 T) : ¬ ContainsH 5 T ↔ IsoTo 5 T c3c3`, where
  `ContainsH n T := ∃ a b mid, IsHamPath n T (a :: (mid ++ [b])) ∧ T a b = true` (path + first→last arc),
  `IsoTo n T S := ∃ π, (∀ i < n, π i < n) ∧ (π injective on {0..n-1}) ∧ ∀ i j, i < n → j < n → i ≠ j → S (π i) (π j) = T i j`,
  and `c3c3` is the labeled `C3[1,C3,1]` (`0` beats the 3-cycle `1 → 2 → 3 → 1`, which beats `4`, and `4 → 0`).
  Direction ⇒ is a kernel check over the 1024 codes (each code has a copy of `H_5` or a permutation mapping it
  onto `c3c3`; 8 chunks); direction ⇐ is a proof (isomorphisms push copies forward; `c3c3_avoids_H5`).
* `copiesH_add (hn : 2 ≤ n) (hT : IsTournament n T) : copiesH n T + n * hcCount n T = hpCount n T`, i.e.
  `#copies(H_n) = H - n · hc`. `copiesH` counts HPs whose first vertex beats the last
  (`containsH_iff_copiesH_pos`); `hcCount` counts closing HPs starting at vertex `0` (one per directed Hamiltonian
  cycle). `closing_count` (any relation, `n >= 2`): closing HPs `= n · hcCount`, by an explicit rotation map
  (`rotN`), position uniqueness of `0`, and a duplicate-free enumeration of all rotations. Example `c3c3_counts`:
  `H = 15`, `hc = 3`, no copy of `H_5`.

Theorem A for every even `n ≥ 2` (an odd number of copies, so `H_n` is unavoidable) is `copiesH_odd_of_even` /
`containsH_of_even` (section 6c). Not formalized: `u(n)`, `κ(7) = 12`, `κ(8) = 17`, Theorem A for odd `n > 5`,
Theorem C.

## 6. S15 Theorems 1-2 (THM-4527): the Syracuse drop and the label map

Definitions (`CollatzDrop`): `v2 n`, `oddPart n` by fuel-bounded halving, with `v2_oddPart_of_eq`
(uniqueness of `n = 2^v s`, `s` odd) and `decomp` (existence); `syr A := oddPart (3 * A + 1)`;
`drop A : Int := ((A : Int) - (syr A : Int)) / 2`; `labelF M := (syr (2 * M - 1) + 1) / 2`;
`Admissible m v := 1 ≤ v ∧ ((2 : Int) ^ v - 3) ∣ (6 * m + 1) ∧ 0 < (6 * m + 1) / ((2 : Int) ^ v - 3)`.

* `drop_identity (hA : A % 2 = 1) : 6 * drop A + 1 = ((2 : Int) ^ v2 (3 * A + 1) - 3) * (syr A : Int)`.
* `admissible_iff : Admissible m v ↔ ∃ A : Nat, A % 2 = 1 ∧ 0 < A ∧ drop A = m ∧ v2 (3 * A + 1) = v`, with
  `drop_injective` (an odd `A` is determined by `drop A` and `v2 (3A+1)`): so for each `m ∈ Z`, `A ↦ v_2(3A+1)` is a
  bijection from `{A odd > 0 : K(A) = m}` onto `{v : Admissible m v}` (Theorem 2(b)); `drop_forward`,
  `drop_backward` are the two halves (the preimage is `A = 2m + s`, `s = (6m+1)/(2^v-3)`).
* `admissible_of_neg (hm : m < 0) : Admissible m v ↔ v = 1`;
  `admissible_of_nonneg (hm : 0 ≤ m) : Admissible m v ↔ 2 ≤ v ∧ ((2 : Int) ^ v - 3) ∣ (6 * m + 1)`;
  `mod_four_iff (hA : A % 2 = 1) : A % 4 = 1 ↔ 2 ≤ v2 (3 * A + 1)` (descents are exactly `v >= 2`).
* Label map (Theorem 1): `labelF_even (hN : 1 ≤ N) : labelF (2 * N) = 3 * N`,
  `labelF_four_j_one : labelF (4 * j + 1) = 3 * j + 1`, `labelF_microcosm (hn : 1 ≤ n) : labelF (4 * n - 1) = labelF n`,
  `labelF_eq_sub_drop (hM : 1 ≤ M) : (labelF M : Int) = (M : Int) - drop (2 * M - 1)` (Theorem 1(a)),
  `drop_two_regular` (`K_{2j+1} = j`, `K_{4m} = 5m - 1`, `K_{4m-2} = K_m + 6m - 4`, Theorem 1(c)).
* Worked values: `owner_K_values` (`K_N`, `N = 1..11`: `0, 2, 1, 4, 2, 10, 3, 9, 4, 15, 5`), `owner_S_values`,
  `multiplicity_examples` (`m = 1` once, `m = 2` twice at `v = 2, 4`, `m = 24` three times at `v = 2, 3, 5`).

Not formalized: Proposition 3 (Paley restatement), Propositions 4-5, the mean multiplicity `1.3437…`.

## 6b. Extras beyond the orchestrator's list

* **THM-4525 upper bounds (`KeyBlocks`, `EdimQ7*`, `EdimQ8*`, `EdimQ9*`).**
  `q7set_resolving : Resolving 7 q7set` (19 vertices), `q8set_resolving : Resolving 8 q8set` (26),
  `q9set_resolving : Resolving 9 q9set` (38): the edim lane's annealed sets (`procgen_edim_20261001_run.py`, S10).
  `q7set_is_vertex_set` etc. check size, no repetition, and range. Method: the one-pass keys
  `Σ_{s∈S} base^{d(e,s)}` (`base` 32 or 64; `resolving_of_edgeKeys` needs only that the key is a function of the
  histogram) are compared with precomputed numbers in blocks of 64 edges (`eq_of_blocks`, one declaration per
  block, spread over several modules), and the numbers are shown distinct by a verified merge sort
  (`msortR_perm` proves only "permutation of the input"; strict increase is checked by computation). The sets
  for `Q_10..Q_12` are not formalized: the key blocks read `edgeList d` from its start, and the late blocks of
  `Q_10` already cost about 1.2 GB.
* **THM-4524 C1 for `Z_n` (`Circulant`).** `circulant_arcCount_even (hn : n % 2 = 1) (D) (hu : u < n) (hv : v < n) :
  arcCount n (circ n D) u v % 2 = 0`, with `circ n D a b = D ((b - a) mod n)`. The proof is the note's
  involution `P ↦ reverse (φ P)`, `φ(x) = u + v - x`, and a lemma that a fixed-point-free involution on a finite
  set has even size (`even_of_involution`). It holds for every circulant digraph, not only tournaments; only
  cyclic groups are covered (not `Z_3 × Z_3` etc.). `qr7_allEven` derives the `QR_7` case without enumeration.
* **THM-4524 A2 (`SwitchSum`).** `two_loop_sets`: for a tournament `T` and an ordering `π` there is `L₀` with
  `IsPath (switch L T) π ↔ SameLoops N L L₀ ∨ SameLoops N L (compl L₀)`; `count_loop_sets`: exactly 2 of the
  `2^N` bit lists; `switch_sum (hN : 1 ≤ N) (hT : IsTournament N T) :
  ((allBools N).map (fun bs => hpCount N (switch (nthB bs) T))).sum = 2 * fact N`, with `hpCount_complete :
  hpCount N complete = fact N` and `length_allBools : (allBools N).length = 2 ^ N`.
* **THM-4524 §3 (`HPExist`).** `exists_hamPath` (every tournament has a Hamiltonian path, by insertion) and
  `no_covered_universal_arc (hn : 3 ≤ n) (hT) (hcover : every arc u → v, u ≠ v, lies on some HP) :
  ¬ ∃ u v, u < n ∧ v < n ∧ u ≠ v ∧ T u v = true ∧ ∀ p, IsHamPath n T p → usesArc u v p = true`.
* **L2 corollaries (`Hypercube`).** `hist_not_palindrome` (`d` even, `|S|` odd: no histogram is a palindrome) and
  `no_antipodal_collision` (then no edge collides with its antipodal edge).

## 6c. Rédei's theorem (`Redei`)

`redei : ∀ (n : Nat) (T : Nat → Nat → Bool), IsTournament n T → hpCount n T % 2 = 1`: every tournament has an odd
number of Hamiltonian paths (Rédei 1934). Axioms: `propext, Quot.sound`.

The proof is by induction on `n` and deletes the top vertex `x = n`, so no relabeling is needed (`T - x` is `T` on
`{0, …, n-1}`):

1. `hp_split`, `hp_decomp`: each Hamiltonian path of `T` has `x` first, last, or between exactly one pair `a, b`.
2. `count_head`, `count_last`: the first two classes are the HPs `Q` of `T - x` with `x → first Q`, resp. `last Q → x`.
3. `count_interior`: the paths through `a → x → b` correspond (removing `x`, `eraseN`) to the HPs through `a → b` of
   `forceArc T a b` (`T - x` with the pair `{a, b}` oriented `a → b`), times `[a → x ∧ x → b]`.
4. If `b → a` in `T`, single-arc shaving (`shave_arc`) for `forceArc T a b` and for `T` gives
   `c_F(a → b) + H(T') = c'(b → a) + H(F)`; both `H` are odd by induction, so `c_F(a → b) ≡ c'(b → a)`.
   If `a → b`, then `F = T'` on `{0, …, n-1}`.
5. So `H(T) ≡ Σ_Q ([x → first Q] + [last Q → x] + #{consecutive pairs of Q separated by x})` (`weighted_arcs`,
   `rsum_comm`), and each summand is odd (`odd_contribution`: the side of `x` changes an odd number of times
   along `0, s_1, …, s_n, 1`). Hence `H(T) ≡ H(T - x) ≡ 1`.

Corollaries (all `propext, Quot.sound`): `allOdd_mod_four_of_tournament` (all arcs odd ⇒ `n ≡ 1, 2 mod 4`),
`allEven_odd_of_tournament` (all arcs even ⇒ `n` odd), `oddArcs_parity_of_tournament` (#odd arcs `≡ n - 1`),
`exists_odd_arc_of_even` (for even `n ≥ 2` some arc is odd), `shave_keeps_odd_iff` (THM-4524 §3: deleting `u → v`
keeps `H` odd iff `c(u → v)` is even), `copiesH_odd_of_even` and `containsH_of_even` (THM-4526 Theorem A, even
case: a tournament on an even number `n ≥ 2` of vertices contains an odd number of copies of `H_n`).

## 6d. Constant-`H` switching classes (`ConstantH`, THM-4524 A3(e))

`constant_switching_class (hN : 1 ≤ N) (hT : IsTournament N T) (hconst : ∀ L : Nat → Bool, hpCount N (switch L T) = hpCount N T) :
2 ^ (N - 1) * hpCount N T = fact N ∧ ∃ k, N = 2 ^ k`. Proof: `switch_sum` gives `2^N · H = 2 · N!`; `H` is odd
(`redei`), so `v_2(N!) = N - 1` (`v2_oddPart_of_eq`); Legendre's formula `legendre_two : v2 (fact n) + s2 n = n`
(proved from `s2_succ : s_2(n+1) + v_2(n+1) = s_2(n) + 1` and `v2_mul`) gives `s_2(N) = 1`, so `N` is a power of two
(`pow_two_of_s2_eq_one`). The case `N = 4` occurs: `constant_class_example : ∀ L, hpCount 4 (switch L (tourC 16)) = 3`
(a kernel check over the 16 loop sets).

## 6e. THM-4525 L3 (`AltSum`) and THM-4524's parity-break number (`ParityBreak`)

* **L3.** For the edge `e = {u, u + 2^i}` (`bit u i = 0`) and `χ(w) = (-1)^{Σ_{j ≠ i} w_j}` (`chi d i w`, with
  `sgn n = (-1)^n`): `edgeDist_projection : edgeDist d u (u + 2 ^ i) s = rsum d (fun j => if j = i then 0 else
  (if bit u j = bit s j then 0 else 1))` (the projection lemma); `alt_sum : (S.map (fun s => sgn (edgeDist d u (u + 2 ^ i) s))).sum
  = chi d i u * (S.map (chi d i)).sum`; `alt_sum_hist` (the same with `Σ_{r even} H_e(r) - Σ_{r odd} H_e(r)` on the
  left); `alt_sum_two_values` (within one direction only `± Σ_s χ(s)` occur).
* **`ρ = 2` for all-even tournaments.** `parity_break_two (hn : 2 ≤ n) (hT : IsTournament n T) (heven : all arc
  counts even) : (∀ u v, u < n → v < n → hpCount n (deleteArc T u v) % 2 = 1) ∧ ∃ u v w, … T u v = true ∧ T v w = true ∧
  hpCount n (deleteArc (deleteArc T u v) v w) % 2 = 0`. Ingredients: `out_arcs_sum` (`Σ_w c(v → w) + end(v) = H`, so
  `end(v)` is odd), `end_sum` (`end(v) = Σ_u #(paths ending u → v)`), `arc_split` (`c(u → v) = #(ending u → v) +
  Σ_w c(u → v → w)`), and the inclusion-exclusion `delete_two` (`H(T - e - f) + c(e) + c(f) = H(T) + c(u → v → w)`).

## 7. Kernel cost: lessons (supersede the first instance's notes)

* Structural recursion compiles to `brecOn`, which allocates a tuple per step of kernel reduction, and
  `BEq`/`DecidableEq` instance chains are expensive. On the `Q_6` check the first (list-of-histograms) checker
  peaked at **1.6 GB**; with recursor-based functions (`Nat.rec`/`List.rec`, `Nat.beq`, `Nat.ble`, `cond`) and
  a one-pass numeric key the module needs about 2 s and stays near 670 MB. Each fast function is proved equal to
  its readable version.
* The verified DFS enumerator is cheap: `QR_7` (189 paths, all 49 ordered pairs) and all seven `QR_7 - z` checks
  run together in about 2 s. A chunk of 128 five-vertex tournaments costs about 1.7 s and 660 MB; 256 would cost
  950 MB, hence chunks of 128 (memory grows about 1.1 MB per tournament inside one declaration).
* `--memory=N` (`-M`) tracks RSS closely (with `-M 550` a 605 MB check aborted) and makes a runaway check fail
  cleanly; it is set to 1400 in the lakefile. Lake has no job limit, so `verify.py` builds module by module.
* Kernel caches are per declaration, but the RSS of one Lean process keeps growing across declarations (about
  50 MB per heavy declaration). The `Q_8` certificate in a single file hit the 1.4 GB cap after 13 blocks; split
  into modules of 3-4 blocks, every part stays below 900 MB.
* Comparing the keys of one block with literals is cheap; the cost of a late block is dominated by walking the
  spine of `edgeList d` (core `flatMap`/`filterMap`) up to that block, which is why `Q_10` (5120 edges) was not
  done. Distinctness is cheap with the verified merge sort: 1024 keys in under 1 s, 2304 keys in about 2 s.
* A literal list with thousands of entries overflows the elaborator's recursion depth; the keys are stored as
  64-element block literals (`noncomputable def … : List (List Nat)`) and flattened.

## 8. Files

* `04-computation/lean/ProcgenSelfieEdim/`:
  * `lakefile.toml`, `lean-toolchain`, `lake-manifest.json`, `.gitignore` (excludes `.lake/` and build outputs);
  * `ProcgenSelfieEdim.lean` (the root) and `ProcgenSelfieEdim/*.lean`: `ListLemmas`, `Hypercube`, `EdimQ6`,
    `CountingBound`, `Stabilizer`, `KeyBlocks`, `EdimQ7Data`, `EdimQ7P0-1`, `EdimQ7`, `EdimQ8Data`, `EdimQ8P0-3`,
    `EdimQ8`, `EdimQ9Data`, `EdimQ9P0-11`, `EdimQ9`, `Tournament`, `HamPath`, `ArcParity`, `TournamentCode`,
    `SelfieFinite`, `AltSum`, `Shaved`, `CycleCount`, `Circulant`, `HPExist`, `SwitchSum`, `Redei`, `ConstantH`,
    `ParityBreak`, `CollatzDrop`;
  * `gen_certificates.py` (regenerates the `EdimQ7*`, `EdimQ8*`, `EdimQ9*` modules from the lane's sets; the
    Lean files are trusted only through the kernel checks);
  * `AxiomAudit.lean`, `verify.py`, `README.md`, `verification.json`, `verification.log`.
* This note.
* Scratch (not for commit): `scratch/procgen_lean/` (axiom tests, build helpers, cost experiments).

No git state was changed. No HYP/THM file was created or edited.
