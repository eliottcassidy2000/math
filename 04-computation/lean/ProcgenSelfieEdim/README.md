# ProcgenSelfieEdim: Lean formalization of THM-4524, THM-4525 and S15's Collatz / shaving results

**Status: PROVED in Lean 4.30.0 (core + bundled `Std`, no Mathlib), within the scopes below.**
Every audited theorem depends only on `propext` and/or `Quot.sound` (many on no axiom at
all). There is no `sorry`, no `native_decide` / `decide +native`, and no custom axiom.
Finite facts are checked by kernel reduction (`decide +kernel`).

Lane: procgen lean lane, session `collatz-procgen-20260922`, 2026-10-01. The companion note
is `05-knowledge/results/procgen_lean_20261001_formalization.md`.

## What is proved

| Module | Content | Main theorems |
|---|---|---|
| `Hypercube` | `Q_d` on `{0, …, 2^d-1}`, Hamming distance, edges, `d(e,s) = min(d(u,s), d(v,s))`, histograms `H_e(r)`, edge-multiset resolving; a sound checker | `hist_antipode` (**L2**, all `d`), `hist_not_palindrome`, `no_antipodal_collision` (L2 corollaries: `d` even, `|S|` odd), `resolving_of_check` |
| `EdimQ6` | the paper's set `{1,3,13,15,17,22,29,31,33,37,44,45,51,53,57}` (mask `0x02283022a042a00a`) | `paperSet_resolving` (**THM-4525 (a)**: `edim_m(Q_6) <= 15`) |
| `CountingBound` | counting + pigeonhole | `not_resolving_of_length_le_six` (**L4**: no list of `<= 6` vertices resolves `Q_6`, so `edim_m(Q_6) >= 7`) |
| `Stabilizer` | graph automorphisms of `Q_d` are Hamming isometries; histograms are stabilizer-invariant | `trivial_stabilizer` (**L1**, `d >= 2`), `paperSet_trivial_stabilizer` |
| `AltSum` | THM-4525 **L3** | `edgeDist_projection` (`d({u, u+2^i}, s) = #{j ≠ i : u_j ≠ s_j}`), `alt_sum` (`Σ_s (-1)^{d(e,s)} = χ(u) Σ_s χ(s)`), `alt_sum_hist` (histogram form), `alt_sum_two_values` |
| `KeyBlocks`, `EdimQ7*`, `EdimQ8*`, `EdimQ9*` | certificates for the lane's explicit sets: keys checked in blocks of 64 edges (one declaration and module part each), distinctness by a verified merge sort | `q7set_resolving` (`edim_m(Q_7) <= 19`), `q8set_resolving` (`<= 26`), `q9set_resolving` (`<= 38`) |
| `Tournament` | tournaments on `{0, …, N-1}`, switching, the tiling model (base path `N-1 → ⋯ → 0`) | `path_arc`, `selfie_fibre`, `loops_ne_compl`, `selfie_onto` (**A1**: exactly 2-to-1, fibres `{L, Lᶜ}`), `reversedCount_parity`, `selfie_parameter_count` |
| `HamPath` | Hamiltonian paths; a DFS enumerator proved to list every HP exactly once | `mem_hamPathsR`, `nodup_hamPathsR`, `hpCount_unique`, `arcCount_unique` |
| `ArcParity` | double counting and parity | `sum_arcCount` (`Σ c(e) = (N-1) H`), `oddArcs_parity`, `allOdd_mod_four`, `allEven_odd` (with `H` odd as a hypothesis; unconditional versions in `Redei`), `arcTotal`, `shave_arc` (`H(T-e) = H(T) - c(e)`) |
| `TournamentCode` | every labeled tournament is `tourC c` for some code `c < 2^C(n,2)` | `tourC_code` |
| `SelfieFinite` | finite facts (kernel checks) | `qr7_counts` (`H(QR_7) = 189`, every arc on 54 HPs), `qr7del_counts` / `qr7del_allOdd` (for every `z`, `QR_7 - z` has `H = 45`, 23 HPs on the 3 `z`-antipodal arcs, 13 on the others), `no_allOdd_three_to_five` |
| `Circulant` | THM-4524 **C1** for the cyclic groups | `circulant_arcCount_even` (every circulant digraph on `Z_n`, `n` odd: every arc on an even number of HPs), `qr7_allEven` |
| `HPExist` | insertion | `exists_hamPath`, `hpCount_pos` (every tournament has an HP), `no_covered_universal_arc` (no tournament on `n >= 3` vertices has every arc on some HP and some arc on every HP) |
| `SwitchSum` | THM-4524 **A2** | `hpCount_complete` (`N!` orderings), `two_loop_sets`, `count_loop_sets`, `switch_sum` (`Σ_L H(switch_L T) = 2 · N!`, over the `2^N` loop sets) |
| `Redei` | **Rédei's theorem** by induction on `n`: classify the HPs of `T` by the position of the top vertex `x`; with `T_ab` = `T - x` forced `a → b`, `H(T) = #{x first} + #{x last} + Σ_{a,b} [a→x→b] c_{T_ab}(a→b)`, and mod 2 every HP of `T - x` contributes an odd amount | `redei` (`hpCount n T % 2 = 1` for every tournament), `hp_decomp`, `count_head`, `count_last`, `count_interior`; unconditional corollaries `allOdd_mod_four_of_tournament`, `allEven_odd_of_tournament`, `oddArcs_parity_of_tournament`, `exists_odd_arc_of_even`, `shave_keeps_odd_iff`, `copiesH_odd_of_even` / `containsH_of_even` (THM-4526 Theorem A, even `n`) |
| `ConstantH` | THM-4524 **A3(e)** | `constant_switching_class` (if `H` is constant on a switching class: `2^(N-1) H = N!` and `N = 2^k`), `legendre_two` (`v_2(n!) + s_2(n) = n`), `constant_class_example` (`N = 4`: every switch of `tourC 16` has `H = 3`) |
| `ParityBreak` | THM-4524 §3 parity-break number | `parity_break_two` (all arcs even, `n >= 2`: no single deletion makes `H` even, two consecutive ones do), `delete_two` (`H(T-e-f) + c(e) + c(f) = H + c(e,f)`), `out_arcs_sum` (`Σ_w c(v→w) + end(v) = H`), `end_sum`, `arc_split` |
| `Shaved` | copies of `H_n` (HP plus first→last arc) | `every_four_contains_H4`, `avoids_H5_iff` (a 5-tournament avoids `H_5` iff it is `≅ C3[1,C3,1]`) |
| `CycleCount` | rotations of closing Hamiltonian paths | `closing_count`, `copiesH_add` (`#copies(H_n) + n · hc = H`), `containsH_iff_copiesH_pos` |
| `CollatzDrop` | Syracuse map `S(A) = oddpart(3A+1)`, drop `K(A) = (A - S(A))/2`, label map `F(M) = (S(2M-1)+1)/2` | `drop_identity` (`6K+1 = (2^v-3) S`), `drop_forward`, `drop_backward`, `drop_injective`, `admissible_iff`, `admissible_of_neg`, `admissible_of_nonneg`, `mod_four_iff`, `labelF_even`, `labelF_four_j_one`, `labelF_microcosm`, `labelF_eq_sub_drop`, `drop_two_regular`, worked values |

Conventions:

* Vertices are natural numbers; a tournament is a Boolean relation `T` with exactly one of
  `T a b`, `T b a` for distinct `a, b < n` (`IsTournament`); two relations are the same
  labeled tournament when they agree on all ordered pairs of distinct vertices (`SameOn`).
* "The number of objects with property `P`" is the length of a duplicate-free list whose
  members are exactly those objects; `length_eq_of_same_members`, `hpCount_unique`,
  `arcCount_unique` show the number does not depend on the list.
* The tiling model's base path is `N-1 → ⋯ → 0` (the note's `N → ⋯ → 1`, shifted by one).
* `QR_7 - z` is relabeled `0, …, 5` in increasing order (`lift`).

## What is not proved here

* The lower bound `edim_m(Q_6) >= 15` (about 1.4e10 search leaves) and the 229-orbit
  census: out of scope. Formalized: `>= 7` (L4) and `<= 15` (the paper's set).
* The explicit sets for `Q_10, Q_11, Q_12` (sizes 48, 65, 76): the same certificate scheme
  applies, but checking keys against the full `edgeList d` costs too much kernel time
  (estimated 8 min or more for `Q_10`).
* L5; THM-4524 A3 (a)-(d), B1, C2, E1, the dead-arc formula and the census; C1 for
  non-cyclic abelian groups; THM-4526 Theorem A for odd `n` (only `n = 5` exhaustively), Theorems B, C.

## Trust boundary and kernel cost

Finite statements are reduced by the kernel (`decide +kernel`); the elaborator does not
evaluate them, and no native code is trusted. The checkers are written with the
recursors `Nat.rec` / `List.rec` and `Nat.beq` / `Nat.ble`, because structural recursion
(compiled to `brecOn`) and `DecidableEq` instance chains cost about twice the kernel
memory. Every checker is proved equal to, or sound for, its readable counterpart
(`hamR_eq`, `keyR_eq`, `keyBR_eq`, `countR_eq`, `usesArcR_eq`, `arcCountR_eq`,
`checkArcs_sound`, `hasEvenArc_sound`, `hasH_sound`, `findIso_sound`, `msortR_perm`, ...).
Large checks are split into declarations of 128 tournaments or 64 edges, and the `Q_8`/`Q_9`
certificates into several modules: kernel caches are per declaration, but the process
memory of one Lean run accumulates across declarations.

## Reproduction

```text
python3 verify.py
```

The `EdimQ7*`, `EdimQ8*`, `EdimQ9*` modules are generated by `python3 gen_certificates.py 7 8 9`
from the edim lane's explicit sets; they are trusted only through their kernel checks.

The verifier scans the sources for proof escapes, checks that every theorem appears in
`AxiomAudit.lean`, runs `lake clean`, builds the modules one at a time (at most one Lean
process at once), runs `lake env lean AxiomAudit.lean`, rejects any axiom outside
`{propext, Quot.sound}`, and writes `verification.json` and `verification.log`.
The lakefile passes `--threads=1 --memory=1400` to Lean, so a runaway kernel check aborts
instead of exceeding about 1.4 GB.

Measured on an Apple M2 (8 GB): see `verification.json` (`build_seconds`,
`peak_child_rss_mb`).
