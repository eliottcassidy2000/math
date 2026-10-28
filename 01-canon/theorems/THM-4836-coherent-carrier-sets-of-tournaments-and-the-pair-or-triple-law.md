---
id: THM-4836
title: "Coherent carrier sets and the pair-or-triple law. For a tournament T on n >= 4 vertices let Sigma(T), a subset of {L, C, S}, record which of the three S3 coherent carriers T contains: L a directed Hamiltonian path, C a directed (n-1)-cycle (plus the remaining vertex), S a source or a sink. Then L is always present (Redei), C iff T is strong or T = (strong (n-1)-set) + source/sink (Camion-Moon), and S iff T has a source or sink. The classes {L,C}, {L,C,S}, {L,S}, {L} have s_n, 2s_(n-1), 2t_(n-1) - t_(n-2) - 2s_(n-1) and [x^n] S_(>=3)(x)^2 T(x) types (s strong types, t all types). Hence every tournament type maps to a PAIR or a TRIPLE of carrier kinds exactly when n <= 5. The first singletons are C3 => C3 (n = 6), then 3 types at n = 7 and 16 at n = 8: precisely the non-strong tournaments whose first and last strong components have >= 3 vertices."
status: >
  PROVED from cited classical theorems: Rédei (every tournament has a Hamiltonian path); Camion and Moon (every
  strong tournament on m >= 3 vertices is Hamiltonian and vertex-pancyclic). The counts come from the unique ordered
  decomposition into strong components. FINITE-EXACT confirmation: every type for 4 <= n <= 7
  (coherent_carrier_sets.py; pre-registered P3, pushed in efdf0b0a5f before the run, PASSED exactly). Not yet
  audited.
  The n = 4 dictionary (TT4 -> {L,S}, strong -> {L,C}, C3 +- source/sink -> {L,C,S}; tetrahedral duality acting as
  complementation C <-> S on these sets) is from the S3 note, restated.
  Concurrent work: session81 (GRAPH_CARRIERS.md) reads the owner's "pair or triple" as an encoding-capacity question
  for independent (n-1)-edge graph types, which fails at n = 3, 9, 14. This theorem is a different map: a tournament
  type goes to the set of coherent carrier kinds it contains.
session: kind-pasteur-2026-10-27-S8
source: 05-knowledge/results/marked_carriers_euler_characteristic_tiers_20261027.md
scripts:
  - 04-computation/tournament_carriers/coherent_carrier_sets.py (+ .out, .json)
related: [THM-4835, HYP-9435, trijection_four_nodes_three_edges_20261020, session80_tournament_symmetry_20261027, session81_counterexample_families_20261028, Redei 1934, Camion 1959, Moon 1966]
---

# THM-4836 — Coherent carrier sets and the pair-or-triple law

## Statement

Let T be a tournament on `n >= 4` vertices. Write `s_m` and `t_m` for the numbers of strong and of all tournament types
on m vertices, with `t_0 = 1`.

1. **L, always.** `L ∈ Σ(T)` always.
2. **C.** `C ∈ Σ(T)` iff one of the following holds:
   - T is strong;
   - `T = (strong tournament on n-1 vertices) ⇒ v` or `v ⇒ (strong tournament on n-1 vertices)`.
3. **S.** `S ∈ Σ(T)` iff T has a source or a sink.
4. **Class counts.** The numbers of types are:

   | `Σ(T)` | number of types |
   |---|---|
   | `{L,C}` | `s_n` |
   | `{L,C,S}` | `2 s_{n-1}` |
   | `{L,S}` | `2 t_{n-1} - t_{n-2} - 2 s_{n-1}` |
   | `{L}` | `t_n - s_n - 2 t_{n-1} + t_{n-2} = [x^n] S_{>=3}(x)^2 T(x)` |

   Here `S_{>=3}(x) = Σ_{m>=3} s_m x^m` and `T(x) = Σ_{m>=0} t_m x^m`.
5. **Singletons.** `Σ(T) = {L}` iff T is not strong and its first and last strong components (in the condensation
   order) both have at least 3 vertices.
6. **Pair-or-triple law.** Every type has `|Σ(T)| >= 2` iff `n <= 5`. The singleton counts are 1, 3, 16 at
   `n = 6, 7, 8`.

## Proof

1. Rédei.
2. A directed `(n-1)`-cycle omitting v makes `T - v` Hamiltonian, hence strong.
   - If T is strong and `n >= 4`, Moon's vertex-pancyclicity gives cycles of length `n-1`.
   - If T is not strong and `T - v` is strong for some v, the condensation of T has exactly the two components
     `T - v` and `{v}`.
3. Definition.
4. The counts follow from the unique decomposition of T as an ordered sum `T_1 ⇒ ... ⇒ T_k` of strong components.
   - T has a source iff `|T_1| = 1`, and a sink iff `|T_k| = 1`. So the number with a source or a sink is
     `2 t_{n-1} - t_{n-2}`.
   - The types `{L,C,S}` are `1 ⇒ (strong n-1)` and `(strong n-1) ⇒ 1`, giving `2 s_{n-1}` (n ≥ 4, so these two
     are distinct).
   - The remaining source/sink types are `{L,S}`.
   - The strong types are `{L,C}`.
   - The rest have `k >= 2` with `|T_1|, |T_k| >= 3`. They are counted by `S_{>=3}(x) · T(x) · S_{>=3}(x)` (any
     middle ordered sum), which is part 5.
5. Part 5 is the complement of parts 2–3 among non-strong T.
6. For `n <= 5` a non-strong T cannot have both end components of size `>= 3`. Then use `S_{>=3}^2 = x^6 + 2x^7 + 13x^8 + ...`
   and `t = 1, 1, 1, 2, ...`. `□`

## Data (FINITE-EXACT, all types)

| n | {L,C} | {L,C,S} | {L,S} | {L} | t_n |
|---:|---:|---:|---:|---:|---:|
| 4 | 1 | 2 | 1 | 0 | 4 |
| 5 | 6 | 2 | 4 | 0 | 12 |
| 6 | 35 | 12 | 8 | 1 | 56 |
| 7 | 353 | 70 | 30 | 3 | 456 |
| 8 (formula) | 6008 | 706 | 150 | 16 | 6880 |

## Reading

- **Duality at n = 4.** Tetrahedral duality swaps `TT4 ↔ strong` and acts on carrier sets as the complementation of
  carriers `K_3 + K_1 ↔ K_{1,3}` (`C ↔ S`), fixing `P_4` (L) (S3 Theorems 2-3).
- **For n ≥ 6.** The law fails exactly where the composition calculus of session80 (`A = S + Σ v_k S^k`, ordered sums
  of strong components) produces a two-sided strong padding: a tournament whose two ends are strong of size ≥ 3 has
  neither a source/sink nor a strong `(n-1)`-subtournament.
- **The coarsest member of the ladder.** `Σ(T)` is the coarsest member of the marked-carrier ladder of THM-4835. It
  records only whether the marked types over T are nonempty.
