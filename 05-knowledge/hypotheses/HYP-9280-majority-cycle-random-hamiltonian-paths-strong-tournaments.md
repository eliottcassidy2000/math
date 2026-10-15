---
id: HYP-9280
title: "Majority-cycle conjecture: in every strong tournament, the arcs whose tail usually precedes the head (in a uniform random Hamiltonian path) form a strongly connected spanning subdigraph; hence the majority tournament is strong, has at least n - 2 directed 3-cycles, and each gives a pair ordered with probability in (1/2, 2/3] - which would prove the tournament 1/3-2/3 conjecture (HYP-9277 part 1). Weakest local sufficient condition: no vertex precedes every other vertex with probability > 1/2. Sharp family: transitive plus the reversed extreme arc, value (2^(n-3) + 1)/(2^(n-2) + 1)."
status: >
  OPEN (the conjecture). PROVED: the reductions (F strong => M strong => >= n-2 cyclic majority triples => a balanced
  pair in (1/2, 2/3]; the vertex lemma at every vertex => F has a cycle => the same; "no majority source" => M not
  transitive => the same; the stochastic symmetry (SS) => the vertex lemma); Pr[w first] < 1/2 for every non-source;
  (SS) at k = 0; the vertex lemma for in-degree <= 2 and, dually, out-degree <= 2; the weak form (F has a directed
  cycle) when T has at most one inverted arc.
  FINITE-EXACT: the conjecture (F strong), the vertex lemma, D_in > 0 and the sharp T_n minima for every strong
  tournament with n <= 9 (this session's labelled sweep n <= 7 and majority sweep n = 8; the S36 audit's independent
  isomorphism-class sweep n = 8, 9). (SS) holds for all strong n <= 8 and is REFUTED at n = 9 (exactly 5 classes).
  NUMERICAL (annealing, n = 8..12): no counterexample; at most 3 inverted arcs; the sharp minima hold.
session: opus-2026-10-14-S36
source: 05-knowledge/results/strong_balance_pisano_golden_eleven_20261014.md
scripts:
  - 04-computation/experiments/tournament_majority_cycle_20261014.py
  - 04-computation/experiments/tournament_inverted_arcs_structure_20261014.py
  - 04-computation/experiments/tournament_strong_balance_rules_20261014.py
  - 04-computation/experiments/tournament_majority_search_20261014.py
  - 04-computation/experiments/tournament_vertex_lemma_20261014.py
  - 04-computation/experiments/tournament_vertex_lemma_search_20261014.py
  - 04-computation/experiments/tournament_inneighbour_symmetry_20261014.py
  - 04-computation/experiments/tournament_ss_counterexamples_n9_20261014.py
  - 04-computation/experiments/tournament_pointwise_pairing_test_20261014.py
  - 04-computation/experiments/s36_audit_20261014/ (independent audit code, n <= 9)
related: [HYP-9277, HYP-9279, THM-4655, THM-4526, THM-4527]
---

# HYP-9280: the majority cycle of random Hamiltonian paths

**Setting.** `T` is a tournament and `p(x, y)` is the probability that `x` precedes `y` in a uniformly random Hamiltonian path. `H(T)` is odd (Rédei), so `p ≠ 1/2`.
- The **majority tournament** `M(T)` orients `x → y` when `p(x, y) > 1/2`.
- An arc `x → y` of `T` is **forward** if `p(x, y) > 1/2` and **inverted** otherwise.
- `F(T)` is the digraph of forward arcs, so `F(T) = T ∩ M(T)`, and `M(T)` is `T` with its inverted arcs reversed.

**Concurrent work.** mac-mini-2026-10-14-tourn attacked the same owner question in parallel.
- THM-4655 proves the 1/3–2/3 property for every strong "transitive-plus-one" tournament, a family that contains `T_n` below, with the explicit pair `(z, u_2)`.
- HYP-9279 conjectures (A) `δ >= 4/11` for strong `n >= 4`, and (B) that every vertex of a strong tournament has a balanced partner.
- (B) is a different strengthening from this one. This session's rule R3 (`tournament_strong_balance_rules_20261014.py`) is the same statement as (B), with 0 failures for `n <= 7`; mac-mini verified it for all strong `n <= 9`.

## Conjecture

For every strong tournament `T` on `n >= 3` vertices:
1. `F(T)` is strongly connected.

Consequently:
2. `M(T)` is strong.
3. `M(T)` has at least `n − 2` directed 3-cycles (Moon).
4. Every directed 3-cycle `x → y → z → x` of `M(T)` contains a pair with `1/2 < p <= 2/3`. Indeed every linear order satisfies at most two of the three cyclic precedences, so `p(x,y) + p(y,z) + p(z,x) <= 2`.

So the conjecture would prove HYP-9277 part 1, the tournament 1/3–2/3 conjecture. The component reduction transfers it to all non-transitive tournaments.

## A ladder of sufficient conditions

Each statement below implies the next, and the last already gives a balanced pair in `(1/2, 2/3]`.
1. **(SS), stochastic symmetry at every vertex.** Let `d = indeg(w)` and let `N` be the number of in-neighbours of `w` preceding `w`. (SS) asks `Pr[N >= d − k] >= Pr[N <= k]` for all `k`. **False at `n = 9`** (below).
2. **The averaged form at every vertex.** `D_in(w) = Σ_(u→w) (2p(u, w) − 1) > 0`, i.e. `E[N] > d/2`.
3. **The vertex lemma at every vertex.** Every vertex `w` with an in-neighbour has an in-neighbour `u` with `p(u, w) > 1/2`. Then `F(T)` has minimum in-degree `>= 1`, so it contains a directed cycle.
4. **No majority source.** No vertex of a strong tournament precedes every other vertex with probability `> 1/2`. Equivalently, every `w` has some vertex `x`, an in-neighbour or not, with `p(x, w) > 1/2`.
5. **`M(T)` is not transitive.** A transitive tournament has a source, so (4) gives (5). A tournament that is not transitive has a directed 3-cycle, and item 4 of the conjecture gives the balanced pair.

**Remarks on the ladder.**
- (2) ⇒ (3) because `E[N] >= d/2` forces `max_u p(u, w) >= 1/2`, and `≠ 1/2` because `H` is odd. Parity does not upgrade `E[N] >= d/2` to `E[N] > d/2` when `d` is even. Strictness of `D_in > 0` is a FINITE-EXACT fact for `n <= 9`: its minimum is `1/H(T_n)`.
- (3) ⇒ (4): if `w` precedes everybody with probability `> 1/2`, then in particular every in-neighbour `u` has `p(u, w) < 1/2`. `w` has an in-neighbour because a strong `T` has no source.
- (4) only has to be checked at one vertex: the source of `M(T)`, if `M(T)` has one. At `T_n` it is exactly as tight as the vertex lemma. At `w = 2` the only vertex that precedes `w` with probability `> 1/2` is the in-neighbour `1`, with margin `1/H`.
- Other weaker forms, each sufficient: `F(T)` contains a directed cycle; some directed 3-cycle of `T` has no arc with `p(tail first) < 1/3`.

## Proved pieces

- **Reductions.** These are the implications above. A tournament with a strongly connected spanning subdigraph is strong, and a strong tournament has at least `n − 2` directed 3-cycles (Moon).
- **(SS) at `k = 0`: `Pr[N = d] >= Pr[N = 0]`.**
  - Take a path `w x_1 … x_(n−1)` with `w` first, and let `x_j` be the last in-neighbour of `w` in it (one exists, since `d >= 1`).
  - Move `w` to just after `x_j`. The slot is valid: `x_j → w`, and the next vertex, if any, is an out-neighbour, because `x_j` was the last in-neighbour.
  - The map is injective, since `w` is put back at the front by deleting it. Its image consists of paths with `N = d`.
- **First-vertex bound.** If `w` has an in-neighbour, then `Pr[w first] < 1/2`.
  - The `k = 0` injection maps the paths with `w` first into paths where `w` is not first. So `Pr[w first] <= 1/2`.
  - Equality is impossible because `H` is odd.
  - (The pre-audit wording used "the end is always a valid slot", which is false; the `k = 0` injection replaces it.)
- **Vertex lemma for in-degree `<= 2`.**
  - `d = 1`: `N ∈ {0, 1}` and `Pr[N = 1] >= Pr[N = 0]`, strictly by parity. So the unique in-neighbour precedes `w` with probability `> 1/2`.
  - `d = 2`: `E[N] = Pr[N >= 1] + Pr[N = 2] >= (1 − Pr[N = 0]) + Pr[N = 0] = 1`. So `p(u_1, w) + p(u_2, w) >= 1`, and one of them is `> 1/2` by parity.
  - Reversal gives the dual: a vertex of out-degree `<= 2` has an out-neighbour following it with probability `> 1/2`.
  - Consequently, in an acyclic `F(T)` every `F`-source has `T`-in-degree `>= 3` and every `F`-sink has `T`-out-degree `>= 3`. Likewise a majority source has in-degree `>= 3`.
- **At most one inverted arc (weak form only).** If `T` is strong with at most one inverted arc `e`, then `F(T) ⊇ T − e` contains a directed cycle, so `M(T)` is not transitive and `T` has a balanced pair in `(1/2, 2/3]`.
  - If `T − e` were acyclic, `T` would be a transitive tournament plus one reversed arc, and strong only if `e` reverses the extreme pair.
  - In that tournament `e` is forward, with `p = 2^(n−2)/(2^(n−2) + 1)`.
  - This gives the weak form. It does not show that `F(T)` is strongly connected, which is the conjecture itself.

## Evidence

**This session.**
- Every labelled strong tournament with `n <= 7`: `M(T)` is never transitive; `F(T)` is always strongly connected; `T` always has a 3-cycle of forward arcs.
  - Inverted arcs: none for `n <= 5`. Exactly one in 1,440 tournaments at `n = 6` and in 43,680 at `n = 7`.
  - The degree types `(outdeg y, indeg x)` of inverted arcs `x → y` are `(1, 1)` (all 1,440 at `n = 6`; 23,520 at `n = 7`), and `(2, 4)` and `(4, 2)` (10,080 each at `n = 7`).
  - Type `(1, 1)` is always a "squeezed triangle" `y → z → x`, with `z` the unique out-neighbour of `y` and the unique in-neighbour of `x`.
  - The least `p` over inverted arcs is `7/15` at `n = 6` and `9/25` at `n = 7`.
  - Script: `tournament_inverted_arcs_structure_20261014.py`.
- All 236,522,496 labelled strong tournaments at `n = 8`: `M(T)` is never transitive and a forward 3-cycle always exists. Inverted arcs: 0 in 231,227,136, 1 in 4,932,480, 2 in 362,880. Script: `tournament_majority_cycle_20261014.py`.

**The S36 audit** (independent code, canonical isomorphism-class generation matching OEIS; `s36_audit_20261014/`).
- `n = 8`, all 6,008 strong classes; `n = 9`, all 178,133 strong classes (64,026,088,576 labelled):
  - `M(T)` is never transitive, and `F(T)` is always strongly connected.
  - The vertex lemma and `D_in > 0` hold at every vertex.
  - The minimum number of 3-cycles of `M(T)` is `n − 2`.
- Inverted arcs at `n = 9`: at most 3 (4 classes).
- The `n <= 7` numbers above were reproduced, including the 10,420,928 violations of the false rule in "What fails".

**Annealing, `n = 8..12`** (`tournament_majority_search_20261014.py`). The search maximised inverted arcs and minimised 3-cycles of `M`.
- At most 3 inverted arcs were found, which matches the exhaustive `n = 9` maximum.
- The fewest 3-cycles of `M` found was `n − 2` for every `n`.
- `F` was strongly connected at every optimum.
- The HYP-9277 Claim A counterexamples (`n = 8, 9`) have a single inverted arc, of probability `13/45` and `17/81`, and `M(T)` has 10 and 12 three-cycles. The audit confirmed these.

## The sharp family `T_n` and the stochastic symmetry

**`T_n`.** `T_n` is the transitive tournament with its extreme arc reversed; it has `H = 2^(n−2) + 1`. For every `n <= 9` (FINITE-EXACT; `n = 8, 9` by the audit), the following extremes over strong `n`-tournaments are attained only by `T_n`:
- max over 3-cycles of `T` of the smallest arc probability: minimum `(2^(n−3) + 1)/(2^(n−2) + 1)`, i.e. `2/3, 3/5, …, 33/65, 65/129`;
- min over vertices of the best in-neighbour precedence (the vertex-lemma value): the same values;
- min over vertices of `D_in`, and of `D_out`: `1/H(T_n) = 1/(2^(n−2) + 1)`;
- max over vertices of `Pr[w first]`: `2^(n−3)/(2^(n−2) + 1)`.

**Not every sharp constant is at `T_n`.** Let the balance constant be `δ(T) = max over pairs of min(p, 1 − p)`, so that 1/3–2/3 means `δ >= 1/3`. Over strong tournaments its minimum is:
- `2/5` at `n = 4`, where `T_4` is the only strong class;
- `4/11` at `n = 5, 6, 7`;
- `62/157` at `n = 8, 9`.

For `n >= 5` none of these extremes is `T_n` (this session's rules script for `n <= 7`, the audit for `n = 8, 9`, and mac-mini's HYP-9279(A)).

**Conjectured sharp forms** (each equality only at `T_n`):
- every strong tournament has a 3-cycle with all arc probabilities `>= (2^(n−3) + 1)/(2^(n−2) + 1)`;
- every vertex has an in-neighbour that precedes it with at least that probability;
- `D_in(w) >= 1/(2^(n−2) + 1)` for every vertex.

The vertex-lemma annealing (`tournament_vertex_lemma_search_20261014.py`, `n = 8..12`) finds the minimum of both the vertex-lemma margin and `D_in` to be exactly `1/(2^(n−2) + 1)`, attained at `T_n` for `n = 10, 11, 12`.

**The stochastic symmetry (SS).**
- It holds at every vertex of every strong tournament with `n <= 8`: this session for `n <= 7` (11,742,416 (tournament, vertex) pairs at `n = 7`), the audit for `n = 8`.
- It **fails at `n = 9`**, in exactly 5 isomorphism classes (1,572,480 labelled tournaments).
  - Every failure is at a vertex of in-degree 7, hence out-degree 1, with `k = 3`.
  - In all five, `E[N] > 7/2`, the best in-neighbour has `p > 1/2`, and `M(T)` has 3-cycles.
  - This session's two sampled examples, `adj = [502, 508, 504, 241, 480, 448, 384, 256, 8]` and `[250, 508, 497, 500, 480, 448, 384, 256, 1]`, are the audit's classes 2 and 4. They have `E[N] = 307/85` and `419/117` and best in-neighbour `p = 53/85` and `23/39`.
  - So the vertex lemma and the conjecture survive, and only the stochastic-dominance route dies. Script: `tournament_ss_counterexamples_n9_20261014.py`; audit `verify_n9_ss_classes.py`.
- (SS) at `k <= 2` would still give the vertex lemma for in-degree `<= 6`, since `E[N] = Σ_(j>=1) Pr[N >= j]`. Every `n = 9` failure has `k = 3`.
- (SS) cannot be proved set by set. The pointwise strengthening pairs each set `X` of vertices before `w` with its complement and asks `wt(X^c) >= wt(X)` whenever `|X ∩ I| < |X^c ∩ I|`. It fails in 2,427 of 16,426 random tests (`tournament_pointwise_pairing_test_20261014.py`).

## What fails

- "Every arc has `p(tail first) >= 1/3`" is false at `n = 8` (HYP-9277 Claim A).
- "Best in-neighbour precedence `+ Pr[w first] >= 1` at every vertex" is false. It fails at 10,420,928 (tournament, vertex) pairs at `n = 7`, and in every strong tournament with `4 <= n <= 7`.

## Literature and contrast with posets

For linear extensions of a poset, the analogous relation is the linear-extension-majority (LEM) digraph.
- Its cycles are classical: Kislitsyn (1968); Fishburn (1974, 1976); Ewacha–Fishburn–Gehrlein (1990); De Loof–De Baets–De Meyer on cycle-transitivity of mutual rank probabilities.
- Gupta (arXiv:2607.23926, UNREFEREED) finds LEM cycles of length at most 8 at order 14.
- For posets the LEM relation can be transitive, so the cycle route cannot prove 1/3–2/3 there.

No prior tournament version of the majority relation was found, so the majority tournament `M(T)` of random Hamiltonian paths appears to be new. The conjecture says strong tournaments never have a transitive majority. That is the precise sense in which the tournament problem should be easier than the poset problem.
