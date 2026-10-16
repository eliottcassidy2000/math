---
id: HYP-9280
title: "Majority-cycle conjecture: in every strong tournament, the arcs whose tail usually precedes the head (in a uniform random Hamiltonian path) form a strongly connected spanning subdigraph; hence the majority tournament is strong, has at least n - 2 directed 3-cycles, and each gives a pair ordered with probability in (1/2, 2/3] - which would prove the tournament 1/3-2/3 conjecture (HYP-9277 part 1). Weakest local sufficient condition: no vertex precedes every other vertex with probability > 1/2. Sharp family: transitive plus the reversed extreme arc, value (2^(n-3) + 1)/(2^(n-2) + 1)."
status: >
  OPEN (the conjecture). PROVED: the reductions (F strong => M strong => >= n-2 cyclic majority triples => a balanced
  pair in (1/2, 2/3]; the vertex lemma at every vertex => F has a cycle => the same; "no majority source" => M not
  transitive => the same; the stochastic symmetry (SS) => the vertex lemma); Pr[w first] < 1/2 for every non-source;
  (SS) at k = 0; the vertex lemma for in-degree <= 2 and, dually, out-degree <= 2; the weak form (F has a directed
  cycle) when T has at most one inverted arc.
  PROVED (post-audit): at a vertex w of out-degree 1 whose out-neighbour o beats every vertex of U = V - {w, o}
  (equivalently indeg(o) = 1), C(u, w) - C(w, u) = H(T[U]) for EVERY in-neighbour u (so all are witnesses).
  FINITE-EXACT: the conjecture (F strong), the vertex lemma and the sharp T_n minima for every strong tournament with
  n <= 9 (this session's labelled sweep n <= 7 and majority sweep n = 8; the S36 audit's independent isomorphism-class
  sweep n = 8, 9). For every strong tournament of feedback arc number <= 2 (n <= 16), <= 3 (n <= 14), <= 4 (n <= 12),
  <= 5 (n = 10) (post-audit sweep; T_n, the only strong feedback-1 tournament, checked separately): F strong, the vertex
  lemma, no majority source, and the T_n minima of the vertex-lemma value, the 3-cycle certificate and Pr[w first] (the
  last two by the second audit); the D_in and D_out minima FAIL there from n = 11 (D_in < 0). REFUTED: (SS) (n = 9, exactly 5 classes) and the averaged form D_in(w) > 0
  (n = 11: the double chain D(11, 5) has E[N] - d/2 = -23/578 at its out-degree-1 vertex; both verified by explicit
  path enumeration). OPEN with no failure found: (SS) at k = 1, which would give the vertex lemma for in-degree <= 4.
  NUMERICAL (annealing, n = 8..12): weak - the D_in annealing missed the n = 11, 12 failures (its n = 8, 9 values are
  the exact minima). Second audit of the continuation: MISTAKE-617.
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
  - 04-computation/experiments/tournament_near_transitive_sweep_20261014.py
  - 04-computation/experiments/tournament_double_chain_dfs_check_20261014.py
  - 04-computation/experiments/tournament_outdegree_one_vertex_lemma_20261014.py
  - 04-computation/experiments/tournament_witness_rules_20261014.py
  - 04-computation/experiments/tournament_ss_k1_refinement_test_20261014.py
  - 04-computation/experiments/tournament_ss_k1_local_move_matching_20261014.py
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
2. **The averaged form at every vertex.** `D_in(w) = Σ_(u→w) (2p(u, w) − 1) > 0`, i.e. `E[N] > d/2`. **False at `n = 11`** (post-audit section).
3. **The vertex lemma at every vertex.** Every vertex `w` with an in-neighbour has an in-neighbour `u` with `p(u, w) > 1/2`. Then `F(T)` has minimum in-degree `>= 1`, so it contains a directed cycle.
4. **No majority source.** No vertex of a strong tournament precedes every other vertex with probability `> 1/2`. Equivalently, every `w` has some vertex `x`, an in-neighbour or not, with `p(x, w) > 1/2`.
5. **`M(T)` is not transitive.** A transitive tournament has a source, so (4) gives (5). A tournament that is not transitive has a directed 3-cycle, and item 4 of the conjecture gives the balanced pair.

**Remarks on the ladder.**
- (2) ⇒ (3) because `E[N] >= d/2` forces `max_u p(u, w) >= 1/2`, and `≠ 1/2` because `H` is odd. Parity does not upgrade `E[N] >= d/2` to `E[N] > d/2` when `d` is even.
- Rungs (1) and (2) are false in general, but the implications stand. The live statements are (3), (4) and (5).
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
- every vertex has an in-neighbour that precedes it with at least that probability.

Both hold in every sweep below. The third form proposed before the post-audit sweep, `D_in(w) >= 1/(2^(n−2) + 1)`, is REFUTED at `n = 11`. The vertex-lemma annealing (`tournament_vertex_lemma_search_20261014.py`, `n = 8..12`) had reported the `D_in` minimum as `1/(2^(n−2) + 1)`. At `n = 8, 9` that is the exact minimum (FINITE-EXACT, above), and at `n = 10` nothing refutes it. At `n = 11, 12` it was only the best value the search found, not the minimum: the search missed the double chains.

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

## Post-audit continuation (same session)

**The averaged form fails (REFUTED, FINITE-EXACT).**
- **The family.** The double chain `D(n, j)` is the transitive order on positions `0..n−1` (higher beats lower) with the arcs `0 → j` and `j → n−1` reversed. It is strong for `1 <= j <= n−2`.
- **The failing vertex.** At `w = 0`, which has out-degree 1 (out-neighbour `j`), the top vertex `n − 1` precedes `w` only rarely.
  - In `D(11, 5)`, `p(10 before 0) = 33/289`. The other in-neighbours have `p = 9/17` (`u = 1..4`) or `161/289` (`u = 6..9`).
  - So `E[N] − d/2 = −23/578 < 0`, while the vertex lemma holds with `p = 161/289`.
  - `D(12, 6)` gives `E[N] − d/2 = −91/561`.
- **Verification.** These values come from explicit enumeration of all Hamiltonian paths (`tournament_double_chain_dfs_check_20261014.py`) and agree with the subset DP.
- **Growth.** Over the double chains, `min D_in` decreases to `−0.864` at `n = 18`. The vertex lemma, "no majority source" and HYP-9280 hold throughout, with the minimum vertex-lemma value at `j = 1` (`T_n`).
- **Where it fails.** In the sweeps below, `D_in <= 0` first appears at `n = 11`, and always at a vertex of in-degree `>= 9` (second audit).
  - `E[N] − d/2 > 0` at every vertex of in-degree `d <= 6` in every sweep.
  - For `3 <= d <= 6` the minimum is `0.0471`. The `d = 3` minimum falls with `n`, from `0.1796` at `n = 8` to `0.0731` at `n = 14` (feedback 3), so there is no uniform slack even there.
  - For `d = 1, 2` it tends to 0, as low as `0.0001` and `0.0057` at `n = 16`. Those degrees are proved.
- **Status at `n = 10`.** `D_in > 0` holds for every strong tournament with `n <= 9`. At `n = 10` it is unknown outside feedback arc number `<= 5`.

**Exhaustive near-transitive evidence (FINITE-EXACT; `tournament_near_transitive_sweep_20261014.py`).**
- **Coverage.** A tournament with feedback arc number `k` is a transitive order with `k` arcs reversed. In an optimal order no reversed arc joins positions at distance 1, since swapping two adjacent positions would lower `k`. So running over all `k`-sets of arcs `lo → hi` with `hi − lo >= 2` (and the smaller `k`) covers every tournament of feedback arc number `<= k`. `T_n`, the only strong tournament of feedback arc number 1, lies outside these runs and is covered by separate checks (`n <= 20`).
- **Ranges.** Every strong tournament of feedback arc number:
  - `<= 2` for `n <= 16`;
  - `<= 3` for `n <= 14`;
  - `<= 4` for `n <= 12`;
  - `<= 5` for `n = 10`.
- **Results.**
  - `M(T)` is never transitive and `F(T)` is always strongly connected.
  - The vertex lemma and "no majority source" hold at every vertex.
  - (SS) never fails at `k = 0, 1`. It fails at `k = 2` from `n = 10` on.
  - `D_in <= 0` occurs from `n = 11` on, in up to 337 tournaments (`n = 14`, feedback 3).
- **No uniform slack.** The minimum vertex-lemma margin decays like `2^(−n)` even at vertices with `min(indeg, outdeg) >= 3` (`0.00041` at `n = 16`). So no argument with uniform slack can work; the proof must be exact, as the injections for degree `<= 2` are.

**Out-degree-1 vertices.** Let `w` have the unique out-neighbour `o`. Let `U` be the other vertices (all in-neighbours of `w`), `U+ = {u ∈ U : o → u}` and `U− = {u ∈ U : u → o}`.
- **PROVED: if `U−` is empty (`o` beats every vertex of `U`, i.e. `indeg(o) = 1`), then `C(u, w) − C(w, u) = H(T[U])` for every `u ∈ U`.**
  - *Proof.* A Hamiltonian path either ends at `w` or contains `w o`.
  - If `w` is last, the path is `Q w` with `Q` a Hamiltonian path of `T[U ∪ {o}]`. There `o` has in-degree 0, so `Q = o Q'` with `Q'` any Hamiltonian path of `T[U]`. That gives `H(T[U])` paths, in all of which `u` precedes `w`.
  - Otherwise the path is `A w o B'` for an arbitrary ordered pair of paths `(A, B')` partitioning `U`, possibly empty. Every such pair is valid, because `U ⇒ w` and `o ⇒ U`. The swap `(A, B') ↦ (B', A)` is an involution that exchanges "`u` before `w`" and "`w` before `u`". ∎
  - So `p(u before w) = 1/2 + H(T[U])/(2H(T))` for every in-neighbour. For `T_n` at its out-degree-1 vertex this is the extremal value `(2^(n−3) + 1)/(2^(n−2) + 1)`.
- **FINITE-EXACT, no failure.** The session script (`tournament_outdegree_one_vertex_lemma_20261014.py`) ran labelled `n <= 7`, feedback 2 at even `n = 10..16`, feedback 3 at `n = 10..12`, and feedback 4 at `n = 10`. The second audit extended this to exhaustive `n = 8` (all 1,866,256 strong configurations), odd-`n` near-transitive families and random `n = 6..12`. In all of these:
  - a witness always exists in `U+`;
  - when `U−` is nonempty, the average of `2p − 1` over `U−` is positive, and `U−` contains a witness.
- **The `U+` average is not the mechanism.** It can be strictly negative: never for labelled `n <= 7`, but in 2,160 of the 1,866,256 configurations at `n = 8`, up to 395 pairs in the session's sweeps, and up to 1,125 in the audit's (`(11, 4)`, with one case exactly 0).
- **Exact form of the `U−` sum** (second audit, verified per path on all 1,833,488 exhaustive `n = 8` cases): `Σ_(u∈U−)(C(u,w) − C(w,u)) = |U−|·N_1 + S_N`. Here `N_1` is the number of paths ending at `w`, `S_N` is the signed `U−` count over paths that start in `U−`, and the swappable class contributes 0. `S_N` alone can be negative, and the minimum `U−` average tends to 0 (`0.0117` at `n = 16`).
- **Conjecture.** For out-degree-1 vertices, the average of `2p(u before w) − 1` over `U−` is positive. This is OPEN; the swap involution reduces it to a comparison of paths ending in `U−` with paths starting in `U−` inside sub-tournaments.

**Which in-neighbour is the witness? Natural rules fail** (`tournament_witness_rules_20261014.py`; refinement: `tournament_ss_k1_refinement_test_20261014.py`):
- the predecessor-weighted average `Σ_u Pr[pred(w) = u]·(2p(u, w) − 1) > 0` first fails, among the tournaments tested (no failure at `n <= 9` in the families checked), at `n = 10`, for the arcs `(0, 4), (1, 9)`: at `w = 1` its value is `−656/167281`, while the vertex lemma holds there (`u = 4`, `p = 213/409`);
- "minimum out-degree", "minimum score in `T[U]`", "maximum score in `T[U]`" and "the most frequent predecessor" each fail by `n <= 7`;
- "beats the fewest out-neighbours of `w`" survives labelled `n <= 7` and feedback 2 to `n = 14` under the most favourable tie-break (a tied choice can fail already at `n = 6`). It fails outright at `n = 8` with feedback 3: reversed arcs `(0,5), (1,6), (5,7)`, `w = 1`, unique choice `u = 5` with `p = 88/179` (second audit).

The per-in-neighbour refinement of (SS) at `k = 1` also fails, in 41.3% of the triples `(T, w, u)` with in-degree `d >= 3` at `n = 7` (11,734,800 of 28,413,840). It would ask that, for each `u`, the paths with `u` the only in-neighbour before `w` be no more numerous than those with `u` the only in-neighbour after `w`. These failures rule out the simplest local witnesses and refinements; they do not exclude every local argument.

**Next target.** (SS) at `k = 1`: `Pr[N >= d − 1] >= Pr[N <= 1]`.
- It held at every vertex of every tournament tested.
- With the proved `k = 0` case it gives `E[N] >= 1 + Pr[N >= 2] >= 3/2` for `d = 3`, and `E[N] >= 2` for `d = 4` (exact LP check by the second audit). So it would give the vertex lemma for in-degree `<= 4`. The reversed statement gives, at a vertex of out-degree `<= 4`, an out-neighbour that follows it with `p > 1/2`. That is not the vertex lemma at such vertices, which is where every known (SS) and `D_in` failure sits.
- NUMERICAL (seeded random sample): unlike `k = 0`, it is usually not realised by the two local move sets tested. A matching from `{N <= 1}` into `{N >= d − 1}` along single-vertex relocations fails for 717 of 740 (tournament, vertex) pairs at `n = 8`, and with two-vertex blocks for 504. Matchings do exist for 23 and 236 of them, and most failures are trivial: some path with `N <= 1` has no relocation into the target at all (`tournament_ss_k1_local_move_matching_20261014.py`).
