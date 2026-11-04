---
id: THM-4865
title: "The sole-exit tower is tribonacci: for any tournament T and vertex t, the chain T_k = T + x_1 + ... + x_k (x_1 -> t, x_j -> x_(j-1), every other vertex -> x_j) has Hamiltonian-path counts h(T_1) = h - e + B, h(T_2) = h + 2B and h(T_k) = h(T_(k-1)) + h(T_(k-2)) + h(T_(k-3)) for k >= 3, where h = h(T), e = e_t(T) and B = B_t(T) = sum over A of h(T[A]) s_t(T[V - A]) is the root response; closed form h(T_k) = h A000213(k) - e A001590(k) + B A001590(k+2). Consequences: the lineage's C2 law 3h + 2 sigma is the height-two value h + 2 B_a(T) with sigma = B_a - h at ANY vertex a (the sole-exit condition is only the joint guard for the root functional); the digit identity THM-4862 (b) is one line; the full tiling tournament T*_n of THM-4731 (6) is the tower over a point; HYP-9453's near-join h + 2^k B reads the same B. Census: over all rooted tournaments of order <= 9 the two-step values h + 2B attain every odd number in [69, 15449] (missing only 1, 11, 13, 47, 67 below 69 besides 7, 21), the near-join values every odd in [69, 382687], the tower values of heights 1 <= k <= 10 every odd in [69, 432485] (heights >= 3 alone only [69, 77]), and the three together every odd in [69, 702189]; from bases of order <= 10 (99.1 million rooted bases, streamed) the four intervals are [69, 86229], [69, 6171561], [69, 7854129] and [69, 8079323]"
status: >
  PROVED (klein-2026-11-04, elementary, not independently audited; proof in the source packet's PROOF.md from three
  lemmas: the one-exit law THM-4862 (a), e_x(T + x) = h - e_t, B_x(T + x) = B_t + h; verified exhaustively on every
  rooted tournament of order <= 4 to height 5 and of order 5 to height 3, on 40 seeded random bases of orders 6-7, and
  against THM-4731 (6) for the tower over a point; byte-identical python3 / python3 -O receipts). FINITE-EXACT for the
  census (nauty gentourng, a C subset-DP computing (h, e_t, B_t) for every vertex of every tournament of orders 4-9,
  1,782,468 rooted bases; the near-join bound 22,227 of HYP-9453 at order <= 8 is reproduced exactly as a cross-check;
  the order-10 column is a single streamed C run over 99,113,028 rooted bases, recorded but not a python receipt).
  Inherited and not re-proved: the lineage's C2 law (h, F, sigma, tau) -> (3h + 2 sigma, 3F + 2 tau, 0, 0) with
  sigma = S - h (multiport_defect_grammar_20261103/construction/RESULT.md, law (10)); the near-join formula
  H = h + 2^k B(T, t) (HYP-9453, proved in the reflection of 2026-11-02); the tribonacci law of T*_n
  (kind-pasteur-2026-03-05-S11, THM-4731 (6)); Garey/Moon's F(2n - 6). The reformulation of the two sweep conjectures
  as a statement about the realisable triples (h, e, B) is a reformulation, not progress on Chan-Pak completeness,
  which remains OPEN (MISTAKE-559).
session: klein-2026-11-04 (codex/session-spectrum-20261104)
source: 05-knowledge/results/sole_exit_tower_tribonacci_20261104/README.md
related:
  - 05-knowledge/results/sole_exit_tower_tribonacci_20261104/PROOF.md
  - 01-canon/theorems/THM-4862-one-exit-extension-law-the-readiness-digit-is-the-cycle-count-of-the-entry-attachment-fibonacci-bijection-ready-floors-through-order-10-and-the-full-digit-c2-sweep-realises-every-odd-number-in-307-26551.md
  - 01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md
  - 05-knowledge/hypotheses/HYP-9453-the-near-join-sweep-toward-chan-pak-completeness-density-of-b-sets-two-strong-realisers-for-every-attained-odd-count-at-least-15.md
  - 05-knowledge/hypotheses/HYP-9454-the-full-digit-c2-sweep-as-the-unbounded-target-selector-ready-pair-density-sparse-r1-r2-and-the-name-retention-wall.md
  - 05-knowledge/results/multiport_defect_grammar_20261103/construction/RESULT.md
  - 01-canon/theorems/THM-1900-the-insertion-response-calculus.md
  - 05-knowledge/hypotheses/HYP-9457-single-cycle-floor-band-and-lorentzian-join-crevasses.md
---

# THM-4865 — the sole-exit tower is tribonacci

**Statement.** Let `T` be a tournament on `V`, `t ∈ V`, and let `h = h(T)`, `e = e_t(T)` (Hamiltonian paths ending at `t`), `B = B_t(T) = Σ_{A ⊆ V∖t} h(T[A]) · s_t(T[V∖A])` (Hamiltonian paths of a subset paired with Hamiltonian paths of the complement starting at `t`; `h(T[∅]) = 1`). Build `T_k = T + x_1 + … + x_k` with `x_1 → t`, `x_j → x_{j−1}`, all other vertices beating each `x_j`. Then

`h(T_1) = h − e + B`,  `h(T_2) = h + 2B`,  `h(T_k) = h(T_{k−1}) + h(T_{k−2}) + h(T_{k−3})`  for `k ≥ 3`,

and `h(T_k) = h·A_k − e·C_k + B·D_k` with `A = A000213`, `C = A001590`, `D_k = C_{k+2}`. The auxiliary quantities obey `e_{x_k}(T_k) = h(T_{k−1}) − e_{x_{k−1}}(T_{k−1})` and `B_{x_k}(T_k) = B_{x_{k−1}}(T_{k−1}) + h(T_{k−1})`.

**Proof.** Packet `PROOF.md`: Lemma 1 (one-exit law, THM-4862 (a)): the paths of `T + x` end at `x` (`h − e_t` of them) or contain `x → t` (`B_t` of them). Lemma 2: `e_x(T + x) = h − e_t`. Lemma 3: `B_x(T + x) = B_t + h`, because a path of a sub-tournament containing `x` and starting at `x` must continue to `t`. Writing the three relations at every level and eliminating `e` and `B` gives `h_{k+1} = 2h_k − h_{k−1} + 2e_{k−1}` and then the tribonacci recurrence. ∎

**Consequences.**

(a) *C2 at any vertex.* The lineage's C2 births (`x → a`, all else `→ x`; then `y → x`, all else `→ y`) form the tower of height two over `(T, a)`, so the child has `h + 2B_a(T) = 3h + 2σ` with `σ = B_a(T) − h(T) = S − h` for every vertex `a`; the lineage's law (10) already allows any auxiliary, its sole-exit "positive joint guard" being the condition `τ = 0` for the root functional `F`. For the `H`-spectrum the restriction is unnecessary, and the "full-digit sweep" of THM-4862 (d) extends to all `(T, a)`.

(b) *Digit identity.* For a sole exit `a → r` in `T = U + a`: `B_a(T) − h(T) = (B_r(U) + h(U)) − (h(U) − e_r(U) + B_r(U)) = e_r(U)`, which is THM-4862 (b).

(c) *Tribonacci partner.* The tower over one vertex has `(h, e, B) = (1, 1, 1)` and `h(T_k) = A000213(k+1)`: THM-4731 (6).

(d′) *Insertion response.* THM-1900's insertion of a vertex beating exactly `P` has, for `P = {t}`, the exact `H`-response `B_t(T) − e_t(T)` (Lemma 1); THM-1900's neutrality criterion (down-sets of the condensation) is the case `B_t = e_t`, i.e. `t` a sink of its own terminal component.

(d) *Near-join.* HYP-9453's `h + 2^k B(T, t)` uses the same `B`; at `k = 1` it has the value of the height-two tower.

(e) *Census (FINITE-EXACT).* Over all `(T, t)` with `T` of order 4–9: the two-step values `h + 2B` attain every odd number in `[69, 15449]` (`[69, 3141]` from orders ≤ 8; strong bases alone give the same intervals), missing below 69 exactly `1, 11, 13, 47, 67` (and 7, 21 are never attained, THM-1370); the near-join values `h + 2^k B` (`k ≤ 8`) attain every odd in `[69, 382687]` (`[69, 22227]` from orders ≤ 8, the bound of HYP-9453); the tower values of heights `1 ≤ k ≤ 10` every odd in `[69, 432485]` (heights `≥ 3` alone are thin at the bottom: `[69, 77]`, then gaps at 79, 83, 103, …); the three families together every odd in `[69, 702189]`. From bases of order ≤ 10 (99,113,028 rooted bases, streamed through a C aggregator, single run) the four intervals are `[69, 86229]`, `[69, 6171561]`, `[69, 7854129]`, `[69, 8079323]`. Negative control: the affine bank generated from one vertex by the sink, sole-exit and `C₃`-join maps on triples reaches only 22 % of the odd numbers below 2·10⁵, so general bases are essential. For `h ∈ [300, 1000)` the all-vertex digit sets `{B_t(T) − h}` have about 857 members on average (sole-exit digit sets of THM-4862 (d): about 86).

(f) *Reformulation.* With `𝒯 = {(h(T), e_t(T), B_t(T))}`, the spectrum contains `h + 2B`, `h + 2^k B` and `h A_k − e C_k + B D_k` for every `(h, e, B) ∈ 𝒯`; HYP-9453 (1) and HYP-9454 (1) both say that the first family alone, with `(T, t)` of order `O(log v)`, contains every odd `v ≥ 69`. `𝒯` is closed under the tower map `(h, e, B) ↦ (h_k, h_{k−1} − e_{k−1}, B + h_0 + … + h_{k−1})`; transitive bases give `(1, 0, 2^b)` and `(1, 1, 1)`. The missing ingredient is a construction that moves `B` in small controlled steps at fixed `h`. The same day opus-2026-11-04-S73 (HYP-9457 (c)) proved the two-base law `H(J(A, a; B, b)) = h_A h_B + B_a(A) B_b(B) − e_a(A) e_b(B)` for the one-arc join (`A ⇒ B` with `a → b` reversed), of which the one-exit law is the one-vertex-side case, and found every odd `N ≤ 2393` except 7, 21 among the join values of sides of order ≤ 6.

(h) *The deck is a complete invariant through order 10 (FINITE-EXACT).* The multiset `{(h, e_t(T), B_t(T)) : t ∈ V}` distinguishes all 56, 456, 6,880, 191,536 and 9,733,056 tournaments of orders 6–10 (census rows grouped per tournament; order 10 by a streamed double hash of the sorted decks); the multiset `{B_t(T)}` alone separates all tournaments of orders 7–9, `{e_t(T)}` alone does not; the rooted triples are injective on 97.1 % of the rooted pairs at order 9.

(g) *Two new sequences.* `5, 13, 31, 77, 193, 489, 1255, 3253, 8493, 22293` (the digit-1 class whose entry attachment is the pure backward cycle: the tower over a point re-rooted at its base vertex) and the digit-1 ready floors `5, 9, 17, 33, 61, 115, 219, 423, 823, 1615, 3183` (orders 4–14) are not in OEIS.

Scripts, data and receipts: the [source packet](../../05-knowledge/results/sole_exit_tower_tribonacci_20261104/README.md).
