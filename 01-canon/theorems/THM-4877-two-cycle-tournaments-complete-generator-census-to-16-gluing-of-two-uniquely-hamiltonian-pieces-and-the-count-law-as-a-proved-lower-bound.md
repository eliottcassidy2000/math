---
id: THM-4877
title: "Tournaments with exactly two Hamiltonian cycles. (1) Complete generator: for every K >= 1 and n >= 5, every tournament of order n with hc <= K is a one-vertex extension of one of order n - 1 with 1 <= hc <= K by a vertex whose out-neighbourhood is a union of at most K cyclic intervals of a Hamiltonian cycle of the parent, with hc(child) = sum_{w in O, u not in O} P(w -> u) (Moon 1966 + the extension lemma); hence exhaustive censuses: two-cycle counts 5, 18, 60, 191, 589, 1775, 5257, 15360, 44394, 127171, 361595 at orders 6-16 (= (n-2)F(2n-9) - L(2n-10) at every order), three-cycle counts through order 13, four-cycle through 12. (2) Gluing theorem: if A is uniquely Hamiltonian with cycle s a_1 ... a_k s and s a sole entry, B uniquely Hamiltonian with cycle s b_1 ... b_l s and s a sole exit, and the two are identified at s with every cross arc b_j -> a_i except a_k -> b_1, the result has exactly two Hamiltonian cycles; likewise along a common arc s -> s' of two consecutive sole entries / sole exits. (3) Distinct rooted pairs give distinct gluings, so with THM-4868 (F) the two-cycle count satisfies D_n >= sum over k + l = n - 1 of F(2k-3) F(2l-3) + the same sum over k + l = n - 2 = (n-2)F(2n-9) - L(2n-10) for all n >= 6, with equality exactly when every two-cycle tournament of order n is a gluing: exhaustive through order 16, directly verified through order 12. Also FINITE-EXACT: #cut arcs + #non-critical vertices = n for every two-cycle tournament through order 14 (n + 1 is THM-4868's uniquely Hamiltonian identity), one or two sole exits, one or two unique-cycle deletions."
status: >
  PROVED (klein-2026-11-04, elementary, not independently audited): (1) the completeness lemma, from Moon's deletion
  theorem (1966) and the extension lemma of opus-2026-11-04-S74 section 2 (its interval part is the block lemma of the
  sweep packet, PROOFS.md section 6); (2) the gluing theorem for both junction types; (3) injectivity and the lower bound,
  the rooted counts F(2m - 5) (sole-exit-rooted classes of UH(m)) and F(2m - 7) (classes of UH(m) with two sole exits,
  which are consecutive) being THM-4868 (F), (C) (kind-pasteur-2026-11-04-S10, PROVED), and the generating-function
  identity checked to n < 60 and provable by hand. FINITE-EXACT: the censuses (nauty labelg canonical forms; the
  generator reproduces Garey's F(2n - 6) at every order 6-16 and HYP-9459's exhaustive k = 3, 4 columns through 10);
  the converse of (2) (every two-cycle tournament is a gluing) through order 16 by the count equality and through
  order 12 class by class; the identity #cut + #NC = n and the degree facts through order 14; HYP-9459 (2), (3)
  confirmed at orders 11-13 / 11-12. OPEN: the converse for all n (equivalent to the count law for all n; it implies
  HYP-9459 (1)). Receipts byte-identical python3 / python3 -O; generator runs recorded in logs.
session: klein-2026-11-04 (codex/session-spectrum-20261104)
source: 05-knowledge/results/two_cycle_tournaments_gluing_20261104/README.md
related:
  - 05-knowledge/hypotheses/HYP-9459-few-cycle-tournaments-by-minimum-out-degree-f-2n-9-and-2f-2n-11-and-the-out-degree-1-lemma-for-one-and-two-cycles.md
  - 01-canon/theorems/THM-4862-one-exit-extension-law-the-readiness-digit-is-the-cycle-count-of-the-entry-attachment-fibonacci-bijection-ready-floors-through-order-10-and-the-full-digit-c2-sweep-realises-every-odd-number-in-307-26551.md
  - 01-canon/theorems/THM-4868-the-noncritical-run-of-a-uniquely-hamiltonian-tournament-local-roundness-one-run-cut-arcs-plus-noncritical-equals-n-plus-1.md
  - 01-canon/theorems/THM-4867-every-uniquely-hamiltonian-tournament-is-a-one-arc-join-shelor-1977.md
  - 05-knowledge/results/uniquely_hamiltonian_one_arc_cut_shelor_20261104.md
  - 05-knowledge/results/two_branch_rearming_sigma_sweep_20261104/README.md
---

# THM-4877 — two-cycle tournaments: complete generator, census to order 16, gluing, and the count law as a proved lower bound

Notation: `hc` = number of Hamiltonian cycles; `UH(m)` = uniquely Hamiltonian tournaments of order `m`; `D_n` = number of tournaments of order `n` with exactly two Hamiltonian cycles; sole exit / sole entry = vertex of out-/in-degree 1; `P_T(w → u)` = number of Hamiltonian paths of `T` from `w` to `u`.

**(1) Completeness lemma.** Let `T` have order `n ≥ 5` and `hc(T) = h ≥ 1`. Then some vertex `v` has `T′ = T − v` strong with `1 ≤ hc(T′) ≤ h`, and the out-neighbourhood `O` of `v` meets every Hamiltonian cycle of `T′` in at most `h` cyclic intervals; conversely `hc(T′ + v) = Σ_{w ∈ O, u ∉ O} P_{T′}(w → u)` for any strong `T′` and proper nonempty `O`.

*Proof.* `T` is strong, so Moon (1966) gives `v` with `T − v` strong. Each Hamiltonian cycle `C′` of `T′` has an arc `u → w` with `u ∉ O`, `w ∈ O`, and `u → v → w` inserted there is a Hamiltonian cycle of `T`; distinct such switches and distinct `C′` give distinct cycles of `T`, so `hc(T′) ≤ h` and `C′` carries at most `h` switches, i.e. `O ∩ C′` is at most `h` intervals. A Hamiltonian cycle of `T′ + v` is `v → w ⇝ u → v` with `w ⇝ u` a Hamiltonian path of `T′`. ∎

So the tournaments with `hc ≤ K` of order `n` are exactly the extensions, filtered by `hc ≤ K`, of those of order `n − 1` with `1 ≤ hc ≤ K` by out-neighbourhoods with at most `2K` in/out switches along one Hamiltonian cycle of the parent. The packet's generator (`tcgen.c`, `tcgen2.c`, path-count matrix per parent, nauty `labelg` canonical forms) implements this; its censuses are therefore exhaustive:

| order | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `hc = 1` | 8 | 21 | 55 | 144 | 377 | 987 | 2584 | 6765 | 17711 | 46368 | 121393 |
| `hc = 2` | 5 | 18 | 60 | 191 | 589 | 1775 | 5257 | 15360 | 44394 | 127171 | 361595 |
| `hc = 3` | 11 | 38 | 124 | 390 | 1193 | 3574 | 10537 | 30676 | | | |
| `hc = 4` | 6 | 25 | 91 | 315 | 1049 | 3390 | 10705 | | | | |

`hc = 1` is Garey's `F(2n − 6)` at every order; `hc = 2` equals `(n − 2)F(2n − 9) − L(2n − 10)` at every order `6 ≤ n ≤ 16`. Every two-cycle class of orders 6–16 has a sole exit and a sole entry (HYP-9459 (1) through 16, previously through 10); `u_3(n) = F(2n − 9)` holds at orders 11–13 (233, 610, 1597) and `u_4(n) = 2F(2n − 11)` at 11–12 (178, 466), as HYP-9459 (2), (3) predicted.

**(2) Gluing theorem.** Let `A ∈ UH(k + 1)` have cycle `s a₁ … a_k s` with `s` a sole entry of `A`, and `B ∈ UH(l + 1)` have cycle `s b₁ … b_l s` with `s` a sole exit of `B` (`k, l ≥ 2`). Let `Y = A ⊕_s B` be the tournament on `V(A) ∪ V(B)` (the two `s` identified) with the arcs of `A` and `B` and `b_j → a_i` for all `i, j` except `a_k → b₁`. Then `hc(Y) = 2`, the cycles being `s a₁ … a_k b₁ … b_l` and `a₁ … a_k s b₁ … b_l`. The same holds for the arc junction: `A ∈ UH(k + 2)` with cycle `s s′ a₁ … a_k s` and `s, s′` its two sole entries, `B ∈ UH(l + 2)` with cycle `s s′ b₁ … b_l s` and `s, s′` its two sole exits, identified along `s → s′`, cross arcs `b_j → a_i` except `a_k → b₁`.

*Proof.* Put `X = {s, a₁, …, a_k}`, `Z = {b₁, …, b_l}`. The arcs from `X` to `Z` are `s → b₁` and `a_k → b₁` only, so a Hamiltonian cycle `C` of `Y` enters `Z` once, at `b₁`, follows a Hamiltonian path `Q` of `B − s` from `b₁` to some `b_t`, and returns to `X`. Every `b_t`, `t ≠ 1`, beats `s` (`s` has the single out-arc `s → b₁` in `B`), so `Q + (b_t → s → b₁)` is a Hamiltonian cycle of `B`, whence `Q = b₁ … b_l` by uniqueness. The remainder of `C` is a Hamiltonian path `R` of `A` from some `x₀` with `b_l → x₀` to `x₁ ∈ {s, a_k}`. If `x₁ = s`, the last arc of `R` is `a_k → s` (sole entry) and `R − s` is a Hamiltonian path of `A − s` from some `a_i`, `i < k`, to `a_k`; since `s → a_i`, it closes to a Hamiltonian cycle of `A`, so `R − s = a₁ … a_k` and `C = a₁ … a_k s b₁ … b_l` (`b_l → a₁` is a cross arc). If `x₁ = a_k`, then `s ∈ R` has no predecessor on `R` (its only in-neighbour `a_k` is last), so `R` runs from `s` to `a_k` and closes with `a_k → s` to the cycle of `A`: `R = s a₁ … a_k`, `C = s a₁ … a_k b₁ … b_l` (`b_l → s` is an arc of `B`). Both are Hamiltonian cycles of `Y`. For the arc junction replace `s` by `s → s′` (the arcs from `X ∪ {s′}` to `Z` are `s′ → b₁`, `a_k → b₁`; `s′`'s only in-neighbour is `s`, `s`'s is `a_k`), and the same two cases give `a₁ … a_k s s′ b₁ … b_l` and `s s′ a₁ … a_k b₁ … b_l`. ∎

**(3) Lower bound.** The pair of Hamiltonian cycles of a gluing determines its three changed arcs, hence the junction, the two runs and the rooted pieces, and an isomorphism of gluings preserves in/out-degrees, hence maps entry piece to entry piece; so distinct rooted pairs give non-isomorphic gluings, and vertex and arc gluings are non-isomorphic. By THM-4868 (F), (C) the rooted classes `(A, s)` with `A ∈ UH(k + 1)` and `s` a sole entry number `F(2k − 3)`, and the classes of `UH(k + 2)` with two (consecutive) sole entries number `F(2k − 3)`; converse symmetry gives the same for `B`. Hence for `n ≥ 6`

`D_n ≥ E(n − 1) + E(n − 2)`, `E(s) = Σ_{k + l = s, k, l ≥ 2} F(2k − 3) F(2l − 3)` (ordered pairs) `= (n − 2) F(2n − 9) − L(2n − 10)`,

the last equality because `Σ_{k ≥ 2} F(2k − 3) xᵏ = x²(1 − x)/(1 − 3x + x²)`, so `Σ_n (E(n − 1) + E(n − 2)) xⁿ = x⁵(1 + x)(1 − x)²/(1 − 3x + x²)²`, the generating function of OEIS A377857 shifted by two (packet receipt `gf_identity`, `n < 60`). Equality holds at order `n` iff every two-cycle tournament of order `n` is a gluing; the census gives equality for all `n ≤ 16`, and the packet's structure script exhibits the decomposition (segment swap `P J S R / P S J R` with a one- or two-vertex junction `J`, uniquely Hamiltonian pieces, the single forward cross arc) on every class of orders 6–12. The classification implies HYP-9459 (1); its converse for all `n` is the open half of the count law.

**Also FINITE-EXACT (orders 6–14, every two-cycle class).** Exactly one or two sole exits and one or two sole entries; every non-critical vertex `v` has `hc(Y − v) ∈ {1, 2}` and exactly one or two of them have `hc(Y − v) = 1`; `3 ≤ |NC| ≤ n − 2`; and `#cut arcs + |NC| = n` (for `hc = 1` THM-4868 gives `n + 1`; for `hc = 3`, orders 6–13, the excess over `n` is 0 or 1). Through order 12, `NC` is a cyclic interval of exactly one of the two cycles and every cut arc lies on both.

**Related.** THM-4862 (b), (g) (digit-2 classes = two-cycle tournaments through the sole-exit bijection; the law as a conjecture from ten terms), THM-4865, THM-4867 / opus-S74 (extension lemma, Shelor), THM-4868 (kind-pasteur: roundness, sole exits, `F(2n − 5)`, `n + 1`), HYP-9459, A377857 (Florez 2024).
