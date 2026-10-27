---
id: THM-4820
title: "Eisenstein merge debt and the catalyst threshold: merging unit charges 1, omega, omega^2 by a binary tree in which each merge costs the positive part of beta(x,y) = Q(x+y) - Q(x) - Q(y), every tree with root x costs floor(Q(x)/3) plus nonnegative integer slacks, so the catalytic minimum is g(x) >= floor(Q(x)/3), and neutral triples 1 + omega + omega^2 are free catalysts; an odd n prime to 3 with only inert prime factors and (n^2+2)/3 not of the form v^2 + 3t^2 has no tight split, so g(n) = (n^2+2)/3 for n = 11 (41), and the charges (11+k, k, k) cost 41 exactly when k >= 36 (the owner's threshold, reproduced); g is floor(Q/3) or one more on every point of norm <= 2000, with twenty exceptional norms; thresholds k*(n) for n <= 14"
status: >
  PROVED (mac-mini-2026-10-27-slab; independent audit pending at filing): the count formula for beta; neutral triples
  are free catalysts and f(n,0,0) = n(n-1); g(n) = min_k f(n+k,k,k); the potential bound with its tightness criterion;
  the tight-split criterion for rational integers and hence g(n) >= (n^2+2)/3 for n = 11, 23, 25, 29, 41, 47; the lower
  bound k*(n) >= (n^2 - n - g(n))/3. FINITE-EXACT: a cost-41 tree on (47,36,36), checked by an independent
  value-arithmetic verifier, so f(11+k,k,k) = 41 for all k >= 36; f(46,35,35) = 42 by exact dynamic programming over
  every count vector, so the threshold is sharp; the table f(11+k,k,k), k <= 60; the lattice table of g on norms
  <= 2000 (window-restricted trees, norms <= 2500) and on rational n <= 86 (window 7500); k*(n) for n <= 14. RECONSTRUCTION: the owner's statement has no
  source in the repository at filing (searched: every text file, the dyadic, paid-ROOT, renewal and necklace branches).
  This theorem fixes the reading (unit Eisenstein charges, free binary merges, per-merge positive part of the
  Session75 dyadic pair interaction), under which the stated numbers hold exactly. OPEN: HYP-9426 (3)-(5).
session: mac-mini-2026-10-27-slab
source: 05-knowledge/results/regular_slab_prescribed_black_core_eisenstein_catalyst_20261027.md
scripts:
  - 04-computation/experiments/regular_transitive_slab_20261027/catalyst/merge_debt.c (+ merge_debt_11.out, tree_47_36_36.txt, thresholds.out)
  - 04-computation/experiments/regular_transitive_slab_20261027/catalyst/verify_tree.py (+ verify_tree.out)
  - 04-computation/experiments/regular_transitive_slab_20261027/catalyst/lattice_debt.c (+ lattice_W2500.txt, lattice_W7500.txt, rational_criterion.out)
  - 04-computation/experiments/regular_transitive_slab_20261027/catalyst/tight_splits.py (+ tight_splits.out)
related:
  - 05-knowledge/results/session75_overnight_black_root_spectral_20261027/session75-dyadic/README.md
  - 05-knowledge/results/paid_root_cycle_budget_20261026/README.md
  - 05-knowledge/results/external_preprints_mub6_local_smoothing_eth_merge_20261027.md
  - 01-canon/theorems/THM-4813-collatz-numerator-coordinate-forcing-cocycle-gauss-sum-intertwiner-and-sharp-spectral-exception-budgets.md
  - 01-canon/theorems/THM-4285-universal-two-eisenstein-norm-shell-and-exact-contact-trichotomy.md
  - 05-knowledge/hypotheses/HYP-9426-slab-exactness-black-endpoint-counts-eisenstein-exceptional-set-and-catalyst-thresholds.md
---

# THM-4820 — Eisenstein merge debt: the floor Q/3, free neutral catalysts, and 41 at k ≥ 36

**Setting.** Let ω² + ω + 1 = 0, and let Q(x + yω) = x² − xy + y² be the norm of Z[ω], the A₂ root lattice.
- **Charges.** A charge vector (a, b, c) is a multiset of a unit charges 1, b charges ω and c charges ω². Its value is a + bω + cω² = (a − c) + (b − c)ω.
- **Merge trees.** A merge tree is a binary tree whose leaves are the unit charges. Each internal node carries the sum of its children's values.
- **Cost.** A merge of y and z costs β(y, z)₊, where β(y, z) = Q(y + z) − Q(y) − Q(z). This is the pair interaction of Session75's dyadic ledger §1. The positive part per merge is the middle quantity of that ledger's chain D⁺_native ≤ Σ max(β(children), 0) ≤ Σ over pairs.
- **The two minima.** f(a, b, c) is the minimum cost of a tree on that multiset. g(x) is the minimum cost of a tree with root value x whose leaves may be any of the six units.

## (A) Elementary facts (PROVED)

**(A1) Count form.** For count vectors u and u′ with |u| and |u′| charges, β = 3 u·u′ − |u||u′|. This holds because ⟨ωⁱ, ωʲ⟩ = 1 if i = j and −1/2 otherwise, with β = 2⟨y, z⟩.

**(A2) Neutral triples are free catalysts.** f(u + (1,1,1)) ≤ f(u). The merges 1 + ω (β = −1) and then (1 + ω) + ω² = 0 (β = −2) are free, and merging the resulting 0 into the root costs β = 0.

**(A3) Without catalysts.** f(n, 0, 0) = n(n − 1). Every pile is a positive integer, and β(a, b) = 2ab > 0 telescopes to Q(n) − n.

**(A4) The catalytic minimum is a limit.** g(n) = min_k f(n + k, k, k), and f(n + k, k, k) is non-increasing in k by (A2).
- **Negative units are free.** A leaf −1, −ω or −ω² is the free merge ω + ω², 1 + ω² or 1 + ω respectively, each with β = −1.
- **Root value n forces the count shape.** A tree with root value n has a count vector of the form (n + k, k, k).

So f(n + k, k, k) = g(n) for all k ≥ some threshold k*(n).

## (B) The potential bound (PROVED)

Norms of Z[ω] are ≡ 0 or 1 (mod 3), since x² − xy + y² ≡ (x + y)² (mod 3). Put Φ(x) = ⌊Q(x)/3⌋. Write Q(y) = 3a + r and Q(z) = 3b + s with r, s ∈ {0, 1}. Then

    Φ(y+z) − Φ(y) − Φ(z) = ⌊(r + s + β)/3⌋ ≤ β₊

- if β ≥ 1, then (2 + β)/3 ≤ β;
- if β ≤ 0, then the floor is at most 0.

Units and 0 have Φ = 0. Telescoping over a tree, each internal value appears once as a parent and once as a child:

    cost(tree) = Φ(root) + Σ_merges slack,   slack = β₊ − ⌊(r + s + β)/3⌋ ∈ {0, 1, 2, …}.

Hence g(x) ≥ ⌊Q(x)/3⌋ and f(n + k, k, k) ≥ ⌊n²/3⌋ for every k.

**Tight merges.** A merge has slack 0 exactly in two cases:
- β = 1 with r = s = 1, a lossless paid merge at 60°;
- −(r + s) ≤ β ≤ 0.

If Q(x) ≡ 1 (mod 3), a tight split x = y + z with y, z ≠ 0 is either **perpendicular** (β = 0) or **obtuse with β = −1**. The other cases are impossible: β = 1 forces Q(x) ≡ 0, and β = −2 forces r = s = 0, which allows no negative β.

## (C) Rational integers without a tight split (PROVED)

Let n be odd with 3 ∤ n. Write y = a + bω and put u = 2a − b, so that Q(y) = (u² + 3b²)/4 and ⟨y, n⟩ = nu/2. Then:

    y ⊥ (n − y)         iff  (u − n)² + 3b² = n²,
    β(y, n − y) = −1    iff  (u − n)² + 3b² = n² + 2.

The parity condition u ≡ b (mod 2) holds automatically for odd n.

**When the two splits exist.**
- **Perpendicular split.** A nontrivial perpendicular split (b ≠ 0) exists iff n has a prime factor ≡ 1 (mod 3). If every prime factor of n is ≡ 2 (mod 3), the elements of norm n² are the associates of n, and none of them gives b ≠ 0 with u integral.
- **Obtuse split.** In n² + 2 = s² + 3b² we must have 3 | s, so the obtuse split exists iff M := (n² + 2)/3 = b² + 3t². Since M ≡ 1 (mod 8), this fails iff some prime ≡ 2 (mod 3) divides M to an odd power.

**Conclusion.** If n has only prime factors ≡ 2 (mod 3) and M is not of the form b² + 3t², then n has no tight split. In a tree of cost ⌊n²/3⌋ every merge is tight. Follow merges with a zero child, which are tight and keep the value. The first merge reached that has two nonzero children would be a tight split of n, and there is none. Hence

    g(n) >= ⌊n²/3⌋ + 1 = (n² + 2)/3 = M.

The hypotheses hold for n = 11, 23, 25, 29, 41 and 47:

| n | M = (n²+2)/3 | factorization of M |
|---|---|---|
| 11 | 41 | prime |
| 23 | 177 | 3·59 |
| 25 | 209 | 11·19 |
| 29 | 281 | prime |
| 41 | 561 | 3·11·17 |
| 47 | 737 | 11·67 |

In every row an inert prime divides M to an odd power.

**The catalyst threshold.** For n = 11 the cost 41 *is* the obstruction: a tree of cost 40 would need 41 = b² + 3t², and 41 ≡ 2 (mod 3) is inert.
- **The certificate.** `tree_47_36_36.txt` is a merge tree of cost exactly 41 on (47, 36, 36), with 118 merges and 119 unit leaves.
  - It was found by `merge_debt.c`.
  - It was re-verified by `verify_tree.py`. The verifier recomputes every β from Eisenstein values, not from the count formula, and checks the leaf multiset and the root value 11.
- **All k ≥ 36.** With (A2) and the lower bound, f(11 + k, k, k) = 41 for every k ≥ 36.
- **Sharpness.** The exact dynamic program f(u) = min over u = u₁ + u₂ of β₊ + f(u₁) + f(u₂), over every count vector ≤ (46, 35, 35), gives f(46, 35, 35) = 42.

So the minimum merge cost for charges (11 + k, k, k) is 41 exactly when k ≥ 36. The full row, f(11 + k, k, k) for k = 0, …, 36, is:

    110 109 106 103 100 97 94 92 89 87 84 82 80 78 75 73 71 69 67 65 63 61 59 57 55 53 52 50 49 48 47 46 45 44 43 42 41

The value stays at 41 through k = 60.

## (D) The lattice table (FINITE-EXACT inside the stated window)

`lattice_debt.c` relaxes g over trees whose intermediate values have norm at most W. The runs at W = 400, 900 and 2500 agree on every n ≤ 16. For every point of norm at most 2000 (run at W = 2500):

- **Two values only.** g_W ∈ {⌊Q/3⌋, ⌊Q/3⌋ + 1}. Trees realize these values, so g ≤ ⌊Q/3⌋ + 1 is PROVED on this range, and g = ⌊Q/3⌋ is exact at every non-exceptional point.
- **Where the +1 occurs.** It occurs only at Q ≡ 1 (mod 3), at 20 of the 328 such norms: 121, 229, 313, 421, 484, 529, 625, 733, 841, 853, 1033, 1093, 1213, 1429, 1597, 1609, 1681, 1693, 1825, 1936.
  - **Window-free exceptions.** The rational-integer norms 11², 23², 25², 29², 41² have no tight split at all, so (C) makes them exceptional for every tree.
  - **Window-dependent exceptions.** The other fifteen have tight splits, but only through an exceptional part. For example, 229 splits as (−11) + (−6 − 12ω), and the part −11 is exceptional. Their "+1" is therefore window-dependent.
- **Rational integers.** For n ≤ 50, g(n) = ⌊n²/3⌋ + 1 exactly at n = 11, 22, 23, 25, 29, 41, 44, 46, 47, 50. The even ones are 2 or 4 times an odd exception.
- **The larger window (W = 7500, n ≤ 86).** It agrees with W = 2500 on every n ≤ 50 and adds the exceptions 51, 58, 59, 82, 83, 85.
  - **Two new kinds.** 51 has norm ≡ 0 (mod 3), the first such exception. 59 has a tight split, so it is exceptional only by propagation.
  - **The criterion of (C).** It predicts 83 and 85 and excludes 53, 55, 71. All five predictions hold.
  - The sequence g(n) for n = 1, …, 16 is 0, 1, 3, 5, 8, 12, 16, 21, 27, 33, 41, 48, 56, 65, 75, 85.

## (E) Catalyst thresholds (PROVED lower bound; FINITE-EXACT values)

A tree with L = n + 3k leaves and root n has Σβ = n² − L, and Σβ ≤ cost. So any tree attaining g(n) needs

    k ≥ (n² − n − g(n))/3,   which is ≈ (2/9) n².

The exact thresholds, from `merge_debt.c` with K = 4n + 8, are:

| n | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| k*(n) | 0 | 1 | 2 | 4 | 7 | 10 | 14 | 19 | 24 | 31 | 36 | 44 | 52 | 61 |

For n = 11 the lower bound gives 23, against the true 36. Orders n ≥ 15 need K > 4n + 8.

The catalyst lowers the cost from n(n − 1) to about n²/3. Asymptotically it saves two thirds of the cost, and it needs on the order of n² neutral triples to do so.

## Reading (typed)

- **IDENTITY (lattice algebra).**
  - 1, ω and ω² are the weights of the A₂ fundamental representation, and 1 + ω + ω² = 0 is the colour singlet. This is the external note §4 dictionary, in which arcs map to roots and the n = 3 roots are the six Eisenstein units.
  - The constant 3 is the index of √−3·Z[ω] in Z[ω].
  - The only lossless paid merges are the β = 1 merges of two norm-≡1 piles. Each such merge spends one unit of debt and gains one unit of Φ, that is, three units of norm.
- **IDENTITY (arithmetic of the exceptions).** Exceptional status is an inertness statement in Z[ω]. 11 is inert, so it has no perpendicular split, and 41 is inert, so 11 has no obtuse split.
- **ANALOGY (Session75 dyadic ledger).** Session75's sharp four-child witness has its own mod-3 threshold: B copies have zero global excess iff 3 | B. Its debt is the whole-level positive part, which is smaller than the per-merge positive part used here. The two thresholds are different quantities on the same pair interaction.
