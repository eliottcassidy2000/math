---
id: THM-4820
title: "Eisenstein merge debt and the catalyst threshold: merging unit charges 1, omega, omega^2 by a binary tree in which each merge costs the positive part of beta(x,y) = Q(x+y) - Q(x) - Q(y), every tree with root x costs floor(Q(x)/3) plus nonnegative integer slacks, so the catalytic minimum is g(x) >= floor(Q(x)/3), and neutral triples 1 + omega + omega^2 are free catalysts; an odd n prime to 3 with only inert prime factors and (n^2+2)/3 not of the form v^2 + 3t^2 has no tight split, so g(n) = (n^2+2)/3 for n = 11 (41), and the charges (11+k, k, k) cost 41 exactly when k >= 36 (the owner's threshold; derived concurrently and independently in Session75's EISENSTEIN_INVOICES note); g is floor(Q/3) or one more, exactly, on every point of norm <= 2000, with twenty exceptional norms, while g(101) = floor(101^2/3) + 2; thresholds k*(n) for n <= 15"
status: >
  PROVED (mac-mini-2026-10-27-slab; independently audited 2026-10-27, all numbers PASS, labels corrected:
  04-computation/experiments/regular_transitive_slab_20261027/audit_catalyst/REPORT_THM4820.md): the count formula for
  beta; neutral triples are free catalysts and f(n,0,0) = n(n-1); g(n) = min_k f(n+k,k,k); the potential bound with its
  tightness criterion; the tight-split criterion for rational integers and hence g(n) >= (n^2+2)/3 for n = 11, 23, 25,
  29, 41, 47; f(11+k,k,k) = 41 for every k >= 36 (from (C), (A2) and the verified certificate); the window-exactness
  lemma of (D); the lower bound k*(n) >= (n^2 - n - g(n))/3. FINITE-EXACT: f(46,35,35) = 42 by exact dynamic
  programming over every count vector, so the threshold is sharp (re-certified by the audit, 500,209,345 pairs); the row
  f(11+k,k,k), k <= 62; the exact values of g on every point of norm <= 2000 and on rational n <= 86; k*(n) for n <= 15.
  CONCURRENT: Session75's 05-knowledge/results/session75_final_mathematics_20261027/EISENSTEIN_INVOICES.md (committed
  minutes after this theorem's first push; not seen before filing) reaches the same reading, the same potential, the
  same 41 and k >= 36, plus an exceptional-set lifting lemma and the infinite family 11^(1+820k); the owner's statement
  evidently comes from that work. OPEN: HYP-9426 (3)-(5) ((3) as first drafted is REFUTED by g(101)).
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
  - 05-knowledge/results/session75_final_mathematics_20261027/EISENSTEIN_INVOICES.md
  - 04-computation/experiments/regular_transitive_slab_20261027/audit_catalyst/REPORT_THM4820.md
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
- β = 1 with r = s = 1, a lossless paid merge (two single units at 60°, or larger piles at other angles);
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

The value stays at 41 through k = 62 (the audit's box (76, 62, 62)).

## (D) The lattice table (FINITE-EXACT; exact values)

`lattice_debt.c` relaxes g over trees whose intermediate values have norm at most W. The windowed values are exact.

**Window lemma (PROVED; found by the audit).** Every subtree of a tree of cost C costs at most C. By (B), every node v
of the tree therefore satisfies ⌊Q(v)/3⌋ ≤ C, that is, Q(v) ≤ 3C + 1. So the windowed value g_W(x) equals g(x) whenever
W ≥ 3 g_W(x) + 1. At W = 2500 this covers every point of norm < 2500, and at W = 7500 every rational n ≤ 86. A
window-free computation of the zero-slack-buildable points, in the audit, agrees point by point.

For every point of norm at most 2000 (W = 2500):

- **Two values only.** g ∈ {⌊Q/3⌋, ⌊Q/3⌋ + 1}.
- **Where the +1 occurs.** It occurs only at Q ≡ 1 (mod 3), at 20 of the 328 such norms in (1, 2000]: 121, 229, 313,
  421, 484, 529, 625, 733, 841, 853, 1033, 1093, 1213, 1429, 1597, 1609, 1681, 1693, 1825, 1936. Every point at each of
  these norms is exceptional.
  - **No tight split at all.** The rational-integer norms 11², 23², 25², 29², 41² have no tight split, so (C) applies.
  - **By propagation.** The other fifteen have tight splits, but every tight split has an exceptional part. For example,
    229 splits as (−11) + (−6 − 12ω), and −11 is exceptional.
- **Rational integers.** For n ≤ 50, g(n) = ⌊n²/3⌋ + 1 exactly at n = 11, 22, 23, 25, 29, 41, 44, 46, 47, 50. The even
  ones are 2 or 4 times an odd exception.
- **To n = 86 (W = 7500).** The window adds the exceptions 51, 58, 59, 82, 83, 85.
  - **51** has norm ≡ 0 (mod 3), the only such exception below norm 7500. It has no tight split at all: for norm ≡ 0,
    tight means β ∈ {−2, …, 1}, and none of the norm equations 2601, 2599 = 23·113, 2603 = 19·137, 2605 = 5·521 has a
    nontrivial solution (audit).
  - **59** has a tight split and is exceptional by propagation.
  - **Where (C) applies.** (C) applies to 83 and 85, which are exceptional. It does not apply to 53, 55, 59 or 71; of
    these, computation shows 53, 55 and 71 are not exceptional and 59 is.
- **The first excess of 2 (audit).** g(101) = 3402 = ⌊101²/3⌋ + 2.
  - **Lower bound.** A tree of cost ≤ Φ + 1 would need a first split with β ∈ {−4, …, 2}, and every one of the
    corresponding norm equations for 101 fails.
  - **Upper bound.** An explicit tree of cost 3402, independently verified.

  On every norm below 10,201, g ∈ {⌊Q/3⌋, ⌊Q/3⌋ + 1} (W = 10,500). The same arithmetic gives an excess of at least 2 for
  78 odd n ≤ 3000, beginning 101, 115, 197, 235, 257, 289.
- **The sequence g(n), n = 1, …, 16:** 0, 1, 3, 5, 8, 12, 16, 21, 27, 33, 41, 48, 56, 65, 75, 85.

## (E) Catalyst thresholds (PROVED lower bound; FINITE-EXACT values)

A tree with L = n + 3k leaves and root n has Σβ = n² − L, and Σβ ≤ cost. So any tree attaining g(n) needs

    k ≥ (n² − n − g(n))/3.

The exact thresholds, from `merge_debt.c` with K = 4n + 8 (and K = 90 for n = 15), are:

| n | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| k*(n) | 0 | 1 | 2 | 4 | 7 | 10 | 14 | 19 | 24 | 31 | 36 | 44 | 52 | 61 | 70 |

For n = 11 the lower bound gives 23, against the true 36. Orders n ≥ 16 need K > 4n + 8.

The catalyst lowers the cost from n(n − 1) to g(n) ≥ ⌊n²/3⌋. So the saving is at most two thirds asymptotically. That it is asymptotically two thirds, and that k*(n) grows like a constant times n², both need an upper bound g(n) ≤ n²/3 + o(n²), which is CONJECTURE (HYP-9426 (5)). The data: g(n) − ⌊n²/3⌋ ∈ {0, 1} for n ≤ 100, and g(101) − ⌊101²/3⌋ = 2.

## Reading (typed)

- **IDENTITY (lattice algebra).**
  - 1, ω and ω² are the weights of the A₂ fundamental representation, and 1 + ω + ω² = 0 is the colour singlet. This is the external note §4 dictionary, in which arcs map to roots and the n = 3 roots are the six Eisenstein units.
  - The constant 3 is the index of √−3·Z[ω] in Z[ω].
  - The only lossless paid merges are the β = 1 merges of two norm-≡1 piles. Each such merge spends one unit of debt and gains one unit of Φ, that is, three units of norm.
- **IDENTITY (arithmetic of the (C) class).** For the exceptions of (C), exceptional status is an inertness statement in Z[ω]: 11 is inert, so it has no perpendicular split, and 41 is inert, so 11 has no obtuse split. In general it is not: 229 is a split prime, 59 has tight splits, and 51 = 3·17. Session75's lifting lemma ('h + 1_E is a subsolution iff every zero-slack merge with parent in E has a child in E') is the general mechanism, and its parity-shell criterion yields the infinite family 11^(1+820k).
- **ANALOGY (Session75 dyadic ledger).** Session75's sharp four-child witness has its own mod-3 threshold: B copies have zero global excess iff 3 | B. Its debt is the whole-level positive part, which is smaller than the per-merge positive part used here. The two thresholds are different quantities on the same pair interaction.
