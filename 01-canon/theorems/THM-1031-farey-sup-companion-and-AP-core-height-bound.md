---
id: THM-1031
title: The Farey sup-companion to THM-826, and the AP-core height bound — for the interval core {1,…,k} the WIDEST safe arc is δ({1,…,k};λ) = (1 − λ(k+1))/k, attained at the Farey pair (0/1, 1/k); this is the L^∞ analogue of kind-pasteur's L^1 profile theorem (THM-826 sums the same gap-remnants that this maximizes). CONSEQUENCE via THM-1001: a tight n-set whose max-peel core is the AP {1,…,n−1} satisfies max(A) ≤ 2(n−1). Also recorded: the covering lemma (a tight n-set contains a multiple of every q ≤ n), which is the search prune for the height census.
status: PROVED (δ formula: THM-826's two lemmas + the Farey adjacency facts i+j>k and the simultaneous minimization at (1,k); verified exactly n=4..14. Covering lemma: three lines. AP-core height bound: immediate from THM-1001.) The GENERAL height bound max(A) ≤ 3n remains a CONJECTURE (HYP-7390) — this proves only the AP-core case.
source: mac-mini-2026-07-18-S114 (owner: prove the height bound; extend other agents' ideas)
renumber_note: first pushed as THM-1030 at 10:32:10, 62s BEFORE opus's unrelated THM-1030 (killer-collapse-is-the-sieve) at 10:33:12. By the first-pusher protocol the number was mine, but opus's file was already wired into HYP-7420, their session log and close-out letter, whereas this one was newly created and referenced only by its own commit — so I ceded and renumbered here. Precedence recorded only so the ledger is unambiguous; no dispute.
depends_on:
  - THM-826   # kind-pasteur's Farey profile theorem (the L^1 statement this is the L^infty companion to)
  - THM-1001  # the safe-interval element bound w <= 2L/delta(A\{w})
external: LRC(n) SETTLED.
related:
  - HYP-7390  # the height-bound conjecture max(A) <= 3n
  - THM-819   # the primitive harmonic law (THM-826's first segment)
  - THM-763   # the exponential unconditional height bound this aims to replace
  - HYP-2930..2934  # the Farey mediant/operator/bridge threads
---

# THM-1031 — the Farey sup-companion, and the AP-core height bound

**One line.** kind-pasteur's THM-826 computes how *much* of the circle stays lonely;
this computes how *wide* its biggest piece is — the same Farey gap-remnants, maximized
instead of summed — and that sup is exactly what the height bound needs.

## (A) The covering lemma (PROVED)

> A tight `n`-set `A` (`M(A) = 1/(n+1)`) contains a multiple of **every** `q ≤ n`.

*Proof.* Fix `q ≤ n` and `a` with `gcd(a,q)=1`. Since `M(A) = 1/(n+1)`, at `t = a/q`
some `v ∈ A` has `‖va/q‖ ≤ 1/(n+1)`, i.e. `|va|_q ≤ q/(n+1) < 1` because `q ≤ n`.
So `|va|_q = 0`, i.e. `q | va`, and `gcd(a,q)=1` gives `q | v`. ∎

(This is the corpus's *covering* condition, here derived directly from tightness; it is
the prune that makes the height census feasible.)

## (B) The sup-companion of the Farey profile (PROVED)

Setup as in THM-826: `G(λ) = {t : ‖jt‖ ≥ λ, j = 1..k}`, and by its Lemma 1 (nesting)
and Lemma 2 (no intrusion), for `λ < 1/(k+1)` the good set is the **disjoint** union of
gap remnants, one per consecutive Farey pair `a/i < b/j` of `F_k`, of length
`(1 − λ(i+j))/(ij)`. THM-826 sums them. Maximizing them instead:

> **`δ({1,…,k}; λ) := ` widest safe arc ` = (1 − λ(k+1))/k`,**
> attained at the Farey pair `(0/1, 1/k)`.

*Proof.* Each remnant is `(1 − λ(i+j))/(ij)`, decreasing in both `ij` and `i+j`. For
consecutive Farey fractions of order `k` one has `i + j > k`, so `i+j ≥ k+1`; and the
pair `(i,j) = (1,k)` attains `i+j = k+1` **and** minimizes `ij = k` (any other adjacent
pair has `i,j ≥ 2`, hence `ij ≥ 2(k−1) > k` for `k > 2`). Since `(1,k)` minimizes both
arguments simultaneously, it maximizes the remnant. ∎

*Verified exactly* for `λ = 1/(n+1)`, `k = n−1`, `n = 4..14`:
`δ = 1/((n−1)(n+1))` in every case.

## (C) The AP-core height bound (PROVED)

Let `A` be tight with `L = 1/(n+1)` and suppose the max-peel core is the AP,
`A\{max A} = {1,…,n−1}`. THM-1001 gives `max(A) ≤ 2L/δ(A\{max A})`, and (B) with
`k = n−1`, `λ = L` gives `δ = 1/((n−1)(n+1))`. Hence

> **`max(A) ≤ 2·\frac{1}{n+1}·(n−1)(n+1) = 2(n−1)`.**

Combined with the S108 finite check (`{1,…,n−1, w}` is tight **iff** `w = n`), this pins
the AP-core case completely: `A = {1,…,n}`.

**Sharpness note.** The bound is *nearly* attained and correctly excludes near-misses:
the `n=4` sporadic `{1,3,4,7}` has `max = 7 > 6 = 2(n−1)`, and indeed its core `{1,3,4}`
is **not** the AP — consistent, and a check that the hypothesis is doing real work.

## (D) Why this does not yet give the general height bound

For a core `C` containing 1, the near-origin remnant gives
`δ(C) ≥ (1−L)/max(C) − L`, which is positive **iff** `max(C) < n` — i.e. only for the
AP-like (compact) cores. Once the core is spread (`max(C) ≥ n`), speed 1's danger arc of
radius `L` swallows the origin gap and `δ` must come from an interior Farey pair, whose
location depends on `C`'s actual denominators. That is exactly the residual: the general
`max(A) ≤ 3n` conjecture (HYP-7390) needs a lower bound on the *interior* gap structure
of arbitrary cores, not just the interval core THM-826/this theorem resolve.

*Artifacts:* `04-computation/lrc_height_threshold_census_macmini_S114.py` (+out).
Credits: **kind-pasteur THM-826** (the Farey profile whose gap-remnant atoms this
maximizes — the L^1/L^∞ pair), THM-819, THM-1001, and the HYP-2930..2934 Farey threads.
