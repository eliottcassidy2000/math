---
id: THM-856
title: THE HUNTER TREE BOUND CROSSES THE SEVEN-COMB WALL — second-order Bonferroni over a spanning tree is coercive at m′ = 7 remaining combs (and at no m′ ≥ 8): with per-comb density exactly 2/13 and pairwise overlaps at their equidistribution value 4/169, the criterion 4(m′−1) ≥ 13(2m′−13) holds iff m′ ≤ 7.5. The union-bound schema of THM-815 Part C is PROVABLY empty at m′ ≥ 7 (no per-comb bound |I ∩ D_x| ≤ αL + β(x) can be coercive: α < 2/13 is false for large x by equidistribution, α = 2/13 sums to 14L/13 ≥ L); the tree bound is its exact one-rung-stronger replacement. Global pairwise comb overlap has the exact Bezout form μ(D_{x₁} ∩ D_{x₂}) = (a+b)²/(169ab) + O(1/ab) ≥ 4/169 (a,b the coprime reduction) — minimized precisely at equal/consecutive speeds, which is exactly the observed failure locus of the tree bound
status: SCHEMA PROVED (Hunter/Kounias inequality + exact density arithmetic + the m′ ≤ 7.5 computation + the union-bound no-go) + VERIFIED on exact radius-7 pilot packets (prefix {1..5}: 4/4 generic packets coercive, Hunter ≥ +0.033; consecutive-speed packets fail the tree bound exactly as the Bezout minimum predicts while remaining uncovered ≥ 0.04 in fact). NOT a closure of the radius-7 chart: two named residuals below
source: opus-2026-07-15-S312 (owner: overnight marathon toward LRC(14) through formula generation; the frontier's item-3 ask "seek a replacement potential using overlap debt")
depends_on:
  - THM-815 Part C   # the recursion whose union bound dies at 7 combs
related: [THM-778 (mechanical words — the near-equal residual's tool), THM-855 F6 (the moment-closure lens that led here), LRC14-FRONTIER item 3]
verification: 05-knowledge/results/seven_comb_resonance_pilot_opus_S312.out, hunter_tree_wall_crossing_opus_S312.out
---

# THM-856 — the Hunter tree bound at the seven-comb wall

## 1. The no-go (why THM-815's schema is empty at m′ ≥ 7)

Any bound of the schema |I ∩ D_x| ≤ αL + β(x) (single comb vs single interval,
β(x) → 0) must have α ≥ 2/13: the comb D_x has density exactly 2/13, and
|I ∩ D_x|/L → 2L/13 as x → ∞ (equidistribution), so α < 2/13 is FALSE for
large x. At α = 2/13 and m′ = 7 remaining combs, the coverage constraint reads
14L/13 + Σβ ≥ L — satisfied identically: NO constraint on the lifts. The same
holds restricted to any fixed safe set E in place of I (the restricted density
also tends to 2/13 — see the periodicity lemma below). First-moment schemas
cannot cross the wall; this is a theorem, not an obstacle. ∎

## 2. The crossing (Hunter's inequality)

Hunter–Kounias: for any events A_i and ANY spanning tree T on the index set,
μ(∪A_i) ≤ Σμ(A_i) − Σ_{(i,j)∈T} μ(A_i ∩ A_j). Applied to A_i = D_{x_i} ∩ E:

> **uncovered ≥ μ(E) − Σ_i μ(D_i ∩ E) + max_T Σ_T μ(D_i ∩ D_j ∩ E).**

With densities 2/13 + O(1/x) and pairwise overlaps at the equidistribution
value (4/169)μ(E) + err, the criterion for coercivity at m′ combs is

> **4(m′−1) ≥ 13(2m′−13) ⟺ 22m′ ≤ 165 ⟺ m′ ≤ 7.5.**

At m′ = 7: tree sum 24/169·μE vs needed 13/169·μE — **an 85% margin**. At
m′ = 8: 28 < 39 — the tree bound has its own wall at eight. ∎ (schema)

## 3. The exact pairwise overlap formula (global)

For speeds x₁ = ga, x₂ = gb with (a,b) coprime:
**μ(D_{x₁} ∩ D_{x₂}) = (a+b)²/(169ab) + O(1/(ab)) ≥ 4/169**, with equality in
the leading term iff a = b (impossible for distinct speeds — so the infimum is
approached by near-equal ratios, e.g. consecutive integers). *Proof sketch:*
simultaneous teeth ⟺ |jb − ka| ≤ (a+b)/13; each admissible m = jb − ka
contributes a coincidence window of length max(0, (a+b)/(13ab) − |m|/(ab));
sum the triangle. gcd g cancels. AM–GM gives ≥ 4/169. ∎
Reading: **overlap is never scarce; it is minimized exactly at near-equal
speeds** — the tree bound's failure locus is a THIN, structured family.

## 4. The periodicity lemma (the finite table for E-restricted masses)

For a prefix safe set E with rational endpoints and a scale-one lift
x = r + 13h: **x·(|E ∩ D_x| − (2/13)μ(E)) is EXACTLY periodic in h** (period
dividing an explicit Λ(E, r); verified: prefix {1,...,5}, r = 6: period 60,
exact). The per-comb data of the infinite radius-7 chart is a finite table of
rationals. *Proof:* the anomaly depends only on the positions of E's endpoints
in the comb's tooth-coordinate, i.e. on x·(endpoint) mod 1, which is periodic
in h with period = denominator/gcd data. ∎

## 5. Pilot verification (prefix {1,2,3,4,5}, exact ℚ)

μ(E) = 7/15, 10 components. Four random radius-7 packets (lifts h ∈ [2,40]):
Hunter bound = +0.033, +0.034, +0.036, +0.038 — ALL COERCIVE (non-coverage
PROVED per packet by one inequality; actual uncovered ≈ 0.13–0.14, i.e. the
independence prediction (11/13)⁷·μE = 0.145 is accurate to 9%). Consecutive
packets {499..505}, {32..38}: Hunter −0.001, −0.017 (fails, exactly the
Bezout-minimum locus) while actual uncovered = 0.055, 0.040 — still far from
tight.

## 6. What remains for a radius-7 closure (named residuals)

1. **The deficit lemma** (the crux): |μ(D_i ∩ D_j ∩ E) − (4/169)μE| ≤
   C(E)/min(x_i, x_j) with explicit C — pilot-consistent (observed
   |err|·x ≤ 0.31). This bounds the smallest lift of any Hunter-refuting
   packet: x_min ≤ X₀(E) ⟹ adjoin it, drop to m′ = 6, where THM-815's
   recursion is coercive again — closing the chart down to finitely many
   small-lift states plus residual 2.
2. **The near-equal residual**: packets with several near-equal speeds evade
   the tree bound. For 7 consecutive speeds the values (N+i)t form an AP with
   step t bounded away from 0 on E — an AP-window/mechanical-word argument
   (THM-778 machinery) should prove uniform non-coverage; pilot shows
   uncovered ≥ 0.04 there.

**This theorem supplies the frontier's requested "replacement potential using
overlap debt": the potential is Hunter's tree functional, its coercivity
window is exactly one rung (m′ = 7), and its failure locus is exactly
characterized (near-equal speeds) with the next tool named per residual.**
