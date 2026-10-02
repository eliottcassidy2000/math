---
id: HYP-9169
title: "The edge multiset dimension of the hypercube satisfies ln edim_m(Q_d) = Theta(d^(1/3)): sparse random landmark sets of size exp(O(d^(1/3))) resolve the edges of Q_d, matching the entropy lower bound exp((0.6215 - o(1)) d^(1/3))"
status: >
  OPEN. The lower half is PROVED (THM-4525, Lemma L5, entropy). The upper half is supported by:
  - the sparse union bound, in which the paper's random-flow forest lemma
    is used at any inclusion probability q: ln M_d / d^(1/3) lies in
    [2.76, 2.80] for d = 11..32 (VERIFIED in double precision with a
    safety margin, not interval-certified);
  - a heuristic count of central binomial atoms.
  Allikvere's Open Problem 2 (arXiv:2608.09983) asks for the growth rate.
  This hypothesis is the precise form suggested by the 2026-10-01 lane.
source: collatz-procgen-20260922 session, edim lane (2026-10-01); promoted with THM-4525
related:
  - 01-canon/theorems/THM-4525-edge-multiset-dimension-of-q6-is-15.md
  - 05-knowledge/results/procgen_edim_20261001_edge_multiset_dimension.md
---

# HYP-9169 — ln edim_m(Q_d) = Theta(d^(1/3))

**Lower half (PROVED, THM-4525 L5).** Choose an edge `e` uniformly at random.
- Its histogram `H_e` must have entropy `log2(d 2^(d-1))`.
- By subadditivity, that entropy is at most the sum of the level entropies.
- Each level count has mean `m C(d-1,r)/2^(d-1)`, so its entropy is at most the geometric-law value `g(mean)`.
- Only about `sqrt(2 d ln m)` levels have mean at least 1, and each contributes at most about `log2 m` bits.

Hence `ln m >= (c - o(1)) d^(1/3)`, with `c = (ln 2 / sqrt 2)^(2/3)`.

**Upper half (OPEN).** Take a random set `S` of size `m`. A typical pair of edges collides with probability about the product of `sqrt(2 d ln m)` central binomial atoms, each of order `m^(-1/2)`. A union bound over about `d^2 4^d` pairs then succeeds once `(ln m)^(3/2) >> sqrt d`, that is, once `ln m >> d^(1/3)`.

Making this rigorous needs a uniform, `q -> 0` version of the paper's forest-lemma estimates (its Sections 7–9). The lane's numerics stay flat at `ln M_d / d^(1/3) ~ 2.77` for `d = 11..32`.

**Explicit sets.** Small explicit sets (THM-4525: `Q_7 <= 19`, ..., `Q_12 <= 76`) are far below density `1/2`. This answers the qualitative half of Open Problem 4: sparse constructions are substantially better.

**Proof proposed (2026-10-02; independent audit owed, status unchanged until then).**
[`edge_multiset_dimension_growth_20261002.md`](../results/edge_multiset_dimension_growth_20261002.md) proves both halves:
0.8146 = kappa_- <= liminf ln edim_m(Q_d)/d^(1/3) <= limsup <= kappa_+ = (3 sqrt2 ln 2)^(2/3) = 2.0526. The upper half
uses sparse random sets with a uniform, q -> 0 version of the forest lemma (Theorem A); the lower half sharpens L5
(Theorem B). It also proves edim_m(Q_d) < 2 exp(5 d^(1/3)) for every d >= 6 and certifies the sparse table for
11 <= d <= 64 in interval arithmetic. The parallel lane note
[`procgen_edim2_20261001_growth_and_uniform_bounds.md`](../results/procgen_edim2_20261001_growth_and_uniform_bounds.md)
(audit also owed) proves the same two constants independently, with different forests, so the two derivations
can serve as cross-checks for the audit.
