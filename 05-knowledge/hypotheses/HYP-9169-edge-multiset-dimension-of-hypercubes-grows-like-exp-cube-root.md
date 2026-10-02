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

**UPDATE 2026-10-01 (THM-4534, edim2 lane; orchestrator-audited): the Theta form is PROVED.**
`0.8146 <= liminf ln edim_m(Q_d)/d^(1/3) <= limsup <= 2.0526`.
- **Upper bound.** Sparse random landmarks with an analytic union bound: Fourier atom bound plus star forests. The constant
  `C* = (3 sqrt2 ln2)^(2/3)` is exactly the reach of the forest-lemma method.
- **Lower bound.** A sharper evaluation of the entropy inequality gives `c* = (3 ln2/(2 sqrt2))^(2/3)`.
- **What remains OPEN.** Whether the limit exists, and its value. The conjecture is now the sharper statement that
  `lim ln edim_m(Q_d)/d^(1/3)` exists, possibly equal to `C*`, the heuristic threshold of uniformly random landmark sets.

**Second, independent derivation (2026-10-02, branch `claude/project-thread-96889h`; blind independent audit found
no false claim, wording fixes applied).**
[`edge_multiset_dimension_growth_20261002.md`](../results/edge_multiset_dimension_growth_20261002.md) proves the same
two constants with different forests (Lemma H and the toward-the-centre forest for the upper half, Theorem A; the
same L5 inequality evaluated more sharply for the lower half, Theorem B). It adds the lower half of the Open Problem 4
dichotomy (Theorem D(2)), the explicit bound edim_m(Q_d) < 2 exp(5 d^(1/3)) for every d >= 6 (computer-assisted at
its finite inputs), and a sparse table certified in interval arithmetic for every 11 <= d <= 64. The companion note
[`edge_multiset_dimension_op3_20261002.md`](../results/edge_multiset_dimension_op3_20261002.md) settles Open Problem 3
with one estimate at density 1/2 for every d >= 10 (audit owed).
