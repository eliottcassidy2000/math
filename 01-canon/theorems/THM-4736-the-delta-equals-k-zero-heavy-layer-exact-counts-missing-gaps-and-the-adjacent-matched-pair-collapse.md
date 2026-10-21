---
id: THM-4736
title: "The next zero-heavy boundary: the mass-deficit-k layer. For pair-free weighted cores Q[w] with weighted feedback f <= k and all weights <= k (zero-heavy), the classes of total weight M_k - k (M_k = 2k^2 + 3k) and, for nontrivial strong cores, S_k - k (S_k = 2k^2 + k) number exactly Z_{k,k} = (2k-1)!! (C(3k,k) - 1) and Z^s_{k,k} = I_k C(3k-2,k) - (2k-3) I_{k-1} (I_k = matchings on 2k ordered points crossing every proper cut: 1, 2, 10, 74, 706). The deficit identity forces a perfect matching of k backward arcs on 2k unit endpoints and leaves two families: (A) every gap present with gap deficits summing to k, (B) one gap missing and every other gap of weight k; (B) loses exactly the records whose deleted gap separated two endpoints matched to each other (they become a homogeneous pair), (2k-1)!! of them in general and (2k-3) I_{k-1} in the strong case, and all remaining records are pairwise non-isomorphic. This explains the packet's 42 at k = 2 (45 naive) and continues its exact layers (3)-(4) one step past the unit-gap boundary; the strong values are 11, 344, 15490 at k = 2, 3, 4"
status: >
  PROVED (pre-audit): the family reduction from the deficit identity, the collapse count, the distinctness of the
  surviving records (the packet's intrinsic-order argument extended to one missing gap and to one unit gap plus
  one reduced gap), the weighted feedback being exactly k (the sharp bound M_{k-1} < M_k - k forbids f <= k-1).
  FINITE-EXACT (zero_heavy_delta_k.out, zero_heavy_delta_k_strong.out, zero_heavy_delta_k_strong_k4.out):
  all records built as weighted tournaments, homogeneous pairs detected, weighted feedback computed by a subset DP
  (k <= 3), weight-preserving isomorphism by VF2 -- all cores: 42 (k = 2, equal to the packet's complete compiler),
  1245 (k = 3); strong: 11, 344, 15490 (k = 2, 3, 4), each equal to the displayed formulas, with exactly 3, 15 and
  1, 6, 50 homogeneous-pair collapses.
  CITED: charge_palettes_boundary_20261020g/zero_heavy.md (sharp bounds, deficit identity (6), exact layers (3)-(5)
  through delta = k-1, the k = 2 value 42), defect_boundary_20261020/sparse_theory.md (pair-free weighted cores,
  weighted feedback, the capped compiler), THM-4463.
  NOT PROVED: the layers delta >= k+1 (two missing gaps, shared endpoints, non-matching defect graphs enter);
  any unrestricted A000568 recurrence.
source: mac-mini-2026-10-20-debtcones (owner's prompt: the next zero-heavy boundary)
depends_on:
  - 05-knowledge/results/charge_palettes_boundary_20261020g/zero_heavy.md
  - 05-knowledge/results/defect_boundary_20261020/sparse_theory.md
related:
  - 05-knowledge/results/zero_heavy_ports_20261020g/README.md
  - 05-knowledge/hypotheses/HYP-9348-endpoint-layers-black-core-connectivity-and-the-next-zero-heavy-layers.md
scripts:
  - 04-computation/experiments/zero_heavy_20261020/zero_heavy_delta_k.py
  - 04-computation/experiments/zero_heavy_20261020/zero_heavy_delta_k_strong.py
outputs:
  - 04-computation/experiments/zero_heavy_20261020/zero_heavy_delta_k.out
  - 04-computation/experiments/zero_heavy_20261020/zero_heavy_delta_k_strong.out
  - 04-computation/experiments/zero_heavy_20261020/zero_heavy_delta_k_strong_k4.out
---

# THM-4736 -- one step past the unit-gap boundary

**PROVED (pre-audit) + FINITE-EXACT** as itemized in the status.

## 0. Setting (zero_heavy.md)

`Q[w]` is a pair-free weighted core (quotient by maximal homogeneous transitive chains, THM-4463)
with chain weights `w_i`; `f(Q[w]) = min_pi sum_{backward i -> j} w_i w_j` over total orders
(sparse_theory.md (3)). Fix `k >= 1`, `f <= k`, all `w_i <= k` (zero-heavy), `v = sum w_i`.
Choose an optimal order with `r` backward arcs (`r <= k`), endpoint count `t`, endpoint weight
`W`, non-endpoint vertex count `g <= t + 1` (`g <= t - 1` for nontrivial strong cores) and gap
deficit `U = sum_gaps (k - w_gap)` over the present gaps. The packet proved the sharp bounds
`v <= M_k = 2k^2 + 3k`, `v <= S_k = 2k^2 + k` (strong), the exact deficit identity

```text
M_k - v = (2k+1) q + (k + r - W) + k (2r - t) + k (t + 1 - g) + U,   q = k - r,           (6)
```

(all terms non-negative; `t + 1 - g` replaced by `t - 1 - g` in the strong case) and the exact
class counts `Z_{k,delta} = (2k-1)!! C(2k+delta, delta)` and `Z^s_{k,delta} = I_k C(2k+delta-2, delta)`
for `0 <= delta <= k-1`. At `delta = k` its compiler gave `42` for `k = 2` where (3) would give `45`.

## 1. The two families at deficit k

Set `M_k - v = k` in (6). `q >= 1` costs `2k + 1 > k`, so `r = k`: the defect graph has `k`
edges with `sum w_i w_j <= k`, hence every endpoint has weight `1` and `W = t`. If `t <= 2k-1`
then `k(2r - t) >= k` with `k + r - W = 2k - t >= 1` left over: impossible. So `t = 2k`, the
backward arcs form a **perfect matching on 2k unit endpoints**, `k + r - W = 0`, and

```text
k (t + 1 - g) + U = k:      (A)  g = 2k + 1 (every gap present), U = k;
                            (B)  g = 2k (exactly one gap missing), U = 0 (every other gap of weight k).
```

Family (A) = weak compositions of `k` into the `2k + 1` gap deficits with every part `<= k - 1`
(a part equal to `k` is a missing gap, i.e. family (B)): `C(3k, k) - (2k + 1)` per matching.
Family (B): `2k + 1` choices of the missing gap. Together `C(3k, k)` records per matching,
`(2k-1)!! C(3k,k)` in all -- the naive continuation of (3). In the strong case there are no
exterior gaps: `2k - 1` interior gaps, matchings crossing every proper cut (`I_k`), so
`C(3k-2, k) - (2k-1)` records in (A') and `2k - 1` in (B') per matching, `I_k C(3k-2,k)` in all.

## 2. The collapse

A record of family (B) whose missing gap lies between two endpoints `e_i, e_{i+1}` that are
**matched to each other** is not pair-free: every other vertex precedes or follows both, with all
its arcs forward, so `{e_i, e_{i+1}}` is a homogeneous pair (it is the clone chain of weight 2
oriented `e_{i+1} -> e_i`, i.e. the record is really a different, lighter-defect core). Every
other record is pair-free (Section 3). The number of (matching, adjacent matched pair) incidences
is `(2k-1) (2k-3)!! = (2k-1)!!` for all matchings on `2k` ordered points, and `(2k-3) I_{k-1}`
for crossing matchings: deleting an adjacent matched pair from a crossing matching leaves a
crossing matching on `2k - 2` points (the two cuts beside the pair are crossed by other pairs,
which survive), the pair can sit in any of the `2k - 3` interior positions of the smaller
sequence, never at an end (an end pair `{1,2}` would leave the cut `2|3` uncrossed), and the
correspondence is bijective. A missing *exterior* gap collapses nothing. Hence

```text
Z_{k,k}   = (2k-1)!! C(3k, k) - (2k-1)!!            = (2k-1)!! (C(3k,k) - 1),
Z^s_{k,k} = I_k C(3k-2, k) - (2k-3) I_{k-1}.
```

Values: `Z_{k,k} = 42, 1245, 51870, ...` (`k = 2, 3, 4`); `Z^s_{k,k} = 11, 344, 15490, 1520100, ...`
(with `I = 1, 2, 10, 74, 706`). The `k = 2` value is the packet's compiler count; the census
reproduces every number through `k = 3` (all) and `k = 4` (strong) with exactly the predicted
collapses `3, 15` and `1, 6, 50`.

## 3. Distinctness and exact feedback

*Weighted feedback is exactly k.* The displayed order costs `k`, and `f <= k - 1` is impossible
because the weight `M_k - k = 2k^2 + 2k` exceeds `M_{k-1} = 2k^2 - k - 1` (and `S_k - k = 2k^2`
exceeds `S_{k-1} = 2k^2 - 3k + 1`), the sharp bounds of the packet. So the layer is inside the
`f = k` stratum automatically (the census confirms: zero wrong-feedback records).

*Pair-freeness and non-isomorphism.* Vertices of weight `>= 2` are all gap vertices, carry
only forward arcs among themselves, and therefore have an intrinsic transitive order; a
weight-preserving isomorphism maps the `j`-th such vertex of one record to the `j`-th of the
other, so two records with different weight words on the large vertices are non-isomorphic.
Between consecutive large vertices (and before the first / after the last) lies a **unit
interval** of weight-one vertices. In family (B) and in family (A) without unit gaps every
unit interval is a single endpoint, except the merged interval of (B), which holds the two
endpoints beside the missing gap: they are distinguished by their matching partners (each a
singleton interval, hence fixed) unless matched to each other -- the collapsed case. In family
(A) a unit gap (deficit `k - 1`) is possible; for `k >= 3` at most one (two would cost
`2k - 2 > k`), accompanied by one gap of weight `k - 1 >= 2`; the packet's unit-gap analysis
(the two- or three-vertex interval, its roles fixed by partners outside or by the `C3`
rotation which leaves the record unchanged) applies verbatim. For `k = 2` two unit gaps can
occur; this case is the exhaustive `42` census. In every case the matching and the complete
weight word, including the position of the missing or reduced gaps, are recovered from the
class, so distinct records give distinct classes, and the records are pair-free (a homogeneous
pair would have to lie in one unit interval and be adjacent in the order: only the collapsed
case). Family (A) and family (B) records have different vertex counts (`4k + 1` vs `4k`) and
are never isomorphic. `#`

## 4. What changes at deficit k + 1

At `M_k - v = k + 1` the identity (6) admits `t = 2k - 1` with a shared endpoint
(`k(2r - t) = k`, `k + r - W = 1`), two missing gaps, a missing gap plus gap deficits, and
unit intervals holding up to three endpoints; the matching is no longer forced and the
collapse rule has several shapes. That layer is the next target (HYP-9348 (d)).

## 5. Files

`04-computation/experiments/zero_heavy_20261020/zero_heavy_delta_k.py` (`.out`: `k = 2, 3`,
all cores) and `zero_heavy_delta_k_strong.py` (`.out`: `k = 2, 3`; `_k4.out`: `k = 4`, feedback
DP skipped for speed, distinctness and collapses checked). Records are built as weighted
tournaments (endpoints weight 1, gap vertices weighted), homogeneous pairs detected directly,
isomorphism by `networkx` VF2 with node-weight matching after an invariant bucket.
