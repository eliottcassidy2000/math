---
id: THM-4479
title: "Collatz's distance to provability in the strategy cube tends to zero at the sharp rate 2^(-(1-h)k): flipping exactly the undecided residues is provable, and expanding necklaces force at least 2^(hk)/(3k^2) flips (HYP-9138 proved)"
status: >
  PROVED + INDEPENDENTLY AUDITED + FINITE-EXACT.
  In the strategy cube of THM-4474 (a sign on each odd residue mod 2^k,
  T_sigma(n) = n/2 or (3n + sigma(n mod 2^k))/2), let delta_k be the least
  number of odd residues at which a class-(i) strategy differs from
  Collatz. Class (i) means bounded-lookahead provability, equivalently
  every cycle of the parity graph has odd density < log_3 2.
  (1) Flipping exactly Bad_k, the residues on which Collatz has no k-step
  descent, gives a class-(i) strategy sigma_k. Its critical density is
  F_k, the best lower rational approximation of log_3 2 with denominator
  <= k, attained by an upper Christoffel word. Every n > n_0(k) falls
  below itself within k steps. Hence delta_k <= |Bad_k| <= 2^(hk), with
  h = h(log_3 2) = 0.9499555.
  (2) delta_k >= N_k >= 2^(hk)/(3k^2), where N_k counts the binary
  necklaces of length k with more than k log_3 2 ones.
  So (1/k) log_2(delta_k/2^(k-1)) -> -(1-h) = -0.0500445, and HYP-9138 is
  SETTLED with a sharp exponent.
  FINITE-EXACT: delta_k = 1, 2, 2, 4, 5, 9, 14, 23 for k = 2..9;
  40 <= delta_10 <= 44; 52 <= delta_11 <= 72. The approximants sigma_k are
  transitive (every positive orbit reaches 1) for k <= 300. DRIFT: the
  5n+-1 provable class is empty at levels 2..6 and nonempty from 7, and
  5n+1's exact distance at k = 7 is 29/64; whether it tends to 0 is OPEN.
  Collatz itself stays in class (iv) at every level. Nothing here bears on
  its truth.
  UPDATE 2026-09-26 (THM-4485): new lower bounds delta_11 >= FVS^odd_11 = 58
  (was 52) and delta_12 >= 95 (was N_12 = 70), from exact feedback sets of
  the expanding cycles; the upper bounds remain 72 (k = 11) and 131 (k = 12,
  pruned).
  UPDATE 2026-09-26 (THM-4482): the upper Christoffel word is one maximizer
  of sigma_k, not the only one. The maximal cycles are all periodic
  concatenations of density-F_k first-descent blocks, a positive-entropy
  family for k >= 4. Other provable level-k strategies can exceed F_k
  (e.g. 3/5 at level 4), up to the numerator bound G_(2^(k-2)).
source: collatz-procgen-20260922 session, cube-distance lane (2026-09-25), proving the session's HYP-9138; audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md
related:
  - 05-knowledge/hypotheses/HYP-9138-cube-distance-to-provability.md (now SETTLED)
  - 01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md (pairing family, upper bound)
  - 01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md (pairing family and arbitrary edits, lower bound by a different mechanism)
  - 01-canon/theorems/THM-4477-cauchy-schwarz-price-bound-and-am-fair-criticality.md (distribution-only barrier 2(1-h))
note: 05-knowledge/results/procgen_cubedist_20260925_distance_to_provability.md
scripts: 04-computation/experiments/procgen_cubedist_20260925_{lib,bad,exact,prune,controls,run}.py and procgen_cubedist_20260925_engine.c (sha256 in the note, section 8; re-checked unchanged)
script_audit: 04-computation/experiments/procgen_cubedist_20260925_orchestrator_check.py
output: 05-knowledge/results/procgen_cubedist_20260925.out
output_audit: 05-knowledge/results/procgen_cubedist_20260925_orchestrator_check.out
certificate: 05-knowledge/results/procgen_cubedist_20260925_k10_certificate.json.gz
output_sha256: 4596290838f30a09589b1c1d0072004793bec4ff10f9c8a13494ff74ac2a2dd5
certificate_sha256: b351d6dccf135c585f9a3fdd00a5c8fdb7a4e7082e868e08094a93de3e3ce198
hash_basis: raw bytes
audit: >
  The orchestrator checked Lemma 1.1 (heredity), Lemma 1.2 (blocks),
  Theorem 1 (a)-(d) and Theorem 2 (necklace packing, binomial and entropy
  estimate) line by line. Independent code
  (procgen_cubedist_20260925_orchestrator_check.py, written from the
  note's statements without reading the lane's scripts) confirms:
  * |Bad_k| by direct simulation equals a ballot DP for k <= 20. It equals
    the note's table, and Bad_k lies in 3 mod 4 with -1, -5 in Bad_k.
    The DP agrees at k = 40, 100, 200, and |Bad_k| <= 2^(hk) for k <= 200.
  * sigma_k has an integer potential certificate at threshold F_k for
    k <= 16, and exact Karp gives rho_max = F_k for k <= 11. For
    k <= 16, F_k is attained by the periodic orbit of c/(2^d - 3^a), e.g.
    319/13 for 5/8, checked to be a cycle of G_(sigma_k).
  * n_0(k) by first-descent-word DFS is 1, 4 and 24 at k >= 2, 5 and 8
    (k <= 20). No n in [2, 2*10^5] fails to descend within k steps, and
    every n <= n_0(k) reaches 1.
  * N_k by enumeration equals N_k by Burnside and the note's table.
    N_k >= 2^(hk)/(3k^2) holds for k <= 200. Rotation of parity words is
    an edge of Collatz's parity graph for k <= 12.
  * delta_k = 1, 2, 2, 4, 5 for k = 2..6, by exhaustive search over flip
    sets of increasing size (a third method after the cube lane's MaxSAT
    and this lane's HiGHS implicit hitting set).
  * delta_7 = 9 and delta_8 = 14, re-derived by an implicit hitting set
    with OR-tools CP-SAT: a third solver after MaxSAT and HiGHS.
  * 5n+-1: no class-(i) strategy at levels 2..5 (exhaustive, 65,812
    strategies). The level-6 hitting-set problem is infeasible (148 no-goods).
    A class-(i) strategy exists at level 7 with rho_max = 3/7, and 5n+1's
    exact distance there is 29/64.
  * The k = 10 certificate: all 1949 stored no-goods are re-derived as
    expanding closed walks, and the 44-flip upper set is class (i) with
    rho_max = 5/8. The script generates its own 237 seeds. HiGHS on this
    independent model (2186 no-goods) proves the minimum hitting set is
    exactly 40, so delta_10 >= 40. CP-SAT with 2 workers returned UNKNOWN
    after 3000 s on the same model.
  The lane's full pipeline was re-run (787 s, 358 MB). Its output is
  identical to the committed .out except for timing fields.
---

# THM-4479 -- the strategy cube's distance to provability, with a sharp exponent

**PROVED + INDEPENDENTLY AUDITED.** Full note, with the exact optima, the pruning experiments, the 3n−1 control and the 5n±1 DRIFT control: [procgen_cubedist_20260925_distance_to_provability](../../05-knowledge/results/procgen_cubedist_20260925_distance_to_provability.md).

## 1. Setting (THM-4474)

* **Strategies.** A level-`k` sign strategy `sigma` puts `+1` or `-1` on each odd residue mod `2^k`. Then `T_sigma(n) = n/2` for even `n` and `(3n + sigma(n mod 2^k))/2` for odd `n`. Collatz is `sigma = +`.
* **Parity graph.** `G_sigma` has nodes `Z/2^k`. From `s` there are two edges, to the two lifts mod `2^k` of `T_sigma(s) mod 2^(k-1)`.
* **Class (i).** Bounded-lookahead provability. By Theorem A of THM-4474 it holds iff every cycle of `G_sigma` has odd density `< c = log_3 2`; a cycle with `a` odd nodes and length `p` is expanding iff `3^a > 2^p`.
* **The distance.** `delta_k` is the least number of odd residues at which a class-(i) level-`k` strategy differs from Collatz. Its Haar fraction is `delta_k / 2^(k-1)`.
* **The undecided residues.** Let `a_j(x)` be the number of odd terms among `x, Tx, ..., T^(j-1)x` (Collatz `T`) and `M_j(x) = 3^(a_j(x))/2^j`. The first `k` parities depend only on `x mod 2^k` (Terras), and `Bad_k = {r mod 2^k : M_j(r) > 1 for j = 1..k}`. For `k >= 2`, `Bad_k` lies in `3 (mod 4)`. For `x` outside `Bad_k`, `d(x) = min{j <= k : M_j(x) < 1}`; equality `M_j = 1` never occurs.

## 2. Theorem 1 -- flip the undecided residues

Let `sigma_k = -` on `Bad_k` and `+` elsewhere, and write `T_k = T_(sigma_k)`.

**Lemma 1.1 (heredity).** If `x` is outside `Bad_k`, then `T^i x` is outside `Bad_k` for `0 <= i < d(x)`.

*Proof.* Suppose `T^i x` lies in `Bad_k` for some `1 <= i < d(x)`. Then `M_i(x) > 1`, and `M_(i+t)(x) = M_i(x) M_t(T^i x) > 1` for `1 <= t <= k`. So `M_j(x) > 1` for every `j <= i + k`, which contradicts `M_(d(x))(x) < 1` with `i < d(x) <= k`. ∎

**Lemma 1.2 (blocks).** Every `x` in `Z_2` falls in exactly one case.
* **(F)** `x` lies in `Bad_k`. Then `T_k x = (3x-1)/2` is even and `T_k^2 x = (3x-1)/4`. This is two steps with multiplier `3/4`.
* **(C)** `x` lies outside `Bad_k`. Then `T_k^j x = T^j x` for `j <= d(x)`: the block follows Collatz along a first-descent word of length `d(x) <= k`, with multiplier `< 1`.

*Proof.* (F): `x = 3 (mod 4)` gives `3x - 1 = 0 (mod 4)`. (C): by Lemma 1.1, every odd point of the block before time `d(x)` carries `sigma_k = +`. ∎

**Theorem 1.** For every `k >= 2`:
* **(a)** Every `T_k`-orbit in `Z_2` is a concatenation of (F)- and (C)-blocks, each of length `<= k` and multiplier `< 1`. So `sigma_k` is class (i).
* **(b)** `rho_max(G_(sigma_k)) = F_k := max{a/d : d <= k, 3^a < 2^d}`, the best lower approximation of `log_3 2` with denominator `<= k`.
  * Upper bound: a cycle carries a periodic point (THM-4474, Lemma 2), whose period is a union of whole blocks. So its density is a mediant of block densities, each `<= F_k`.
  * Attainment: the word with prefix counts `ceil(c j)` for `j < d` and `a` at `j = d` is the upper Christoffel word of slope `F_k = a/d`. It is a first-descent word; its periodic point `c_w/(2^d - 3^a)` descends within `d` steps from every point, so it avoids `Bad_k` and is a `T_k`-cycle.
* **(c)** Every integer `n > n_0(k) = max over first-descent words w of floor(c_w/(2^d - 3^a))` satisfies `T_k^j(n) < n` for some `j <= k`.
* **(d)** `delta_k <= |Bad_k| <= sum_(j > ck) C(k,j) <= 2^(hk)`, so `delta_k/2^(k-1) <= 2^(1-(1-h)k) -> 0`.

*Why the rewiring does not matter.* A flip at `r` moves the target from `T(r)` to `T(r) - 1`, which rewires `G_sigma` globally; this was the obstacle named in HYP-9138. But an orbit meets a flipped residue only at the start of a block, because (C)-blocks avoid `Bad_k` (Lemma 1.1). At a flipped residue the orbit immediately takes a `3/4` descent. No block sees the rewiring.

## 3. Theorem 2 -- the matching lower bound

In parity-word coordinates, Collatz's parity graph is the de Bruijn graph `B(2,k)`: `v_0...v_(k-1) -> v_1...v_(k-1) b`, and a node is odd iff `v_0 = 1`.

**Theorem 2.** `delta_k >= N_k >= (1/k) C(k, ceil(ck)) >= 2^(hk)/(3k^2)` for every `k >= 2`.

*Proof.*
* The rotations of a `k`-word form a closed walk of length `k` in `B(2,k)`. Its density is the word's density, and distinct necklaces are node-disjoint.
* If `sigma` agrees with `+` at every odd node of an expanding necklace, the walk survives in `G_sigma`, because edges out of even nodes never change. So a class-(i) strategy changes a sign at an odd node of every expanding necklace, and these residues are distinct.
* Each necklace has at most `k` words, and `C(k,m) >= 2^(k h(m/k))/(k+1)`. For `k >= 10`, `h(ceil(ck)/k) >= h(c) - 1.4415/k`, because `|h'| <= log_2(0.7309/0.2691) = 1.4415` on `[c, c + 0.1]`. So `N_k >= 2^(hk)/(2.717 k(k+1)) >= 2^(hk)/(3k^2)`. The cases `k <= 9` are checked directly. ∎

**Corollary (sharp exponent).** `2^(1-(1-h)k)/(3k^2) <= delta_k/2^(k-1) <= 2^(1-(1-h)k)`.

## 4. Finite-exact data

| `k` | `N_k` | `delta_k` | `\|Bad_k\|` | Haar `delta_k/2^(k-1)` |
|---|---|---|---|---|
| 2 | 1 | 1 | 1 | 0.500 |
| 3 | 2 | 2 | 2 | 0.500 |
| 4 | 2 | 2 | 3 | 0.250 |
| 5 | 2 | 4 | 4 | 0.250 |
| 6 | 5 | 5 | 8 | 0.156 |
| 7 | 5 | 9 | 13 | 0.141 |
| 8 | 6 | 14 | 19 | 0.109 |
| 9 | 16 | 23 | 38 | 0.090 |
| 10 | 19 | 40–44 | 64 | 0.078–0.086 |
| 11 | 52 | 52–72 | 128 | 0.102–0.141 |

* **The optima are not `Bad_k`.** The optimal flip sets use 1.4–2.3 flips per expanding necklace. Some flips lie outside `Bad_k` (5 of 14 at `k = 8`), and a few are `d -> u` flips at residues `1 (mod 4)`. `Bad_k` is a proof device, not the shape of the optimum.
* **One flip per expanding necklace is not enough.** Flipping only the lowest, or only the highest, ballot rotation of each expanding necklace fails for `k = 5` and for every `7 <= k <= 14`.
* **Transitivity.** `n_0(k)` changes only at the denominators `2, 5, 8, 27, 46, 65, 149, 233` of the `F_k`. For every `k <= 300`, every `n <= n_0(k)` reaches 1 under `T_k`, so every positive orbit of `T_k` reaches 1.
* **SHEET.** Negation carries everything to 3n−1: the 3n−1 distance is the same `delta_k`, and `-Bad_k(3n+1) = Bad_k(3n-1)`.

## 5. DRIFT: the contrast with 5n+1

* The 5n±1 provable class is **empty** at levels 2–6 (UNSAT certificates re-checked by a SAT solver; exhaustive for `k <= 5`) and nonempty from level 7.
* 5n+1's exact distance at `k = 7` is `29/64`. Pruned class-(i) sets give `k` times the Haar fraction about `3.2` for `7 <= k <= 14`, against a necklace lower bound of about `2/k`.
* Theorem 1 has no analogue: 5n+1's undecided set has positive density (about `0.176`), and its flip block `(5r-1)/4` has multiplier `5/4 > 1`.
* **Proposition 4 (PROVED).** Every closed class of a class-(i) 5n±1 strategy has stationary mass at least `(1/2 - log_5 2)/k` on the changed residues.
* Whether 5n+1's distance tends to 0 is **OPEN**. Positive drift changes the rate from exponential to at best polynomial; whether it changes the limit is unknown.

## 6. Meaning

* **One exponent, three mechanisms.** The price of bounded-lookahead provability now has the same sharp exponent `1 - h(log_3 2) = 0.0500445` in three settings:
  * the pairing family: THM-4475 from above, THM-4478 from below;
  * arbitrary fixed-horizon edits: THM-4478;
  * the strategy cube: this theorem.

  The three lower bounds use three different mechanisms:
  * **Integer capacity** (THM-4478): exact affine spacing of actual integer ancestors.
  * **Moments** (THM-4477): distribution-only arguments, which stop at `2(1-h)`.
  * **Necklace packing** (here): available because the cube is periodic, so a finite graph supplies disjoint expanding cycles.
* **What it does not do.** `sigma_k` differs from Collatz exactly on the undecided classes, which is exactly where Collatz's behaviour is unknown. So the theorem says Collatz lies in the Haar closure of the provable strategies, and nothing about Collatz itself. Collatz stays in class (iv) at every level (THM-4474).
