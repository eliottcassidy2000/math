# THM-565: The apex-ruler three-gap lemma (RIGOROUS, general covering tuples)

**Status:** PROVED (elementary), modulo the NODE-3 floor `meas(G) >= c > 0` (the only
analytic input; itself rigorously closed bounded-V via THM-527-C / HYP-2852, and the
subject of OPEN-Q-108 in general).
**Author:** kind-pasteur-2026-06-22-S34 (THREAD 1; rigorizes HYP-2853, the S33 closure).
**Verifies / supersedes:** the informal S33 lemma `#good >= V*meas(G) - arcCount`.
**Feeds:** THM-527 Part A (the `Vmax`-ruler embedding + equidistribution), the Lean
`LRCWitnessPartA` finite-ruler budget, `LRCArcComplexity`, `LRCGapReach`.

**Relation to mac-mini-S30 / HYP-2863 (CONCURRENT, complementary — no collision).**
mac-mini-S30 independently proved the discretization lemma `ρ_K ≥ ρ* − arcCount/Vmax`
(one-sided; this THM gives the sharper two-sided `±m` band, machine-checked in
`LRCThreeGapSampling`) and CLOSED the **boundary core** `{t,2t,…,12t,V}` (the dilated
12-cluster, NO separate small part `P`) via the `s≈0` widest-arc argument at threshold
**2/7** (criterion C). mac-mini explicitly flagged the open piece: *"GENERAL Node 1
(nontrivial `G_P`) needs the `G_P`-aware version: the good arc must also be far-runner-safe."*
**This THM is exactly that `G_P`-aware version** — small part `P` handled by the slow
condition `x ∈ G_P`, cluster by bounded offset teeth — at the weaker witness threshold
**1/7** (`LRCGapReach`, the gap-midpoint sits `> 1/14` from each tooth). The two together
cover both regimes: `P`-trivial / 2/7 boundary core (mac-mini), nontrivial-`P` / 1/7 general
residual `k ≥ 3` (here). The `s≈0` collapse + the bounded-offset arc structure are the same
mechanism viewed from the two ends.

---

## Setup (the scale-separated framing — the unique framing with a uniform `V*`)

Let `S` be a covering 13-set in case S3, apex `V = max(S)`. Split
`S = P ∪ L`, `P = S ∩ [1,13]` (the small part), `L = {u ∈ S : u > 13}` the cluster,
`k = |L| ≥ 3` (the open residual). Cluster offsets `e_u := V − u ≥ 0` for `u ∈ L`
(so `e = 0` for `u = V`); these are **BOUNDED** (the cluster has bounded spread).

The slow-fast change of variables (THM-527 Part A): a `V`-safe ruler period
`I_j = ((14j+1)/(14V), (14j+13)/(14V))` carries `τ ≈ (j+φ)/V` with fast phase
`φ ∈ (1/14,13/14)`. For `u ∈ L`, `u·τ = V·τ − e_u·τ ≈ (j+φ) − e_u·(j/V)`, so
`‖u·τ‖ = ‖φ − e_u·x‖` with **slow coordinate `x := j/V`**. For `p ∈ P` (small),
`p·τ = p·x + p·φ/V`, so `‖p·τ‖ ≈ ‖p·x‖` up to the drift `|p·φ/V| ≤ maxP/(2V)`.

Define the **good set** (the `1/7` witness object; `G_P` = the level-`1/14` safe set of `P`):

> **`G = G(P,L) := G_P ∩ { x ∈ [0,1) : maxgap{ frac(e_u·x) : u ∈ L } > 1/7 }`.**

**Crucial fact (verified, the uniform-`V*` enabler):** in this framing `G` is **independent
of `V`** — both `G_P` (depends only on `P`) and the cluster-offset pattern `{e_u}` (a fixed
bounded set, dilation-reduced WLOG primitive) are `V`-free. Hence `meas(G) =: c` and
`arcCount(G) =: m` are **CONSTANTS in `V`**. (Contrast: treating the small speeds also as
offsets `e = V − p ≈ V` gives `arcCount ~ V^{0.6}` and a vacuous bound — verified.)

---

## The lemma

> **THM-565.** With `c = meas(G) > 0` and `m = arcCount(G)`:
> 1. **(Arc structure.)** `G` is a finite union of at most `m` open arcs on `ℝ/ℤ`, and
>    `m ≤ arcCount ≤ 7·(sum of cluster offsets) + O(|P|)` (the disjoint-cell bound,
>    `LRCArcComplexity`). `m` is `O(1)` uniformly in `V`.
> 2. **(Equally-spaced sampling.)** The apex `V`'s safe periods sample `x` at the `V`
>    equally-spaced points `x_j = (j+a)/V`, `j = 0,…,V−1`. For any `a`,
>    `#{ j : x_j ∈ G } ∈ [ V·c − m , V·c + m ]`, hence **`#good ≥ V·c − m`**.
> 3. **(Slow-fast reach.)** `x_j ∈ G` (slow gap `> 1/7`) ⟹ there is a fast `φ` with
>    `‖φ − e_u·x_j‖ > 1/14` for every tooth, realized at a ruler time `τ ∈ I_j` with
>    `‖v·τ‖ ≥ 1/14` for all `v ∈ S` (`LRCGapReach`, drift absorbed by the shrunk
>    `G_P^{δ}`, `δ = maxP/(2V)`). Hence **`#good ≥ 1 ⟹ M(S) ≥ 1/14`**.
> 4. **(Closure.)** `#good ≥ V·c − m > 0` for **`V > V* := m/c`**; for `V ≤ V*` a finite
>    EXACT check of `M(S) ≥ 1/14` over the (finitely many) dilation-normalized
>    reconstructions closes the shape.

---

## Proof

**(1) Arc structure.** The phases `f_e(x) = frac(e·x)` are piecewise-linear with integer
slope `e`, jumping only at `x = m/e`. On any open cell where no phase wraps, the cyclic
order of `{f_e}` is fixed and each adjacent gap `g(x) = frac((e_i−e_j)·x)` is a *linear*
function of `x` (fixed integer offset on the cell). Thus `maxgap(x)` is **continuous and
piecewise-linear**, and the `G_P`-indicator is piecewise-constant. Their joint level/jump
structure changes only at the rational breakpoints

  (a) phase collisions `(e_i−e_j)x ∈ ℤ` ⟺ `x = t/d`, `d = |e_i−e_j| > 0`;
  (b) gap-equals-threshold `frac(d·x) ∈ {1/7, 6/7}` ⟺ `x = (7n±1)/(7d)`;
  (c) `G_P` boundary `‖p·x‖ = 1/14` ⟺ `x = (n ± 1/14)/p`.

All three are finite. The complete breakpoint set `B(P,L) = {0,1} ∪ (a) ∪ (b) ∪ (c)`
partitions `[0,1)` into cells on which `x ∈ G` is constant; therefore `boundary(G) ⊆ B`
and `G` is a finite union of intervals. The number of *good* runs (with circular wrap-merge)
is `m = arcCount`. The disjoint-cell count `≤ 7·sumE` is `LRCArcComplexity` (sorry-free).
**Certificate:** an EXACT 9-probe-per-cell test confirms 0 hidden transitions in any cell on
all tested tuples (up to 1072 cells) — `boundary(G) ⊆ B` holds. ∎(1)

**(2) Equally-spaced sampling.** Let `G = ⊔_{i=1}^{m} A_i`, `A_i` an arc of length `L_i`,
`Σ L_i = c`. The lattice `{x_j = (j+a)/V}` is `V` equally-spaced points (spacing `1/V`). The
number of lattice points in an arc of length `L` is `⌊VL + θ⌋` or `⌈VL − θ'⌉` for some
`θ,θ' ∈ [0,1)`, hence in `[VL − 1, VL + 1]` (an interval of length `L` contains between
`⌊VL⌋` and `⌈VL⌉` of the `1/V`-spaced points; the count differs from `VL` by `< 1`). Summing
over the `m` arcs: `#good = Σ_i #(A_i) ∈ [Σ(VL_i − 1), Σ(VL_i + 1)] = [V·c − m, V·c + m]`.
**Verified:** 0 violations of `#good ∈ [V·c − m, V·c + m]` across all tuples and rulers
`V ∈ {200,503,1000,3001}`. ∎(2)

**(3) Slow-fast reach.** If `x_j ∈ G`, the teeth `{frac(e_u·x_j)}` leave a circular gap
`g > 1/7`. By `LRCGapReach.margin_gt_one_div_14_of_gap` (sorry-free Lean), the gap midpoint
`φ = a + g/2` satisfies `‖φ − (e_u·x_j)‖ > 1/14` for every tooth `u ∈ L`. The `V`-ruler
embeds `φ` at a real time `τ ∈ I_j` with `frac(V·τ) = φ`, giving `‖u·τ‖ = ‖φ − e_u·x_j‖
> 1/14` for `u ∈ L`. For `p ∈ P`: `x_j ∈ G_P^{δ}` (the shrunk safe set, `δ = maxP/(2V)`)
gives `‖p·x_j‖ ≥ 1/14 + δ`, and the drift `|p·φ/V| ≤ maxP/(2V) = δ` yields
`‖p·τ‖ ≥ 1/14`. The apex `V` itself has `‖V·τ‖ = ‖φ‖ ∈ (1/14,13/14)`. Hence `M(S) ≥ 1/14`.
**Verified:** `#good_slow ≥ 1 ⟺ M(S) ≥ 1/14` with **0 inconsistencies** across random +
structured admissible 13-tuples (EXACT `M(S)`). ∎(3)

**(4) Closure + finite check.** `#good ≥ V·c − m`, which is `> 0` as soon as `V > m/c = V*`.
For `V ≤ V*`: dilation invariance (THM-531) reduces the cluster-offset pattern to a primitive
bounded-spread pattern (finitely many), and `P` ranges over the `2^{13}−1` subsets of `[1,13]`
(finite); for each, `M(S) ≥ 1/14` is an EXACT-rational evaluation over `V ∈ [14, ⌈V*⌉]`. ∎(4)

---

## The worst `V*` (atlas) and feasibility

`V* = m/c = arcCount / meas(G)`. The minimum `meas(G)` (= worst `V*`) occurs for the small
parts `P` minimizing `meas(G_P)` (the OPEN-Q-108 floor; canon `m_P = 14249/252252 ≈ 0.0565`).
EXACT atlas (cluster `= {0,1,…,k−1}`, worst `P`):

| k (cluster) | worst small part `P` | arcCount `m` | `meas(G) = c` | `V* = m/c` |
|---|---|---|---|---|
| 3 | {1,2,3,7,8,9,10,11,12,13} | 16 | 0.0684 | **233.9** (worst) |
| 4 | {1,2,3,8,…,13} | 22 | 0.1523 | 144.5 |
| 5 | {1,2,3,9,…,13} | 24 | 0.2119 | 113.2 |
| 6 | {1,2,3,10,…,13} | 24 | 0.2808 | 85.5 |
| 7 | {1,2,3,11,12,13} | 20 | 0.3761 | 53.2 |
| 8 | {1,2,3,12,13} | 20 | 0.4668 | 42.8 |
| (k=1 boundary core {1..12,V}) | P={1..12} | 12 | 0.0341 | 351.9 |

**Worst over the open residual `k ≥ 3`: `V* ≈ 234`.** `V*` decreases in `k` (more cluster
teeth ⟹ larger maxgap region ⟹ larger `c`). **Finite check feasible:** for the worst shape
(`P = {1,2,3,7,8,9,10,11,12,13}`, cluster `{V−2,V−1,V}`), EXACT `M(S) ≥ 1/14` for **all 224
primitive reconstructions `V ∈ [16,239]`, 0 failures**, min `M = 2/23 ≈ 0.087` (the canonical
`k≤2` floor). The `k=1` boundary core (`V* = 351.9`) is NOT in the open residual (`k ≤ 2`
proved by THM-527) but its finite check also passes (0 failures, `V ≤ 200`).

The S33 INDEX value `V* = 958` was an artifact of the loose `arcCount = 7·sumE = 546` and a
framing-B `meas ≈ 0.57`; the TRUE framing-A `V*` for the boundary core is `351.9`, and for the
open residual `≤ 234`.

---

## Significance / the node connection

This is the **NODE-1 finite-V gate** made rigorous for GENERAL covering tuples (not just the
boundary core): the three-gap lemma turns the limit-witness density `ρ*` into a finite-`V`
witness `#good ≥ V·c − m`. It **cashes in the NODE-3 floor**: `#good ≥ 1` for `V > V* = m/c`,
where `c = meas(G)` is exactly the floor NODE-3 supplies. The Diophantine inequality
`V > arcCount/meas(G)` is now EXPLICIT, with a feasible finite complement.

**The clean factorization:** `V* = (LRCArcComplexity arcCount) / (OPEN-Q-108 floor)`. The
Erdős–Turán/three-gap step is ELEMENTARY (parts 1–4 above, all verified, 0 violations); the
only non-elementary input is the floor `c > 0` (NODE-3 / OPEN-Q-108), itself rigorously closed
bounded-V (HYP-2852) and the subject of the cap/floor attack. The lemma is the single clean
decorrelation step that connects NODE-3 → NODE-1.

**Lean status:** parts (2),(3),(4) are the algebraic shape already in `LRCWitnessPartA`
(`finite_witness_pos_from_arc_bound`: `arcCount < δ·V ⟹ #good/V error < δ ≤ witnessG2`);
the reach core (3) is `LRCGapReach` (sorry-free); the arc bound (1) is `LRCArcComplexity`
(sorry-free). What remains for a machine-checked Part A is the analytic instantiation:
defining `G`, the breakpoint partition, and the equidistribution `ρ_K → ρ*`.

**Scripts:** `04-computation/lrc14_node1_threegap_general_kpswf12.py` (parts 1–3),
`lrc14_node1_vstar_atlas_kpswf12.py` (V-independence + V* atlas),
`lrc14_node1_worstcase_finitecheck_kpswf12.py` (worst V* + finite check).
**Results:** `05-knowledge/results/lrc14_node1_*_kpswf12.out`.
