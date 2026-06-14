---
id: THM-499
title: The spectral-reframe boundary — cycle counts are power sums, but H is strictly finer than the eigenvalue spectrum (cospectral tournaments with distinct H from n=6)
status: PROVED (exhaustive, exact integer spectral signatures + exact H; explicit cospectral/distinct-H witnesses)
source: kind-pasteur-2026-06-13-S5
depends_on:
  - THM-118   # c_k = tr(A^k)/k for k<=5 (cycle counts ARE power sums)
  - THM-498   # cycle-spectrum gaps = power-sum exclusions
related:
  - THM-029   # H-forbidden 7
  - THM-079   # H-forbidden 21
  - THM-133   # spectral-OCF chain
  - HYP-2492  # cycle gaps = skew-spectral exclusions (this REFINES its overreach to H)
---

# THM-499 — where the spectral reframe stops

This session's reframe (THM-498/HYP-2492) showed the cycle-count gaps are
**skew-spectral exclusions**: `c_k = tr(A^k)/k` for `k ≤ 5` (THM-118), so a
forbidden cycle-count is a non-realizable power-sum value. HYP-2492 conjectured
the OCF forbidden values `H ∈ {7,21}` are *also* power-sum exclusions. This file
draws the exact boundary: **they are not** — `H` is a strictly finer invariant
than the eigenvalue spectrum.

## Statement

Let `sig(T) = (tr A, tr A², …, tr Aⁿ)` (the exact integer power-sum signature ⟺
the characteristic polynomial of `A` ⟺ the eigenvalue spectrum).

> **THM-499.**
> 1. `H` (the Hamiltonian-path count = the OCF `I(Ω,2)`) **is** determined by
>    `sig(T)` at `n ≤ 5`, but **is NOT** at `n ≥ 6`: there exist **cospectral
>    tournaments with distinct `H`**. Witnesses (n=6, exhaustive): the spectral
>    class `sig = (0,0,18,28,30,120)` carries `H ∈ {25,29,33}`; `(0,0,12,12,10,48)`
>    carries `H ∈ {13,17}`; `(0,0,21,36,35,159)` carries `H ∈ {33,37}`.
> 2. Even the odd-cycle-count pair `(c3,c5)` does not determine `H` (6 split
>    classes at n=6).
> 3. Consequently the forbidden OCF values `{7,21}` are **NOT power-sum
>    exclusions** — unlike the cycle gaps (`c5=10`), they live in the
>    conflict-graph **disjointness** layer (`α₂` = # vertex-disjoint odd-cycle
>    pairs in `H = 1 + 2α₁ + 4α₂ + …`), a strictly higher invariant than the
>    spectrum.

*Status:* exhaustive over all `2^{C(n,2)}` labeled tournaments at `n=5,6`; `sig`
computed as exact integer traces, `H` exact (subset-DP Hamiltonian-path count);
witnesses explicit. (`04-computation/spectral_reframe_boundary_kps5.py` + `.out`.)

## The exact mechanism (PROVED): `H = 1 + 2(c3+c5) + 4D`, and `D` is the first non-spectral ingredient

At `n ≤ 6` the OCF `H = I(Ω,2) = Σ_k α_k 2^k` truncates exactly: the only odd
cycles are 3- and 5-cycles, and a vertex-disjoint *pair* of odd cycles must be two
triangles filling all 6 vertices (`α_3` needs ≥ 9 vertices). So

> **`H(T) = 1 + 2(c3 + c5) + 4D`**, where `α_1 = c3 + c5` (total odd cycles) and
> `D = α_2 = #{unordered pairs of vertex-disjoint triangles}`.

*Verified exhaustively, exact: holds for ALL `2^{10}` (n=5) and `2^{15}` (n=6)
tournaments, 0 exceptions* (`04-computation/ocf_first_nonspectral_ingredient_kps5.py`).
Now the boundary is mechanical and sharp:

- `c3 = tr(A³)/3` and `c5 = tr(A⁵)/5` are **spectral** (THM-118), so `α_1` is fixed
  within a cospectral class.
- `D = α_2` is **NOT spectral**: it varies within exactly the 3 cospectral classes
  that split `H` at `n=6`. The witness `sig=(0,0,18,28,30,120)` has `c3=6, c5=6`
  (so `α_1=12`, `1+2·12=25`) and `D ∈ {0,1,2}`, giving `H = 25 + 4D ∈ {25,29,33}` —
  exactly the observed split.
- At `n=5`, `D ≡ 0` (two disjoint triangles need 6 vertices), so `H = 1+2(c3+c5)`
  is fully spectral — which is *why* `H` is spectrally determined at `n ≤ 5`.

> **The spectral-reframe boundary IS the onset of `α_2 = D` at `n = 6`.** The OCF
> is spectral exactly up to the first vertex-disjoint pair of odd cycles; that pair
> is the first piece of conflict-graph data the eigenvalue spectrum cannot see.

This also locates it on the OCF's own structure: `α_1` (odd-cycle count) is the
power-sum/spectral layer; `α_2, α_3, …` (the disjointness / independence-polynomial
layers of `Ω`) are the strictly-finer conflict-graph layers, and the forbidden
values `H ∈ {7,21}` (which need `α_1=3 ⟹ α_2≥1`, THM-029/079) live in the `α_2`
layer — *not* the spectrum.

## The two kinds of forbidden value (the refined picture)

| forbidden value | layer | reason it's forbidden |
|---|---|---|
| `c5 = 10` (n=6) | **spectral** (power sum `tr A⁵`) | no tournament spectrum realizes `tr A⁵ = 50` (THM-498) |
| `H ∈ {7,21}` | **conflict-graph / disjointness** (`α₂`) | `α₁=3 ⟹ α₂≥1` etc. (THM-029/079) — NOT a spectral statement |

So "every forbidden value is a power-sum exclusion" is FALSE. The OCF invariant
hierarchy stratifies by *which* data the value depends on: the low cycle counts
are pure spectrum (`tr A^k`), but `H` injects the disjointness structure of the
odd-cycle conflict graph `Ω`, which cospectral tournaments can differ in.

## The invariant map across the boundary (the determinant lens is on the SPECTRAL side)

Sorting the standard tournament invariants by which side of the boundary they sit:

| invariant | spectral? | evidence |
|---|---|---|
| `c3, c4, c5` | **spectral** | `= tr(A^k)/k` (THM-118) |
| `d = det(I+S)` | **spectral** | `det(I+S) = ∏(1+μ_j²)` (skew spectrum); constant on every A-cospectral class (0/28 split at n=6) — verified |
| `H = I(Ω,2)` | **NOT spectral** | `= 1+2(c3+c5)+4D`; `D=α_2` varies within cospectral classes (THM-499) |
| `c7, α_2, α_3, …` | **NOT spectral** | first compound walk (k=6) / disjointness layers |

So the repo's empirical finding that **`d ⊥ H`** (Pearson ≈ 0, the determinant lens,
THM-468/472/473) gains a structural cause: `d` is a *spectral* coordinate (a
function of the skew spectrum `{μ_j²}`, which the A-spectrum fixes at n=6) and `H`
is the *non-spectral* coordinate (it reads `α_2 = D`). They are orthogonal because
they live on **opposite sides of the spectral boundary** — `d` sees only the
spectrum, `H` sees only past it. (Verified `04-computation/cospectral_det_check_kps5.py`.)

## A worked proof on the spectral side: the `c5=10` gap (score-stratified)

The spectral exclusion `c5=10` (THM-498) is now PROVED by a finite certificate (not
the `2^{15}` brute force): stratify by the 23 Landau score sequences at n=6. Every
class with `c3 ≤ 7` has `c5 ≤ 9`; the unique `c3=8` class — the regular tournament
score `(2,2,2,3,3,3)` — achieves `c5 ∈ {6,8,11,12}`, **skipping 10** (and 7, 9). No
score class realizes `c5=10`. ∎ (`04-computation/c5_gap_score_stratification...out`.)
This is "efficiency becomes proof" finished: the spectral-side gap reduces to a
human-checkable finite stratification; the non-spectral H-gaps `{7,21}` do not — they
need the `α_2` argument (THM-029/079), per the boundary.

## Why this matters (the reframe map gains an edge)

- `H` is **strictly finer** than the eigenvalue spectrum: cospectral tournaments
  (eigenvalue twins) are separated by `H`. This sharpens the repo's "H as a
  universal tournament fingerprint" (engineering domain 12): the fingerprint
  succeeds *because* it sees past the spectrum into `Ω`.
- It tells the spectral-exclusion program (HYP-2492, the route to proving `c5=10`
  and beyond) exactly how far it reaches: it can prove the **cycle-count** gaps
  spectrally, but the **H** gaps need a conflict-graph (independence-polynomial)
  argument — the two are genuinely different theorems, and conflating them was the
  HYP-2492 overreach this file corrects.
- The cospectral/distinct-H pairs are concrete new objects: minimal "eigenvalue
  twins, OCF-distinct" tournaments, the n=6 frontier of spectral non-determination
  for the Rédei–Berge / Hamiltonian-path statistic.

**Artifacts:** `04-computation/spectral_reframe_boundary_kps5.py` (+ `.out`).
Reflection: extends `efficiency-becomes-proof-kps4` (the spectral reframe and its
boundary).
