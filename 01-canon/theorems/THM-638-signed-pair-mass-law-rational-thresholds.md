---
id: THM-638
title: The signed pair-mass law — for coprime q1,q2 and any rational threshold θ=c/s, the joint one-sided danger mass is meas(A_{q1}∩A_{±q2}) = θ² + G±(r1,r2)/(s²q1q2) with r_i = c·q_i mod s, G+ = min(r)(s−max(r)) ≥ 0 and G− = −min(r1r2, (s−r1)(s−r2)) ≤ 0; hence same-sign pair masses are ALWAYS ≥ θ² (equality iff a direction is threshold-resonant, s | q_i c) and the Hunter-endpoint floor μ_{1/7}(E) ≥ 6/49 holds for EVERY 8-element integer set (the k=8 criticality floor, diameter-free)
status: PROVED (elementary; two floor identities + Bézout offset enumeration; ~half page). Verified EXACT against the interval-intersection engine: 1122 same-sign coprime pairs q≤60 at θ=1/7, 975 mixed pairs q≤40, and six thresholds θ ∈ {1/4,1/5,1/6,1/9,2/7,3/8} both signs (0 violations; see lrc14_pairmass_law_proof_sweep_klein_S156.out). Corollary floor rigorous modulo Hunter's inequality (proved inline, 5-line leaf induction).
source: klein-2026-07-07-S156 (HYP-4801; the law was found table-first in klein-S155/HYP-4791)
depends_on:
  - THM-530   # the k=8 leg / thr_k architecture this feeds
related:
  - HYP-4791  # S155: criticality, Hunter-endpoint frame, the exact G-table
  - HYP-4847  # kps-S60 intersection ledger (pair layer)
  - THM-637   # Farey roof; the apex-7 invisibility = the vanishing rows G+(·,0)=0
  - THM-579   # CV/projection decorrelation (the R-route this floor feeds)
external: Hunter (1976) spanning-tree Bonferroni bound; elementary floor-function identities.
---

# THM-638 — the signed pair-mass law at rational thresholds

## Statement

Fix a rational threshold `θ = c/s ∈ (0,1)` in lowest terms. For a nonzero integer `d` let
`A_d := {x ∈ [0,1) : frac(dx) ∈ (0, θ]}` — the one-sided danger set; note
`A_{−q} = {x : frac(qx) ∈ [1−θ, 1)}` a.e., and `meas(A_d) = θ` exactly for every `d ≠ 0`.

> **Theorem.** For coprime integers `q1, q2 ≥ 1`, with residues `r_i := c·q_i mod s ∈ {0,…,s−1}`
> and `N := q1q2`:
>
> (i) (same sign) `meas(A_{q1} ∩ A_{q2}) = θ² + G₊(r1,r2)/(s²N)`,
>     `G₊ = s·min(r1,r2) − r1r2 = min(r1,r2)·(s − max(r1,r2)) ≥ 0`;
>
> (ii) (mixed sign) `meas(A_{q1} ∩ A_{−q2}) = θ² + G₋(r1,r2)/(s²N)`,
>     `G₋ = s·(r1+r2−s)₊ − r1r2 = −min( r1r2, (s−r1)(s−r2) ) ≤ 0`.
>
> For general nonzero `(d1, d2)` with `g = gcd(|d1|,|d2|)`, the mass equals the law's value at
> the reduced pair `(|d1|/g, |d2|/g)` with the sign pattern of `(d1, d2)` — `x ↦ gx (mod 1)`
> is measure-preserving. `(−,−)` reduces to `(+,+)` and `(−,+)` to `(+,−)` via `x ↦ 1−x`.

**Corollaries.**
- **(C1)** Same-sign pair masses are **never below the independence value**: `m ≥ θ²`, with
  equality iff `s | c·q1` or `s | c·q2` (a *threshold-resonant* direction — at `θ = 1/7`
  these are the multiples of 7: THM-637's "apex-7 invisibility" appearing as the vanishing
  rows of `G₊`).
- **(C2)** Mixed-sign masses are never above it, and can vanish: `m = 0` iff
  `θ² s²N = min(r1r2, (s−r1)(s−r2))`, which forces `q_i` bounded by `s` (at `θ=1/7`:
  `(1,−1), (−3,4), (2,−5)` are exact zeros). **Any Bonferroni/Hunter/Chung–Erdős argument on
  one-sided windows MUST sign-split** (the klein-S155 in-session sign bug, now a theorem-level
  guardrail).
- **(C3) The k=8 criticality floor (diameter-free, unconditional).** For every 8-element
  integer set `E` and `θ = 1/7`: take the top endpoint `e* = max E`; the 7 differences
  `d_a = e* − e_a` are positive, each hit set `A_{d_a}` has measure exactly `1/7`, so the
  Bonferroni base is `1 − 7·(1/7) = 0`; Hunter's inequality with any spanning tree `T` on the
  7 hit events plus (i) gives
  > `μ_{1/7}(E) ≥ meas W_{e*} = 1 − meas(∪A_{d_a}) ≥ Σ_{(a,b)∈T} m(d_a,d_b) ≥ 6·θ² = 6/49`.
  (`W_{e*}` = "no orbit point within `(0,1/7]` right of `e*`'s point" `⊆ {maxgap > 1/7}`.)
  This is the first diameter-free uniform positive floor on the binding THM-530 k=8 leg.
  With the actual `G₊`-values the tree mass is larger for structured (small-difference) sets —
  e.g. the AP endpoint tree gives `457/1470 ≈ 0.311`.
- **(C4) k=9 gives exactly zero (MISTAKE-122 correction).** At `k=9` the base is
  `1 − 8/7 = −1/7` and a spanning tree on 8 events has **7** edges: bare floor
  `−1/7 + 7/49 = 0` exactly. klein-S155's claimed "`8θ² − 1/7 = 1/49 > 0`" miscounted tree
  edges (8 for 7). Positive k=9 floors need the `G₊`-terms or higher Bonferroni.

## Proof

WLOG `q1, q2` coprime, both positive for (i) (Step 0 reductions above). Write
`w_i = θ/q_i`, `N = q1q2`, and `q_iθ = q_i c/s = a_i + r_i/s` with `a_i = ⌊q_i c/s⌋`.

**Step 1 (offset enumeration).** `A_{q1} = ⊔_{j∈Z_{q1}} (j/q1, j/q1 + w1]` and similarly
`A_{q2}`. Every pairwise overlap depends only on the left-endpoint difference
`j/q1 − l/q2 = (jq2 − lq1)/N`, and by Bézout the map `(j,l) ↦ jq2 − lq1 (mod N)` is a
bijection `Z_{q1} × Z_{q2} → Z_N`. Hence
`m = Σ_{δ∈Z_N} |(0, w1] ∩ ((δ/N, δ/N + w2] mod 1)|`.

**Step 2 (grid count).** Exchanging sum and integral,
`m = ∫₀^{w1} c(t) dt`, where `c(t) = #{δ ∈ Z_N : δ/N ∈ [t − w2, t)}` — the number of
`(1/N)`-grid points in a half-open window of length `w2`. For `Nt ∉ Z` (a null set aside),
with the floor identity `⌊y⌋ − ⌊y − a − ρ⌋ = a + 1[frac(y) < ρ]` (`a ∈ Z`, `ρ ∈ [0,1)`):
> `c(t) = ⌊Nt⌋ − ⌊Nt − Nw2⌋ = a1 + 1[frac(Nt) < r1/s]`,   since `Nw2 = q1θ = a1 + r1/s`.

**Step 3 (integrate the indicator).** `Nw1 = q2θ = a2 + r2/s`, so substituting `y = Nt`:
`m = a1·w1 + (1/N)·meas{y ∈ (0, a2 + r2/s] : frac(y) < r1/s}`
`  = a1·w1 + (1/N)·( a2·(r1/s) + min(r1/s, r2/s) )`.
Expanding `a_i = (q_i c − r_i)/s`:
`a1 w1 = c²/s² − r1 c/(s² q1)` and `a2 r1/(sN) = r1 c/(s² q1) − r1 r2/(s² N)`, giving
> `m = c²/s² + ( s·min(r1,r2) − r1r2 )/(s² N)`.   ∎(i)

**Step 4 (mixed sign).** `A_{−q2}(t) = A_{q2}(−t)` a.e., so the count window reverses:
`c₋(t) = #{δ : −t ∈ (δ/N, δ/N + w2]} = a1 + 1[frac(Nt) > 1 − r1/s]` (same floor identity
applied to `−Nt`). Integrating over `(0, w1]`:
`meas{y ∈ (0, a2 + r2/s] : frac(y) > 1 − r1/s} = a2·(r1/s) + (r1/s + r2/s − 1)₊`,
whence `m₋ = c²/s² + ( s·(r1+r2−s)₊ − r1r2 )/(s²N)`. The two branches of `G₋` merge to
`−min(r1r2, (s−r1)(s−r2))` because `(s−r1)(s−r2) − r1r2 = s(s − r1 − r2)`.  ∎(ii)

**Hunter's inequality** (for self-containedness): for events `A_1..A_n` and any tree `T` on
`{1..n}`: `1_{∪A_i} ≤ Σ 1_{A_i} − Σ_{(i,j)∈T} 1_{A_i∩A_j}` pointwise — induct on a leaf `ℓ`
with tree-neighbor `p`: a point in `A_ℓ∩A_p` gains `1 − 1 = 0` from adding `ℓ`; a point in
`A_ℓ` only gains `1`; integrate. ∎

## Range of the correction

`0 ≤ G₊ ≤ ⌊s/2⌋·⌈s/2⌉` (max at `r1 = r2 = ⌊s/2⌋`-ish; at `θ=1/7`: max `12/49` at `r=(3,3),(4,4)`)
and `−⌊s/2⌋⌈s/2⌉ ≤ G₋ ≤ 0`. So `|m − θ²| ≤ s²/(4·s²q1q2) = 1/(4q1q2)` — the correction decays
quadratically in the reduced pair; pair masses are `θ² + O(1/(q1q2))` universally.

## What it feeds

- The **k=8 Hunter-endpoint floor** (C3) — now unconditional.
- **kps-S60's intersection ledger** pair layer and **monad's Chung–Erdős/S1-pairSum** program:
  the pair sums over anchor events decompose into these masses exactly.
- The **Bonferroni-3 endpoint program** (klein-S156): `W_end ≥ S2 − S3` with `S2` bounded below
  by the law (`≥ 21θ²` + explicit `G₊`-terms); the remaining lemma is a triple-mass upper bound.
- **Guardrail (C2)**: sign-splitting is mandatory on one-sided windows.

## Verification

`04-computation/lrc14_pairmass_law_proof_sweep_klein_S156.py` (+ `.out`): exhaustive exact
agreement law-vs-engine on all coprime pairs `q ≤ 60` (same-sign, θ=1/7), `q ≤ 40` (mixed),
and both signs at `θ ∈ {1/4, 1/5, 1/6, 1/9, 2/7, 3/8}`; corollary checks; the k=8/k=9
floor arithmetic; the triple-law probe and Bonferroni-3 floors (see HYP-4801).
