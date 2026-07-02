# THM-606: The depth-d ladder (the multi-level nested certificate, with feasibility calculus)

**Status:** PROVED (elementary; downward induction; exact-rational depth-3 instantiation verified —
see `04-computation/lrc14_depth_d_ladder_verify_opus_20260702_S36.py` + `.out`)
**Author:** opus-2026-07-02-S36 (HYP-3905)
**Generalizes:** kps `cert_two_level` (HYP-3959; = the d = 2 instance, hypothesis-for-hypothesis)
and THM-603's bridging step; the "multi-level recursion glue" of the 3-item DAG surface.

---

## Setting

Band `h > 0` (here `h = 1/14`). Window `[lo, hi] ⊂ ℝ`. Small speeds `P ⊂ ℤ_{≥0}`. Levels
`ℓ = 1..d`, each with an integer reference speed `V_ℓ > 0`, offsets `offs_ℓ ⊂ ℤ_{≥0}`
(runners `V_ℓ − o`), a phase target `c_ℓ ∈ ℚ`, and a drift budget `μ_ℓ ≥ 0` with **`μ_d = 0`**.

Write `δ_0 := hi − lo` and `δ_ℓ := μ_ℓ / V_ℓ` (so `δ_d = 0`), and let

    Safe(β, w, c)  :⟺  ∀ τ ∈ [lo, hi] : dist(w·τ − c, ℤ) ≥ β

(the semantic content of kps's decidable `arcSafe β w c lo hi`; `Safe(β, s, 0)` for small speeds).

## Hypotheses

- **(H-P)** every `s ∈ P` has `Safe(h, s, 0)`;
- **(H-ℓ)** for `ℓ = 1..d`: every `o ∈ offs_ℓ` has `Safe(h + μ_ℓ, o, c_ℓ)`;
- **(H-len)** for `ℓ = 1..d`: `1 < V_ℓ · (δ_{ℓ−1} − δ_ℓ)`.

## Theorem

Under (H-P), (H-ℓ), (H-len) there is `τ ∈ [lo, hi]` with `dist(s·τ, ℤ) ≥ h` for every `s ∈ P`
and `dist((V_ℓ − o)·τ, ℤ) ≥ h` for every `ℓ` and `o ∈ offs_ℓ`. (At `h = 1/14` and 13 total
runners: a `1/14`-lonely time for the multi-cluster family — for EVERY tuple `(V_1, …, V_d)`
satisfying (H-len).)

## Proof (downward induction; the carried variable is a REAL window start)

For `k = 1..d+1` let

    CLAIM(k): for every a ∈ ℝ with lo ≤ a and a + δ_{k−1} ≤ hi, there is
              τ ∈ [a, a + δ_{k−1}] with dist((V_ℓ − o)τ, ℤ) ≥ h  ∀ ℓ ≥ k, ∀ o ∈ offs_ℓ.

**CLAIM(d+1)** is trivial: the window is `[a, a]` (as `δ_d = 0`); take `τ = a`; no levels remain.

**CLAIM(k) from CLAIM(k+1).** The interval `[V_k·a − c_k, V_k(a + δ_{k−1} − δ_k) − c_k]` has
length `V_k(δ_{k−1} − δ_k) > 1` by (H-len), so it contains an integer `j_k`; set
`t_k := (j_k + c_k)/V_k ∈ [a, a + δ_{k−1} − δ_k]`. Then `[t_k, t_k + δ_k] ⊆ [a, a + δ_{k−1}]
⊆ [lo, hi]`, so CLAIM(k+1) applies at `a' := t_k` and returns `τ ∈ [t_k, t_k + δ_k]` handling
all levels `≥ k+1`. For level `k` itself: `V_k·τ − c_k = j_k + ε` with
`ε := V_k(τ − t_k) ∈ [0, V_k δ_k] = [0, μ_k]`, so for `o ∈ offs_k`

    dist((V_k − o)τ, ℤ) = dist((o·τ − c_k) − ε, ℤ) ≥ dist(o·τ − c_k, ℤ) − μ_k
                        ≥ (h + μ_k) − μ_k = h,

using 1-Lipschitzness of `dist(·, ℤ)` and (H-k) at `τ ∈ [lo, hi]`. ∎

The theorem is CLAIM(1) at `a := lo`, plus (H-P) at the resulting `τ ∈ [lo, hi]`. ∎

**Sanity (d = 2).** (H-len)₁ is `1 < V₁(hi − lo − μ₁/V₁)`, (H-len)₂ is `1 < V₂·μ₁/V₁`, (H-1) is
`arcSafe (h+μ) o c₁`, (H-2) is `arcSafe h o c₂` — *verbatim* the hypotheses of kps's
`cert_two_level`. The Lean wrapper is therefore a structural induction whose step is exactly the
existing two-level proof body with `(lo, hi)` generalized to `(a, a + δ_{k−1})`.

## The structural point: inflation does NOT accumulate

Level `k`'s band is `h + μ_k` — **not** `h + Σ_{m ≥ k} μ_m`. Each level's drift is measured from
its OWN ruler `t_k`, and the nested windows `W_d ⊆ … ⊆ W_k` bound it by `V_k δ_k = μ_k` alone.
Deeper levels move the final time only *within* windows already paid for. This depth-uniformity
is what lets the ladder serve as the induction engine of the all-`n` pipeline (HYP-3860,
THM-603): the certificate cost per level is a function of that level's scale ratio only, so the
induction hypothesis never degrades with depth.

## Feasibility calculus (what parameter tuples exist)

- **(F0) single-cell bound (inherited from the certificate form).** `Safe(β, w, c)` on a window
  of width `δ_0` forces `w·δ_0 ≤ 1 − 2β`: every level's offsets obey
  `max(offs_ℓ) ≤ (1 − 2(h + μ_ℓ))/δ_0`, uniformly in `ℓ` (full-window certificates are the
  price of data-independent decidability; wide-spread deep clusters need shape SPLITTING, not a
  deeper ladder).
- **(F1) telescoping.** (H-len) ⟺ `δ_{ℓ−1} > 1/V_ℓ + δ_ℓ`; hence `δ_0 > Σ_ℓ 1/V_ℓ` is
  necessary, and the canonical schedule `δ_ℓ := (1+η)·Σ_{m>ℓ} 1/V_m` (any `η > 0`) satisfies
  (H-len) with margin. The budgets are then `μ_ℓ = (1+η)·V_ℓ·Σ_{m>ℓ} 1/V_m`.
- **(F2) ratio form.** If consecutive references grow by `ρ` (`V_{ℓ+1} ≥ ρ V_ℓ`) then
  `μ_ℓ ≤ (1+η)/(ρ−1)`: **band inflation per level ≈ 1/(ratio − 1)**. Budget `μ* ⟹ ρ ≥ 1 + (1+η)/μ*`.
  At `h = 1/14` with `μ* = h/2 = 1/28`: `ρ ≈ 29`.
- **(F3) route comparison.** The measure-route gap peel (HYP-3900) needed `Λ ≈ 10⁴` for its
  error budgets; the witness ladder needs `ρ ≈ 30`. A point-witness pays the scale gap once
  (Lipschitz drift), a measure bound pays arc-count × equidistribution — ~300× more separation.
  This is why the witness route's finite windows (ratios below `ρ`) are small enough that the
  star-census reached `V* ≤ 67` (kps HYP-3960), and why the witness route is the right chassis
  for the formal proof.
- **(F4) merging.** Clusters with ratio `< ρ` merge into one level with joint offsets (subject
  to (F0)); the ladder degrades gracefully to smaller `d`. With 13 runners, `d ≤ 13` always,
  `d ≤ 6` for clusters of size ≥ 2 — and (F1) is then automatic for `V_1 ≥ 15`, `δ_0 ~ 1/50`.
- **(F5) consumption.** The shape universe becomes: `(P; scale-ordered clusters with ratios ≥ ρ)`
  — the witness-route mirror of the census-exhaustiveness case tree (opus-S33 doc): ladder =
  the gap-split case; below-ratio configurations = finite windows (kps star-census + opus
  window modules); base = bounded census. Same tree, exponentially cheaper constants.

## Verification

`lrc14_depth_d_ladder_verify_opus_20260702_S36.py`: exact-rational depth-3 instantiation on a
13-runner three-cluster family (P = {1,2}; clusters at ratios ~40), certificates found by exact
search, the ladder walk executed symbolically, all 13 runners verified `≥ 1/14` at the
constructed rational τ, robustness over randomized window starts `a`, and the (F1)/(F2) tables.
See the `.out` for the exact witnesses and margins.

-> HYP-3959/3960 (kps), THM-603, HYP-3860 (all-n pipeline), HYP-3900/3904 (opus), OPEN-Q-108.
