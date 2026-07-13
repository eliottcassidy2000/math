---
source: opus-2026-07-11-S263
status: REDIRECT (Gowers fails; the right tool is E3/additive). Bounding the multi-linear cancellation (S262
  residual) via Gowers norms FAILS: the linear forms v_i·t are all PARALLEL, so the generalized von Neumann
  theorem's general-position hypothesis is violated -- Gowers norms do not bound int prod g(v_i t). The
  multi-linear cancellation is instead governed by the ADDITIVE RELATIONS ±v±w_i±w_j=0 among the speeds -- the
  sumset / Schur = E3 structure (LEM-015), verified (corr(|eps_v|, #relations)=0.527, monotone). So the
  covering-min residual reduces to an E3 / dissociation bound: the coprime core is additively dissociated from
  the non-core => eps small; runner 1 = max E3 = near-AP = S255. The crux is the project's own additive
  invariant, not Gowers.
tags:
  - lrc14
  - covering-min
  - multilinear
  - gowers-norms
  - additive-relations
  - E3
  - redirect
---

# Gowers doesn't apply (parallel forms); the multi-linear cancellation is E₃-additive

**opus-2026-07-11-S263.** Owner: bound the multi-linear cancellation via Gowers norms. Trying it shows Gowers
norms are the *wrong* tool — and points at the *right* one, which the project already has.

## Why Gowers norms do not apply

The S262 residual is the multi-linear correlation of dilates of the band `g` at the **same** point:
`∫ ∏_i g(v_i t) dt = Σ_{Σ v_i k_i = 0} ∏_i ĝ(k_i)`. Gowers `U^s` norms bound a multilinear average
`∫ ∏_i g(ψ_i(t))` via the **generalized von Neumann theorem** — but *only* when the linear forms `ψ_i` are in
**general position** (pairwise linearly independent, bounded Cauchy–Schwarz complexity). Here
`ψ_i(t) = v_i·t` are all **proportional to `t`** — a single 1-dimensional direction, maximally degenerate.
General position fails completely, so the generalized von Neumann / Gowers-norm machinery gives **nothing**.
Gowers norms measure higher-order (polynomial) structure along genuinely multi-dimensional patterns; our
correlation is a **1-parameter, parallel-forms** object, which is a different (and in a sense simpler) beast.

## The right tool: additive relations (E₃)

The correlation `Σ_{Σ v_i k_i = 0} ∏ ĝ(k_i)` is dominated by the **`±1` relations**: `±v ± w_1 ± w_2 (…) = 0`,
i.e. the **additive relations** among `{v}` and non-core subsets (each contributes `~b_1^{|S|+1}`). This is
exactly the **sumset / Schur-triple = E₃ structure**, the project's own additive invariant (LEM-015).

**Verified:** `|ε_v|` correlates with the count of additive relations `v = w_i ± w_j` among non-core:

| #relations | mean `|ε_v|` |
|---|---|
| 0 | 0.021 |
| 2 | 0.027 |
| 4 | 0.065 |
| 6 | 0.073 |
| 8 | 0.086 |

`corr(|ε_v|, #relations) = 0.527`, monotone. So the multi-linear cancellation is driven by **how additively
related `v` is to the non-core speeds** — the E₃ structure, not Gowers. **Runner 1** has the *most* relations
(every consecutive non-core difference `w'−w = 1`), hence `ε_1 ≈ 0.57` — the near-AP maximum (LEM-015: E₃ is
*maximized* at the interval/AP), handled by **S255**. Dissociated core runners (large, coprime, few relations)
have `ε ≈ 0.02`.

## Net — the redirect

Bounding the multi-linear cancellation via **Gowers norms fails** (parallel forms violate general position).
The correct tool is **additive combinatorics**: the multi-linear cancellation *is* the additive-relation
(sumset / Schur / E₃) structure of the coprime core against the non-core — which the project already has as
**LEM-015 (E₃ maximized at the AP)**. So the covering-min residual reduces to an **E₃ / dissociation bound**:
the coprime core is additively **dissociated** from the non-core (few relations `±v±w_i±w_j=0`) ⟹ `ε_v` small
⟹ `coreCover < 1` ⟹ LRC(14), with the **AP / runner 1 (max E₃)** as the S255 exception. The whole S253–S263
arc thus closes a loop: the analytic drilling (balance → dual → mollification → completion identity →
multi-linear → Gowers) lands back on the fleet's **E₃ additive invariant** (LEM-015, opus-S246) — the crux is
the *dissociation of the core from the good set's additive structure*, maximized exactly at the AP that S255
already proves.

**Concrete next target:** an effective E₃/dissociation bound — `Σ_{core v} (#additive relations of v with the
non-core) · b_1^{…} < (6/7)^{core}` — using LEM-015's E₃-extremality (AP is the unique max) to control the
non-extremal covering families, with runner 1 / near-AP absorbed into S255. This is additive combinatorics on
the speed set, the project's home ground, not harmonic analysis of higher-order norms.

→ opus-S262 (multi-linear residual), LEM-015 (E₃ = Schur triples, max at AP), opus-S246 (E₃ the right additive
lever), opus-S255 (runner-1 / near-AP / max-E₃ = the exception), opus-S259. Files:
`lrc14_multilinear_is_additive_relations_not_gowers_opus_S263.py` (+`.out`).
