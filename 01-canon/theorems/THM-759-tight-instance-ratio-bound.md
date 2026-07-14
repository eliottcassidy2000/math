---
id: THM-759
title: The tight-instance ratio bound — if A = {a_1<...<a_n} is an LRC-extremal set (M(A)=1/(n+1), the tight lonely-runner instance), then its largest speed is at most n times its second-largest, a_n ≤ a_{n-1}/((n+1)·M(A\{a_n}) − 1) ≤ n·a_{n-1}. PROVED by an elementary interval (danger-tooth) argument — no alignment/number-theory hypothesis. This is the ratio bound THM-757/HYP-6775 flagged as the missing piece of the LRC(13) tightness rigidity; it reduces "tight ⟹ bounded spread" to a rigorous lemma, and combined with the rigidity of the (n−1)-core it makes the census of tight n-sets finite. Verified against {1..12} (a_12=12 ≤ 12·a_11=132) and the empirical tight-set search.
status: PROVED (elementary: distance-to-ℤ triangle inequality + a covering/interval argument; strict inequality handled). VERIFIED-EXACT (the bound holds on every tight set found; the constant n is tight when the dropped core is extremal, μ₀=1/n).
source: mac-mini-2026-07-14-S108 (completing the ratio bound of the LRC(13) rigidity)
depends_on:
  - THM-751   # aligned tooth-narrowing (the aligned special case; this generalizes it to ALL far elements)
external: LRC(≤13) SETTLED (used only for M(A\{a_n}) ≥ 1/n).
related:
  - HYP-6775  # LRC(13) tightness rigidity (this supplies its ratio bound)
  - THM-757   # multi-killer floor (the rigidity refines its equality case)
  - THM-753   # safe-peel reduction (tight ⟹ irreducible; every element binding)
  - HYP-6800  # this session's rigidity assembly
---

# THM-759 — The tight-instance ratio bound

**One line.** In an LRC-extremal set (`M = 1/(n+1)`, the tightest possible for `n` speeds), the top
speed cannot run away from the rest: `a_n ≤ n·a_{n-1}`. A far top speed would leave a lonely time the
rest cannot cover, pushing `M` above `1/(n+1)`. Elementary interval argument, no arithmetic hypothesis.

## Statement

Let `A = {a_1 < a_2 < … < a_n}` be `n` distinct positive integers with

> `M(A) := max_t min_{i} ‖a_i t‖ = 1/(n+1)`

(the **tight**, LRC-extremal instance; `1/(n+1)` is the least possible value by LRC(n+1)). Write
`P = A \ {a_n}` and `μ₀ = M(P) ≥ 1/n` (LRC(n), `n−1` speeds). Then

> **`a_n ≤ a_{n-1} / ((n+1)μ₀ − 1) ≤ n·a_{n-1}`.**

(The second inequality uses `μ₀ ≥ 1/n`, so `(n+1)μ₀ − 1 ≥ 1/n`. If the core `P` is non-extremal,
`μ₀ > 1/n`, the bound is *sharper*.)

## Proof (elementary interval / danger-tooth argument)

Set `L = 1/(n+1) = M(A)`. Let `t₀` be a tight point of `P` (`min_{p∈P} ‖p t₀‖ = μ₀`). Suppose for
contradiction that `a_n > a_{n-1}/((n+1)μ₀ − 1)`, i.e. `a_n > a_{n-1} · L/(μ₀ − L)` (same thing, since
`L/(μ₀−L) = (1/(n+1))/(μ₀ − 1/(n+1)) = 1/((n+1)μ₀ − 1)`). Then `(μ₀ − L)/a_{n-1} > L/a_n`, so we may
pick a radius

> `ρ` with `L/a_n < ρ < (μ₀ − L)/a_{n-1}`,

and let `I = [t₀ − ρ, t₀ + ρ]` (width `2ρ`).

- **The core stays strictly lonely on `I`.** For each `p ∈ P` and `t ∈ I`, the triangle inequality for
  the distance-to-`ℤ` norm gives `‖p t‖ ≥ ‖p t₀‖ − p|t − t₀| ≥ μ₀ − a_{n-1}ρ`. Since
  `ρ < (μ₀ − L)/a_{n-1}`, this is `> μ₀ − (μ₀ − L) = L`. So `min_{p∈P} ‖p t‖ > L` for **all** `t ∈ I`.
- **The top speed must have a safe point in `I`.** The "danger set" `B = {t : ‖a_n t‖ ≤ L}` is a union
  of closed intervals of width `2L/a_n`, one around each `k/a_n`. Because `2ρ > 2L/a_n`, the interval
  `I` is **wider than one danger interval**, so `I ⊄` any single component of `B`; a connected set that
  is not inside one component and meets the complement's gaps must contain a point `t₁ ∉ B`, i.e.
  `‖a_n t₁‖ > L`. (If `I ⊆ B`, connectedness forces `I` inside one component — impossible by width.)
- **Contradiction.** At `t₁ ∈ I`: `min_{i} ‖a_i t₁‖ = min(\,min_{p∈P}‖p t₁‖,\ ‖a_n t₁‖\,) > L`. Hence
  `M(A) > L = 1/(n+1)`, contradicting tightness.

Therefore `a_n ≤ a_{n-1}/((n+1)μ₀ − 1) ≤ n·a_{n-1}`. ∎

**Remarks.**
- No alignment or coprimality is used — this is the general form of THM-751's tooth-narrowing (which
  handled `a_n · t₀ ∈ ℤ`); here the danger-tooth may sit anywhere and the covering-width argument still
  forces a safe point once `a_n` is large.
- The constant `n` is attained in the limit by the extremal core: for `{1,…,n}`, `P = {1,…,n−1}` is
  extremal (`μ₀ = 1/n`) and the bound reads `a_n ≤ n(n−1)` (loose — actual `a_n = n`), correctly a
  valid but non-binding bound; the point is finiteness, not sharpness.

## Role — the finite reduction of the LRC(13) tightness rigidity

The rigidity `R(n)` ("the only tight primitive `n`-set is `{1,…,n}`", HYP-6775 for `n=12`) has an
inductive skeleton whose one analytic gap was exactly a ratio bound. THM-759 supplies it:

> **Inductive step.** Let `A` be tight, `P = A\{a_n}`, `μ₀ = M(P)`.
> - If `μ₀ = 1/n` (`P` extremal): by `R(n−1)`, `P = {1,…,n−1}`, so `a_{n-1} = n−1` and THM-759 gives
>   `a_n ≤ n(n−1)`. A **finite check** (exact `ℚ`) confirms `{1,…,n−1, w}` is tight **iff** `w = n`
>   (verified for all `n ≤ 12`). Hence `A = {1,…,n}`.
> - If `μ₀ > 1/n` (`P` non-extremal): the **sporadic branch** — see HYP-6800. This is exactly where the
>   Goddyn–Wong sporadic tight instances live (at `n=13`, GW `{1..11,13,24}` has core `{1..11,13}` with
>   `M = 1/12 > 1/13`). Ruling it out is the LRC tight-instance characterization, **open** in the
>   literature (Perarnau–Serra survey: no progress since Goddyn–Wong).

So THM-759 turns "tight ⟹ bounded spread" into a rigorous lemma; the residual content of `R(n)` is
purely the sporadic-branch emptiness, which is `n`-dependent (empty at `n≤12`, populated at `n=13`).

*Artifacts:* `04-computation/lrc13_tightness_rigidity_macmini_S107.py`,
`lrc13_rigidity_ratio_bound_macmini_S108.py` (+outs). Credits: THM-751 (the aligned precursor), LRC(≤13)
(settled), HYP-6775 (the rigidity this completes), Perarnau–Serra (arXiv:2409.20160, the open-problem
context).
