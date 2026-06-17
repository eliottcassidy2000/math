---
id: THM-521
title: The seam-arity law — the Erdős-592 2-adic seam is char-2 (the unique field where the all-ones 2-term Schur sum vanishes); a valuation grading v_p is fully sum-free for a k-term relation ⟺ p=2 ∧ k even, and the diagonal k·a escapes a level ⟺ p∣k; corollary — strong Specker cannot come from any valuation/invariant grading (it would need Larson partial-sum features)
status: PROVED (elementary; both directions below). VERIFIED computationally (`erdos592_seam_arity_law_kps.py`): FULL sum-free ⟺ (p=2 ∧ k even) for k=2..6 over p=2,3,5,7; DIAGONAL k·a escapes ⟺ p∣k for k=2..12 (k=6,12 → {2,3}). Resolves the orphan tangent T778; explains the v₂/v₃ asymmetry of THM-469 algebraically.
source: kind-pasteur-2026-06-16-S4 (a repo-mining session: T778 surfaced in 3 of 4 Explore sweeps as the highest-leverage orphan)
depends_on:
  - THM-469   # sum-free grading dichotomy (v_p level-sets sum-free ⟺ p=2) — this is its k-arity generalization + algebraic reason
  - THM-453   # the gap-grading / triangle = 2-term Schur (gaps add along a chain)
related:
  - HYP-2396  # R(n,2)=2n+1 (the linear wall the p=2 grading runs out at)
  - HYP-2558  # the strong-Specker barrier: invariant gradings can't give a strong witness
  - T778      # RESOLVED here (seam tracks arity — true for the diagonal, refined for the full)
  - OPEN-Q-108
  - reflection: the-seam-is-char-2-and-strong-specker-needs-larson-not-valuations-kps
external: Erdős Problem 592; Specker 1957; Chang 1972; Schur's theorem; Larson's interaction schemes.
---

# THM-521 — the seam-arity law

**Setting.** In the Erdős-592 finite calculus the triangle-free constraint becomes a constraint
on **gaps**: three ordinals `x<y<z` whose pairs are all edges force the gap relation
`g(x,z) = g(x,y) + g(y,z)` — a **2-term Schur relation** `a+b=c`. A translation-invariant
witness grades gaps by a valuation `v_p` (THM-453/469); a level is `L_v(p) = {n≥1 : v_p(n)=v}`,
and for `v=0`, `L_0(p) = {n≥1 : p∤n}` (the "units"). THM-469 found the seam: the `v_p`
level-sets are sum-free ⟺ `p=2`. This theorem gives the **algebraic reason** and the
**arity law**.

## A. The two laws (PROVED)

For the **k-term Schur relation** `a_1+⋯+a_k = b` (all terms in one level):

> **(FULL) `L_0(p)` is `k`-term-sum-free ⟺ `p=2` and `k` is even.**
> **(DIAGONAL) the degenerate relation `k·a = b` is killed by `v_p` (escapes the level) ⟺ `p ∣ k`.**

**Proof of FULL.**
(⟸) `p=2`, `k` even: each `a_i` odd ⟹ `Σa_i ≡ k ≡ 0 (mod 2)` ⟹ `Σa_i` even ⟹ `Σa_i ∉ L_0(2)`,
so no `b∈L_0` closes the relation. Sum-free.
(⟹) Suppose sum-free. If `p` odd: `a_1=⋯=a_k=1∈L_0`. If `p∤k` then `k∈L_0` and `1+⋯+1=k` is
an in-level relation — contradiction. If `p∣k`, take `a_1=2` (a unit, `p` odd), `a_2=⋯=a_k=1`:
`Σ=k+1≡1 (mod p)`, a unit, so `k+1∈L_0` — in-level relation, contradiction. So `p=2`. Then
each `a_i` odd, `Σa_i ≡ k (mod 2)`; if `k` odd, `Σ` is odd `∈L_0` and `1+⋯+1=k∈L_0` is an
in-level relation — contradiction. So `k` even. ∎

**Proof of DIAGONAL.** `v_p(k·a) = v_p(k)+v_p(a)`; for `a∈L_v` this equals `v` iff `v_p(k)=0`
iff `p∤k`. So `k·a` stays in `L_v` iff `p∤k`, i.e. is killed (escapes) iff `p∣k`. ∎

## B. The char-2 reading of the v₂/v₃ asymmetry (the snippet's handoff)

The Erdős-592 triangle is `k=2`. The law says the seam is at the unique prime that is *both*
even-`k`-compatible and divides `k=2` — namely `p=2`, and the two notions **coincide at the
triangle**. The crisp algebraic statement:

> **The seam is `char(𝔽_p)=2`** — the unique field in which the **all-ones 2-term Schur sum
> `1+1` is `0`** (a non-unit). For `p=2`, two units (odds) sum to a non-unit (even) *always*
> (full sum-freeness). For `p` odd, `1+1=2≠0` is a unit, so the minimal Schur triple
> `(1,1,2)` already lives in `L_0` — never sum-free.

This is why `v₂` is the unique sum-free valuation grading (THM-469) and why odd `p` need the
**leading-digit refinement** (THM-469 A2) to recover sum-freeness: only `p=2` has the bare
valuation already kill the even-arity Schur closure.

## C. T778 resolved (the seam tracks arity — with a twist)

The orphan conjecture T778 ("a 3-uniform / higher-arity analogue moves the seam to `p=3`") is:
- **TRUE in the DIAGONAL sense:** the `k`-fold relation's seam is at `p∣k` — `p=3` for `k=3`,
  `p=5` for `k=5` (verified k≤12).
- **REFINED in the FULL sense:** *full* sum-freeness is **exclusive to `p=2` with even `k`**.
  Odd-arity relations admit **no** fully-sum-free valuation grading at any prime. The triangle
  (`k=2`) is special because `k=2` is even AND `2∣2`, so the diagonal and full seams coincide.

So "the seam = the Schur arity" holds for the degenerate/diagonal relation; the strong
invariant-algebra property (full sum-freeness) is a `char-2 ∧ even-arity` phenomenon.

## D. Corollary — the strong-Specker barrier (HYP-2558)

Translation-invariant (valuation-grading) witnesses for ω^n ↛ (ω^n,3)² can use **only the `p=2`
grading** (the unique fully-sum-free one for the even-arity-2 triangle). But the `p=2` grading
runs out at the **linear wall** `t=2n+1` (HYP-2396; the "invariant algebra dies at `t=7`" for
`n=3`, THM-470). Therefore:

> **A *strong* Specker witness (`Q(n,t)` SAT for all `t`) cannot come from any valuation /
> translation-invariant grading.** If one exists, it must use genuinely **non-invariant**
> structure — value-dependent **Larson partial-sum / interaction-scheme** features (the "join
> algebra" of survey §5). This *explains* the `t=7` wall (the invariant algebra is exactly the
> sum-free `p=2` grading, capped at the linear bound) and sharpens the HYP-2396 vs HYP-2363
> fork: either no strong witness exists (HYP-2396), or it lives off the valuation gradings.

## Scope / honesty

PROVED & VERIFIED: both laws (elementary; computational check k≤12). The char-2 reading and the
T778 resolution are exact. The corollary D's barrier ("strong Specker ∉ valuation gradings") is
PROVED *given* that invariant witnesses = valuation gradings (THM-469's framing) and HYP-2396's
linear wall; whether a non-invariant strong witness exists is HYP-2558 (open, = the heart of the
constructive-strong-Specker handoff). This is the algebraic backbone the snippet asked for:
"explain the v₂/v₃ asymmetry" (char-2, §B) and "climb to a t-uniform rung → constructive strong
Specker" (the rung must be non-valuation, §D). Cross-links: THM-469, THM-453, HYP-2396/2558,
T778, the-three-twos (the "2 of pair-composition" = the Schur arity here).
