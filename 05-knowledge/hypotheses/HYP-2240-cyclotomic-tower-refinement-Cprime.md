---
id: HYP-2240
title: The cyclotomic-tower refinement — reaching the optimal 2/(2n-1) via the shell tower (the ramified LRC frontier)
status: OPEN (target); the prime-shell constant-factor bound is PROVED (THM-418)
source: claudebox-2026-06-03-S625
depends_on:
  - THM-418   # prime-shell window lemma (M >= 2/p, p>=2n)  -- the unramified, constant-factor piece
  - THM-415   # the optimal value min M = 2/(2n-1)
related:
  - HYP-2230  # unit-distance / grid-disproof bridge (class field tower = shell tower)
  - THM-407   # doubling shell-transitive iff 2n-1 prime (the unramified case)
  - HYP-2220  # the constructive route to LRC(14): M >= 2/27
---

# HYP-2240 — the cyclotomic-tower refinement

## The gap THM-418 leaves

THM-418 proves, by elementary counting on a prime shell `p ≥ 2n`, the uniform bound `M(S) ≥ 2/p`
(`≈ 1/(2n)`, trivial-strength). The optimal is `2/(2n-1)` (THM-415, `≈ 1/n`). The factor-2 gap is
real and lives at the **extremal shell `m = 2n-1`**, where counting is critically tight. Reaching the
optimum is a question about the **cyclotomic tower** `ℚ(ζ_{2n-1})` and its subfields.

## The conjecture

> **(Tower refinement.)** For every primitive multiple-of-`n` config, a good multiplier giving
> `M ≥ 2/(2n-1)` exists on the shell `m = 2n-1` itself (and, when `2n-1` is composite, on the tower
> of its divisor shells `d | 2n-1` = the cyclotomic subfields `ℚ(ζ_d) ⊆ ℚ(ζ_{2n-1})`).
> Equivalently: the `n-1` runners cannot block every unit `±`-class across the whole tower.

A good multiplier at shell `m` is a unit `u ∈ (ℤ/m)^*` whose `±`-class `{u,-u}` avoids the inverses
of the **unit** runners (non-unit runners — sharing a factor with `m` — are automatically unbanded,
since a non-unit is never `±1`). So the dodge fails at `m` iff the unit-runners cover all `φ(m)/2`
`±`-classes. The **tower descent**: if the runners cover all classes mod `m`, they generically do
*not* cover all classes mod a proper divisor `d | m` (projection loses surjectivity) — so a subfield
shell opens a window. Evidence (`lrc_tower_descent_s625.out`, n=14, `2n-1=27=3³`): over 4000
multiple-of-14 configs the 3-adic subtower `{3,9,27}` supplies the optimal `2/27` for a substantial
share (270 directly on shell 27), the prime shell `29` reliably supplies `2/29`, and only `28/4000`
escape the truncated tower (those use larger shells / Criterion B′). Full coverage to `2/(2n-1)` for
*all* configs is open.

## Why this is the LRC echo of the grid disproof (HYP-2230)

- **Prime shell = unramified prime = constant factor.** THM-418's clean `2/p` is the LRC analogue of
  the disproof's *easy* part: a single split/unramified prime gives a constant-factor supply of
  modulus-1 elements. It is never enough to reach the true exponent.
- **`2n-1` composite = ramified = needs the tower.** When `2n-1` is a prime power (`27=3³`, the first
  even LRC frontier, `n=14`), the single shell is critically tight and the *tower of subfields*
  carries the extra windows — exactly as the disproof needs an infinite **class field tower**
  (Golod–Shafarevich, bounded root discriminant) to pass from constant-factor to `n^{1+ε}`. THM-407
  (doubling shell-transitive iff `2n-1` prime) is the unramified/ramified switch on the LRC side.
- **Bounded root discriminant ↔ uniform-in-`n` control.** The disproof's key is keeping ramification
  bounded as the tower grows. The LRC port: a bound on how the `±`-class coverage can concentrate as
  the shell `2n-1` ramifies — a *uniform* statement that the tower always leaves a free class. That
  uniform control is the missing ingredient between THM-418's `1/(2n)` and the conjectured `1/n`.

## To do
1. Prove the tower-descent lemma: if `n-1` units cover all `±`-classes mod `p^k`, they leave a free
   `±`-class mod `p^{k-1}` (or on a coprime prime shell ≤ a bound). Quantify the "bounded
   discriminant" analogue.
2. Decide the `2n-1` prime case (unramified): does a good multiplier on shell `2n-1` exist for every
   multiple-of-`n` config? (counting gives `≥0`; the extremal attains it — is it always `≥1`?)
3. If both, conclude quantitative `C'(n)`: `M ≥ 2/(2n-1)` for all multiple-of-`n` configs, uniformly
   ⟹ (THM-398) the asymptotic LRC bound `M ≥ 2/(2n-1) > 1/n` on the whole class — the same engine as
   the grid disproof, run finitely.
