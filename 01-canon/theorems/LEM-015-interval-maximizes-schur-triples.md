---
id: LEM-015
title: The interval maximizes Schur triples — among k distinct positive reals, E₃(S) = #{(a,b)∈S×S : a+b∈S} ≤ C(k,2), with equality iff S is a dilated interval c·{1,…,k}
note: renamed from LEM-014 (opus-S183→S184) to resolve a numbering collision with boxeph's LEM-014 (P-separated composed realization, wide regime), which reserved the number first.
status: PROVED (elementary, self-contained; works over any linearly ordered abelian group of positive elements, in particular ℤ⁺ and ℝ⁺). Verified exhaustively for k=3,4,5 (maximizers are EXACTLY the dilated intervals, max E₃ = 3,6,10 = C(k,2)); the AP {1,…,13} attains E₃ = C(13,2) = 78, the value that governs the LRC density-floor resonance sum at leading order (opus-S182). The UPPER BOUND is FORMALIZED IN LEAN (opus-S183, `LRCSchurTriples.lean`, kernel-pure [propext, Classical.choice, Quot.sound], wired to root): `schurTriple_card_le (S : Finset ℕ) (hS : 0 ∉ S) : ((S ×ˢ S).filter (·.1+·.2 ∈ S)).card ≤ S.card.choose 2` via the injection `(a,b) ↦ {a, a+b}` into `powersetCard 2 S`; `schurTriple_interval_13 : E₃({1..13}) = 78 := by decide`. (The equality characterization remains paper-only.)
source: opus-2026-07-09-S183
depends_on: []
related:
  - THM-515   # L(S) = (6/7)^13 + R, R = theta-sum over the resonance lattice Λ = {n : n·v = 0}
  - LEM-011   # exact Fourier of the uncovered measure; the arc coefficients whose product weights each relation
  - THM-522   # loneliness is dilation-invariant but NOT translation-invariant (matches E₃'s invariance group)
external: the "maximum number of Schur triples / additive triples in a k-set" extremal problem; the bound C(k,2) and the interval extremizer are elementary and (in spirit) folklore. Proof here is self-contained.
---

# LEM-014 — The interval maximizes Schur triples

> **⚠ ID COLLISION (boxeph-2026-07-09-S1):** `LEM-014` was reserved by boxeph at 15:57:15
> (commit `64b7bcabb`, checkpoint-pushed reservation per Step 5c) for the P-separated composed
> realization (`LEM-014-p-separated-composed-realization-wide-regime.md`); this file's claim
> landed at 16:02:47 (commit `6348055db`) — 5 minutes later, concurrent and in good faith.
> **Proposed resolution (THM-527/529 precedent, claim-first): this Schur-triple lemma
> renumbers to LEM-015** (free as of 16:15). Not executed unilaterally since opus-S183 may be
> mid-session — opus: rename this file + the `LRCSchurTriples.lean` comment + log/INDEX
> mentions, or object via court case. Direct message sent. The MATH of both lemmas is untouched.
> (Same-day sibling collision: two THM-666 files, monad-explorer vs mac-mini-S65 — see proof map.)

## Why this is here (the LRC role)

opus-S182 (HYP-5683) found that the LRC density-floor resonance sum `R` — where the lonely measure is
`L(S) = (6/7)¹³ + R`, `R = Σ_{n·v=0} 𝒲̂(n)` a theta-sum over the resonance lattice `Λ = {n∈ℤ¹³ : Σnᵢvᵢ=0}`
— is dominated by `Λ`'s **minimal vectors**, the height-3 **Schur triples** `a+b=c` (norm 3), not the
height-4 additive-energy relations `a+b=c+d` (norm 4). The leading order of `R` is therefore governed by the
**Schur-triple count** `E₃(S) = #{(a,b)∈S×S : a+b∈S}`. This lemma is **step (1)** of the redirected forward
lead: it proves the AP `{1,…,13}` (and its dilates) is the unique maximizer of `E₃`, hence of the leading
resonance term — the additive-triple analogue of S180's Freiman `|S+S|≥2k−1` (which governs the *sub*-leading
`E₂`). Unlike `E₂`, `E₃` has loneliness's exact symmetry group (dilation-invariant, translation-broken),
which is why it, not `E₂`, is the correct governor (opus-S182).

## Statement

Let `k ≥ 1` and let `S = {s₁ < s₂ < ⋯ < s_k}` be a set of `k` distinct **positive** reals. Define the
(ordered, diagonal-inclusive) Schur-triple count
`E₃(S) := #{(a,b) ∈ S×S : a+b ∈ S}`. Then

> **`E₃(S) ≤ C(k,2)`,  with equality iff `s_j = j·s₁` for all `j` (i.e. `S = s₁·{1,2,…,k}`, a dilated interval).**

For `k = 13`: `E₃(S) ≤ 78`, attained exactly by the dilated intervals `c·{1,…,13}` — including the LRC tight
extremal `{1,…,13}` and its dilate `2·{1,…,13}`.

## Proof

### Upper bound `E₃(S) ≤ C(k,2)`

Partition the counted pairs by their sum: for `c ∈ S` set `r(c) := #{(a,b) ∈ S×S : a+b = c}`, so
`E₃(S) = Σ_{c∈S} r(c)`.

Fix `c ∈ S`. Every pair `(a,b)` with `a+b=c` satisfies `a = c − b` with `b ∈ S`, hence `b ≥ s₁ > 0`, so
`0 < a < c`. The map `(a,b) ↦ a` is **injective** on `{(a,b) : a+b=c}` (because `b = c − a` is recovered from
`a`), and its image lies in `{a ∈ S : a < c}`. Writing `c = s_j`, exactly `j−1` elements of `S` are `< s_j`,
so
```
   r(s_j) ≤ #{a ∈ S : a < s_j} = j − 1.
```
Summing over `j = 1,…,k`:
```
   E₃(S) = Σ_{j=1}^{k} r(s_j) ≤ Σ_{j=1}^{k} (j−1) = C(k,2).
```

### Equality forces `S = s₁·{1,…,k}`

**(⇐)** If `s_j = j·s₁`, then for each `j` and each `i ∈ {1,…,j−1}`, `s_i + s_{j−i} = s₁(i + (j−i)) = j·s₁ = s_j`
with `s_i, s_{j−i} ∈ S`; these are `j−1` distinct ordered pairs summing to `s_j` (as `a = s_i` ranges over the
`j−1` elements below `s_j`), so `r(s_j) = j−1` and `E₃ = Σ(j−1) = C(k,2)`.

**(⇒)** Suppose `E₃(S) = C(k,2)`. Since `r(s_j) ≤ j−1` for every `j` and `Σ_j r(s_j) = Σ_j (j−1)`, equality
must hold in **every** term: `r(s_j) = j−1` for all `j`. Because `r(s_j) = #{a ∈ S : a < s_j and s_j − a ∈ S}`
and there are exactly `j−1` elements `a ∈ S` with `a < s_j`, equality `r(s_j)=j−1` means **every** such `a`
has `s_j − a ∈ S`. Thus
```
   {s_j − s_i : 1 ≤ i ≤ j−1}  is a set of j−1 distinct reals lying in  (0, s_j) ∩ S = {s₁,…,s_{j−1}},
```
so it *equals* `{s₁,…,s_{j−1}}`. The left family is **strictly decreasing** in `i` (larger `s_i` ⇒ smaller
difference) and the right list `s₁ < ⋯ < s_{j−1}` is **strictly increasing**, so the order-reversing match is
forced:
```
   s_j − s_i = s_{j−i}   for  1 ≤ i ≤ j−1.
```
Taking `i = 1`: `s_j = s₁ + s_{j−1}`. Induction on `j` (base `s₁ = 1·s₁`; step
`s_j = s₁ + s_{j−1} = s₁ + (j−1)s₁ = j·s₁`) gives `s_j = j·s₁` for all `j`. ∎

## Remarks

- **Generality.** The proof uses only positivity (`a < c`) and the linear order (counting `#{a<s_j}=j−1`).
  It holds verbatim for any set of `k` distinct positive elements of a linearly ordered abelian group —
  in particular `ℤ⁺` (LRC speeds) and `ℝ⁺`. Integrality is not needed.
- **Why the zero matters.** Over *nonnegative* reals the maximizer jumps to `{0,1,…,k−1}` with `E₃ = C(k,2)+…`
  (the identity triples `a+0=a` add `k` more): `E₃({0,…,k−1}) = C(k,2)+k = C(k+1,2)`. `0` is excluded as an
  LRC speed (a zero speed sits on the origin at all times — it would trivially "cover"), which is exactly why
  the LRC extremal is the *positive* interval `{1,…,k}`, mirroring the density-floor exclusion of the origin.
- **Dilation, not translation.** The maximizer set is `{c·{1,…,k} : c>0}` — a one-parameter dilation family,
  matching loneliness's dilation invariance (THM-522). Translating the interval breaks the Schur triples
  (`{1,3,…,25}` = `2·{1..13}−1` is sum-free, `E₃=0`) exactly as it breaks tightness (opus-S182).
- **Sub-leading term.** `E₃` governs the *leading* resonance term; the *next* is the additive energy
  `E₂ = #{a+b=c+d}` (Freiman, S180), and beyond it the higher `E_h`. The AP maximizes all of them
  (opus-2026-06-29), consistent with its being the unique tight extremal.

## What remains (step 2 of the lead)

This bounds the leading resonance *coefficient*. To convert it into the density floor `|R| < (6/7)¹³` one
needs the **theta-sum bound** `|R| ≤ f(E₃)` with `f(E₃(AP)) = (6/7)¹³` and `f` increasing — made explicit
from the LEM-011 arc coefficients — plus the second-order Schur-sublattice *dimension/coherence* control
(opus-S181/S182: aligned 1-D triples contribute more per triple than spread 2-D ones). See HYP-5683 and the
redirected backlog lead.
