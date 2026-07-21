---
source: klein-2026-07-21-S400
status: proved structural laws + a reduction (THM-1950); the residual strong base is open
tags: [tournament, skew-determinant, disc, Redei-H, SCC-composition, velocity-addition, SL2, WOWII, HYP-8636]
---

# The skew-determinant composes by velocity addition, and that reduces H ≥ disc to the strong core

Owner: synthesise the repo's progress; work high-leverage open math to the most fundamental level;
make small improvements toward proofs; pull often. Target: death-star-S78's HYP-8636,
`H(T) ≥ disc(T)` (equality iff transitive) — the Rédei Hamiltonian-path count (#P-hard) dominated by
the poly-time skew-determinant `disc = |det(I+K)|/2^{n−1}`, `K=A−Aᵀ`. It is the cleanest open
WOWII inequality: it ties the deepest *combinatorial* invariant (H) to the deepest *algebraic* one
(a determinant).

## What I found (all machine-verified; THM-1950)

The natural first move is THM-1860's: reduce to strongly connected via the SCC decomposition, since
`H` is multiplicative (`H(T)=∏H(Cᵢ)`, a Hamiltonian path crosses the linear SCC order). **But disc
is NOT SCC-multiplicative** — and the way it fails is the whole story. Define the *total
inverse-response* `s(T) = 𝟙ᵀ(I+K)⁻¹𝟙`. (Because `K` is skew, `s = ‖(I+K)⁻¹𝟙‖² ≥ 0`, and `≤ n`.)
Two exact composition laws, both proved (Schur complement + block solve), for the ordered sum
`T₁ ⇒ T₂`:

```
disc(T₁ ⇒ T₂) = disc(T₁)·disc(T₂)·(1 + s₁s₂)/2          [super-multiplicative: the factor ≥ 1]
   s(T₁ ⇒ T₂) = (s₁ + s₂)/(1 + s₁s₂)                     [relativistic velocity addition!]
```

The `s`-law is Einstein velocity addition — `s = tanh(rapidity)`, and rapidities *add* under `⇒`.
So the SCC composition acts on the invariant `s` by **Möbius / SL₂ transformations**, and disc's
composition factor `(1+s₁s₂)/2` is exactly the *denominator* of that Möbius map (halved). The
skew-determinant is not multiplicative because it lives on an `SL₂`-torsor, not a `Gₘ`-torsor.

This is the resonance worth flagging: the fleet has been building the a/b = `⟨x+1, x/2⟩ = BS(1,2)`
monoid picture (THM-1885) and the char_A `= SL₂` binary-form frame (THM-1810/1875/1925). Here the
*same `SL₂`* surfaces from a completely different door — the skew-determinant's behaviour under
strong-component composition. `disc` and `s` are an `SL₂` pair; `char_A` multiplicativity (THM-1925)
is the `A`-side, `disc` velocity-addition is the `K`-side, and they are the two faces of the one
`SL₂` acting on the SCC condensation.

## The reduction it buys

Iterating gives `disc(T) = (∏disc(Cᵢ))·[∏(1+sᵢ)+∏(1−sᵢ)]/2^r`. With the invariant
`P(T) = max(1, s(T))·disc(T)`, peeling the top strong component and one elementary inequality
`max(1,x)max(1,y) ≥ max(1+xy, x+y)/2` close the induction:

> **`H(T) ≥ disc(T)` for all `T`, provided the strong base `H(C) ≥ max(1, s(C))·disc(C)`.**

This is the exact structural twin of THM-1860 (`c₃ ≤ H` reduced to strong + the sum-≤-product
kernel): a proved reduction, a proved algebraic kernel, an honest open strongly-connected base.
Verified in invariant form `H ≥ P(T)` over *all* `2^{C(n,2)}` tournaments `n ≤ 6`, base sampled
`n=7` (6387/6387). Equality iff transitive; the base is tight only at `C₃`.

## The instructive correction (small-case discipline)

My first attempt used `s ≥ 1 for strong`, true `n ≤ 6` (min s = 3,3,2,**1**) — and **false at
`n ≥ 7`** (min strong `s` = 0.667, 0.556, …). The naive kernel `∏sᵢ ≥ Φ` dies with it. The repair is
the two-sided `max(1,·)` kernel, proved via the identity `aᵢ + |bᵢ| = 2` (with
`aᵢ = (1+sᵢ)/max(sᵢ,1)`), which collapses to `2·Σ_{k≡r (2)} eₖ(c) ≤ 2·2^{r−1}` for `cᵢ∈[0,1]`. The
`min s = 1` at `n=6` was a marginal-threshold artifact — the same trap as the H-spectrum `{7,21,35,39}
→ {7,21}` (death-star-S70) and `c₃≤H`'s siblings breaking at `n=10`. The `SL₂`/velocity picture
*predicts* it: composition drives `s` toward the `tanh` band `(−1,1)`, so `s<1` is generic and only
the small-`n` strong atoms sit above 1.

## Residual and next steps

- **Open base (the content):** `H(C) ≥ max(1,s(C))·disc(C)` for strong `C`. Room grows
  (`H/(max(1,s)disc)` = 1, 1.67, 1.875, 3.75, 4.22 for `n=3..7`); `C₃` tight. Attack: a
  Hamiltonian-path injection dominating `s·disc`, or an eigenvalue-product bound on the strong
  spectrum (Perron + isotropic pairs, THM-1858; fixed energy `Σμ_j²=C(n,2)`).
- **Structural leads:** disc/`s` as an explicit `SL₂` representation of the condensation monoid;
  the `s = tanh(rapidity)` additive coordinate as a new tournament invariant (rapidity is *additive*
  over `⇒`, unlike everything else); does the rapidity refine the H-spectrum?

Files: `04-computation/{h_ge_disc_reduction_klein_S400, disc_composition_law_klein_S400,
h_ge_disc_reduction_to_strong_klein_S400}.py` (+ outs). THM-1950; answers/reduces HYP-8636.
Cross-links: THM-1860 (template), THM-1925/1810/1885 (the `SL₂`/a-b frame), THM-1858 (spectrum),
THM-474 (disc/skew-det), death-star-S78 (HYP-8636).
