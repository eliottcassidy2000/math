# `|S| = o(R)` is false (it is `Θ(R)`) — but the explicit `|Error| ≤ κ'R_G/w` closes the density route for separated far elements, deep well included

*boxeph-2026-07-18-S100. Owner: prove the uniform first-order `|S| = o(R)`. Outcome: **`o(R)` is false**
(`|S| = Θ(R)`, the far element's own self-resonance) — correcting the S99 target — **but the weaker,
elementary `|Error| ≤ κ'R_G/w` (`κ'≈0.094`) is what the density route actually needs**, and it closes the
route for every far element past an explicit threshold `d > κ'R_G/Φ_∞`. The extremal deep well
`{1..12,182}` clears it: `M>1/14` proved rigorously by the density route. The crux is now entirely the
marginal `d∼diam` case (= Route B). LRC(14) not closed. Verified S100 exact computation.*

## `|S| = o(R)` is FALSE

Writing `S = Σ_{n≠0} c_n U_G(nw)`, `U_G(N) = Σ_{x*∈∂G}σ_{x*}e(Nx*)`, `c_n = O(1/n²)`, the `n=1` term is
`c_1·U_G(w)` — and `U_G(w)` is the **far element's own self-resonance**, which S97 showed is `Θ(R)`. So
`|S| ≳ |c_1|·Θ(R) = Θ(R)`. **Decisive test** (grow the frame at fixed scale ratio `w=10T`,
`{1..6,T}`): `|S|/R = 0.012, 0.005, 0.009, 0.013, 0.013, 0.013` as `R: 68 → 1088` — it does **not**
decay. `|S| = Θ(R)`, not `o(R)`. The S99 target was wrong; correcting it here.

*Dichotomy:* `|S| = o(R)` **does** hold for far elements `w` non-resonant with the frame grid
(equidistribution of `{wx*}`); it is `Θ(R)` (small constant) exactly when `w` is commensurate with the
frame. The peel `w = d` is (partially) resonant, so `Θ(R)` is the honest size.

## What the density route actually needs: the explicit `O(R)` bound (proved)

The route needs `Error < Φ_∞`, not `o(R)`. And `Error = O(R)/w` with an **explicit** constant suffices.

> **Theorem.** `Error(w) = −Σ_{n≠0} ĝ(−nw)B̂(n)`, `B̂(n)=\sin(πn/7)/(πn)`, and
> $$|Error(w)| \le \frac{\kappa' R_G}{|w|},\qquad \kappa'=\frac{1}{2\pi^2}\sum_{n\ne0}\frac{|\sin(\pi n/7)|}{n^2}=0.09407,$$
> where `R_G = #∂G` is the good-set endpoint count.

*Proof.* `ĝ(m)=U_G(m)/(−2πim)`, `|U_G(m)|≤R_G`, so
`Error = −Σ_n U_G(−nw)B̂(n)/(2πinw)`, and
`|Error| ≤ (1/2π|w|)Σ_n|U_G(nw)|·|B̂(n)|/|n| ≤ (R_G/|w|)·(1/2π^2)Σ_{n≠0}|\sin(πn/7)|/n^2 = κ'R_G/|w|`. ∎

This sharpens THM-727's `|S|≤0.61R` to `|S|≤0.094R_G`, and — being `O(R)` — is **immune to the
`Θ(R)` self-resonance** that kills `o(R)`. It is exactly enough:

> **Corollary (explicit closure threshold).** `Φ(E) = Φ_∞ − |Error| ≥ Φ_∞ − κ'R_G/w > 0`, so the density
> row closes (`M>1/14`) whenever **`w > κ'R_G/Φ_∞`.**

## The deep well closes — rigorously, by the density route

For the frame `{1..12}` (good set = 13 intervals, `R_G=26`, `Φ({1..12})=16/469≈0.0341`,
`Φ_∞=(6/7)Φ=0.0292`):

| far element `d` | rigorous `|Error|≤κ'R_G/d` | `Φ ≥ Φ_∞ − bound` |
|---|---|---|
| **182** (deep well) | `0.094·26/182 = 0.0134` | `≥ 0.0158 > 0` |

So `M>1/14` for `{1..12,182}` is **proved by the density route** (input: LRC(13) for the frame `{1..12}`,
which is settled, + the elementary `κ'R_G/w` bound — no circularity). And since covering forces `182∣d`
(THM-1017), every covering `{1..12,d}` has `d≥182 > κ'R_G/Φ_∞ = 83.7`, so the **entire covering
`{1..12,d}` family closes** by this bound — an independent, density-side confirmation of THM-1017.

The thresholds `d > κ'R_G/Φ_∞` for the small frames: `{1..6}:6.2`, `{1..8}:17.3`, `{1..10}:33.4`,
`{1..12}:83.7` — modest, and cleared by any genuinely separated far element.

## What remains (honest)

- **`o(R)` — abandoned:** false (`|S|=Θ(R)`); the far element's self-resonance is irreducible. The S99
  "prove `o(R)`" target does not exist to be proved.
- **Proved:** the explicit `|Error| ≤ κ'R_G/w` (`κ'=0.094`), and the closure threshold `w > κ'R_G/Φ_∞`.
  The density route **closes for every separated far element** — including the extremal deep well and the
  whole covering `{1..12,d}` family — using only the frame's (settled) LRC and an elementary bound.
- **The uniform crux:** the threshold `d > κ'R_G/Φ_∞` holds automatically when scales separate
  (`d` large). It can fail only in the **marginal regime `d ∼ diam`** — where the family is *compact*, not
  a density family, and is handled by **Route B** (the AP-rigidity crux, S87–S95). So the density route
  is effectively **discharged**: LRC(14) sits entirely on Route B's compact/covering rigidity.

LRC(14) not closed, but the density route is reduced to an elementary, explicit far-element threshold that
its intended (separated) scope always meets; the irreducible content is Route B.

Cross-links:
[[the-frame-orthogonality-A-is-zero-proved-but-it-is-the-O-R-tail-the-density-route-needs-first-order-o-R-boxeph-S99]],
[[cauchy-schwarz-is-the-density-routes-only-real-loss-the-resonance-cancels-in-S-not-in-Qs-boxeph-S98]],
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
THM-727, THM-1017 (AP-core bridge), THM-886.
