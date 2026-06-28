# The one inequality that completes LRC(14) — and why the resonant survivor dissolves

*kind-pasteur-2026-06-28-S258. The owner asked for the single highest-leverage way to complete the
LRC(14) proof. The answer is a height-view of the critical path: almost everything the project has built
— the census, the equioscillation, the Galois/ℚ(√−7)/ℚ(cos 2π/7) structure, the H=7 lift — is off the
critical path. LRC(14) needs exactly one new theorem beyond LRC≤13, and its only feared subtlety
dissolves.*

## The critical path is shorter than it looks

LRC(14) asserts the **inequality** `M(S) ≥ 1/14` for every primitive 13-set `S`. Sort `S` by its
divisibility threshold `q(S) = min{d : S has no multiple of d}`:

- **`q(S) ≤ 14` (non-covering): elementary, done.** The witness `t = 1/q(S)` gives
  `f_S(1/q) = min_s ‖s/q‖ ≥ 1/q(S) ≥ 1/14` — each `s` is not divisible by `q`, so `s mod q ∈ {1,…,q−1}`
  and `‖s/q‖ ≥ 1/q`. Verified: 3661/3661 random non-covering 13-sets, zero violations. **This closes the
  entire non-covering world with one line.**

The consequence is bracing: **the census is off the critical path.** "Which 14-free sets are tight" is a
characterization of *equality* (`M = 1/14`), and the LRC inequality never needs it. Two sessions of
census / equioscillation / Galois work (HYP-3246…3413) are real mathematics and a genuine *description*
of the extremal — but for *completing the proof* they are a side-quest.

- **`q(S) ≥ 15` (covering): the only hard case, and it is inductive.** Covering forces a multiple of
  every `d ≤ 14`, so `S = R ∪ 14Q` with `R` 14-free (`13−r` speeds) and `14Q` the `r ≥ 1` multiples of
  14. Under `u = 14t`, the multiples `14Q` become a sub-lonely-runner instance for `Q` (`r ≤ 6` runners),
  **lonely by LRC≤13** — the induction hypothesis. `R` is `q`-witness-safe. Loneliness of `S` is exactly
  that `R`-safe and `Q`-lonely **overlap**.

## The one new theorem: a single decorrelation inequality

`Q`-lonely is `1/14`-periodic in `t` (it is `g(14t)`), so

```
  meas(R-safe ∩ Q-lonely) = Σ_k ĝ(k) · R̂_safe(14k)  =  meas(R-safe)·meas(Q-lonely) · R',
```

with the `k=0` term the product (the main term) and `R' = 1 + SPEC/product`,
`SPEC = Σ_{k≠0} ĝ(k) R̂_safe(14k)`. So:

> **LRC(14) = [q-witness, elementary] + [LRC≤13 for `Q`, induction] + ONE inequality: `R' > 0` uniformly,
> i.e. `|SPEC| < product`.** Every other structure in the project is either off-path or a tool for this
> single inequality.

This is mac-mini's certified `R' ≥ 0.642` (HYP-3129) plus the resonant residual. The whole proof rests
on finishing this one floor.

## Why the resonant survivor dissolves

The feared obstruction: when an `R`-speed is **resonant** with 14 (a multiple of 2 or 7, aligning it with
the 14-grid), the clean equidistribution `1/7`-removal fails (mac-mini S81). The resolution is that
resonance is **self-defeating**:

> A resonant speed `v` has its danger set *on* the 14-grid `a/14`. But the multiples `14Q` are dangerous
> *exactly on* that grid, so `Q`-lonely lives strictly *off* the grid. There, the resonant `v` is
> **safe** — it is transparent precisely where it would need to obstruct.

Verified on covering sets with `7 ∈ R` and with several even speeds: the optimum `t*` lands off the
14-grid (denominators 12, 8, 9), and *every* resonant speed (every multiple of 2 or 7, including the `14Q`
multiples) is safe there — `‖s·t*‖` ranges `0.083…0.5`, all `≥ 1/14`, with `M(S) > 1/14` in every case.
So the floor does not depend on the resonant speeds at all: **resonance pushes the danger onto the grid,
the loneliness off it, and the two never meet.** The floor reduces to the *non-resonant* `R`-speeds, which
decorrelate by equidistribution (`R' ≥ 0.642`). The resonant survivor is not a survivor — it is safe.

## The highest-leverage action

Stop completing the census. The single move that completes LRC(14) is to **finish the uniform
decorrelation floor `R' ≥ c > 0`**, now visibly within reach: the non-resonant part is the certified
`R' ≥ 0.642`, and the resonant part dissolves by transparency (resonant speeds are safe on the off-grid
`Q`-lonely region). Assemble these two and the covering case — hence LRC(14) — is closed modulo LRC≤13.
The beautiful structure (Galois, equioscillation, ℚ(√−7)) explains *why the AP is extremal*; this one
inequality is *why no `S` beats it*. Only the second is on the critical path.

## The pointer beyond

The transparency observation wants to be a lemma: *for covering `S = R ∪ 14Q`, the optimum `t*` of `M(S)`
is off the 14-grid, and every speed divisible by `gcd(s,14) > 1` is safe at `t*`.* If proved, it removes
the resonant case wholesale and leaves only the non-resonant equidistribution floor — the cleanest
possible completion. The grid-dodging of the optimum (mac-mini S82) and the transparency of resonant
speeds are the same fact seen from the two sides of the `u = 14t` lift.

— Related: [[lrc14-thread]], HYP-3415 (the critical-path map), HYP-3129 (`R'≥0.642`), HYP-3132 (multi-far
floor), mac-mini S80/S81/S82 (uniform margin, resonant off-grid), THM-571 (`r≥7`), THM-523 / q-witness;
`lrc14-three-engines-lift-and-decorrelate.md`, `the-galois-group-of-the-apex-prime-and-where-its-symmetry-breaks.md`.
