# Cauchy–Schwarz is the density route's only real loss: the resonance cancels in `S` but not in `Q_s`

*boxeph-2026-07-18-S98. Owner: compute `Φ_∞(frame)` and check `√C < Φ_∞`; extend the fixed-frame scaling
law. Outcome: the check **fails for the real `{1..12}` frame** (`√C = 0.145 ≫ Φ_∞ = 0.029`) — but the
**true** two-scale error `→ 0` there, so the density route closes. The gap is **entirely** the
`|S| ≤ √Q_s` Cauchy–Schwarz step: the resonance is phase-aligned in `Q_s` (blows up, `Θ(d²)`) but
phase-cancelling in `S` (the true error vanishes). `Q_s = o(r²)` was the wrong (false) target. LRC(14)
not closed. Verified S98 exact-rational computation.*

## `Φ_∞(frame)` — exact, and the check `√C < Φ_∞`

Good set `G(E) = {t : ‖vt‖ ≥ 1/14 ∀ v∈E}`, `Φ(E) = |G(E)|` (exact rational union-of-arcs), and the
far element equidistributing keeps a `1−2/14 = 6/7` fraction, so `Φ_∞(E') = (6/7)·Φ(E')`:

| frame | `Φ(frame)` | `Φ_∞` | `C_k = Q_s(d)/d²` | `√C` | `√C < Φ_∞`? |
|---|---|---|---|---|---|
| `{1..6}` | `16/35 = 0.457` | `96/245 = 0.392` | 0.137 | 0.370 | **yes** (6% margin) |
| `{1..8}` | 0.266 | 0.228 | 0.038 | 0.194 | yes |
| `{1..10}` | 0.138 | 0.118 | 0.033 | 0.183 | **no** |
| `{1..12}` | 0.034 | **0.029** | 0.021 | 0.145 | **no** |

So the Cauchy–Schwarz proxy `√C` (from S97: `Error ≲ √Q_s/d = √C`) passes for small frames but **fails
for the real LRC(14) frame** `{1..12}`: `√C = 0.145` exceeds `Φ_∞ = 0.029` fivefold. `Φ_∞` collapses with
frame size (the good set shrinks toward the tight AP), while `√C` barely moves. **The `√C` bound cannot
close the real row.**

## But the TRUE two-scale error vanishes — CS is the only loss

Computing the honest quantity `Error(d) = Φ(E'∪{d}) − Φ_∞(E')` directly from good-set measures:

- **`{1..6}`:** `max_d |Error| = 2/35 ≈ 0.057` at the *smallest* far element `d=7`; `Error → 0` as
  `d → ∞` (`d=420` gives **exactly 0**). Far below `Φ_∞ = 0.392`.
- **`{1..12}`:** `Error(d) → 0` for far `d` (`d=420, 840`: `1.1×10⁻⁴`, ratio `0.004` of `Φ_∞`). The only
  tight `d` is `d=13` — which completes the AP `{1..13}` (`M = 1/14` exactly, `Φ = 0`, `Error = −Φ_∞`) —
  and `d=13` is **not a far element**, so the peel `w = d ≥ diam` never selects it.

Since klein's THM-727 gives `|Error| = |S|/w` exactly, `Error → 0` means `|S| = w·|Error| = o(w)`. But
S97 gives `√Q_s = √C·w = Θ(w)`. Therefore

> **`|S| / √Q_s → 0`: the Cauchy–Schwarz step `|S| ≤ √Q_s` is asymptotically infinitely lossy at the
> resonant peel.** The resonance is real in the second moment `Q_s = Σ|U_s|²/ℓ²` (all terms positive,
> `Θ(d²)`) but **cancels** in the signed sum `S = Σ_s Σ_ℓ (−1/2πiℓ)U_s(ℓw)ĝ_s(ℓ)` (phases align in
> `Q_s`, oppose in `S`). Squaring too early throws the cancellation away.

## The mechanism: a fixed frame orthogonality `ν̂ ⟂ ĝ`

Under the S97 scaling law `U_s(ℓw) = w·ν̂_s(ℓ) + o(w)` (`ν̂` the fixed frame section-comb), the leading
part of `S` is `w·A` with

> `A = Σ_s Σ_{ℓ≠0} (−1/2πiℓ)·ν̂_s(ℓ)·ĝ_s(ℓ)`.

`Error → 0` forces `A = 0`: **the fixed-frame comb `ν̂` is orthogonal to the good-set weights `ĝ` in the
`(1/ℓ)`-paired inner product.** This is a *fixed, finite, frame-local identity* — the correct replacement
for the false `Q_s = o(r²)`. The residual (`d=420,840 ≈ 10⁻⁴`) is genuinely lower order, consistent with
`A = 0` exactly and `Error = O(d^{-α})`. The deep well `{1..12,182}` fits: `Φ = 0.0239 > 0`
(`M = 14/183 > 1/14`), `Error(182) = −0.0053`, a small resonant deficit off `Φ_∞ = 0.0292`.

## The redirect (this is the actionable outcome)

- **Stop bounding `Q_s`.** `Q_s = Θ(d²)` is true (S97) and `√Q_s > Φ_∞` for the real frame — the
  `Q_s`-route provably cannot close the row. The "one analytic inequality `Q_s = o(r²)`" was both false
  and the wrong object.
- **Bound `S` directly.** The true error vanishes because the leading amplitude `A = ⟨ν̂, ĝ⟩` cancels — a
  fixed frame-local orthogonality. Proving `A = 0` (and the lower-order residual `= o(1)`) closes the
  density row **without any second-moment cancellation**. This is a change of target, not a bound to grind.
- **Where the cancellation lives:** in the *signs* of `ĝ_s(ℓ)` against the comb `ν̂_s(ℓ)`. The
  derivative gain `sin(πn/7e')` (finish-map §4) that "kills `n=0`" is exactly the structure making
  `⟨ν̂, ĝ⟩` alternate — the same clean bilinear part the covering side enjoys.

## Task 2: the scaling law is universal; `C_k` shrinks with the frame

`Q_s(d) = C_k·d²` holds for every frame `{1..k}` (`C_6,8,10,12 = 0.137, 0.038, 0.033, 0.021`); the
deep-well frame `{1..12}` with far element `182` obeys it. `C_k` is a genuine **invariant of a bounded
speed set** (the `ℓ²`-weighted energy of its section-comb `ν̂`). It shrinks as the frame grows — but so
does `Φ_∞`, faster, which is why the *ratio* (not `C` alone) is what the CS proxy gets wrong.

## Net (honest)

- **The check:** `√C < Φ_∞` is **true for `{1..6}` (barely), false for the real `{1..12}`** — the CS
  proxy is too lossy for large frames.
- **The finding:** the true two-scale error `→ 0` regardless (exact good-set computation); the entire gap
  is the `|S| ≤ √Q_s` Cauchy–Schwarz step. The resonance blows up `Q_s` but cancels in `S`. This
  **retires the S96/S97 "resonance wall" as a Cauchy–Schwarz artifact**, not a real obstruction.
- **The redirect:** the density row closes iff the fixed frame orthogonality `A = ⟨ν̂, ĝ⟩ = 0` — a finite
  frame-local identity, replacing the false `Q_s = o(r²)`.

LRC(14) not closed, but the density route's true remaining piece is now correctly identified: a fixed
orthogonality between the scaling-law comb and the good-set weights, reachable only by bounding `S`
directly (never through `Q_s`).

Cross-links:
[[the-self-similar-resonance-is-a-scaling-law-to-a-fixed-frame-base-not-a-genuine-recursion-boxeph-S97]],
[[soft-weyl-closes-density-off-resonance-but-the-far-peel-is-a-lonely-runner-one-level-up-boxeph-S96]],
THM-727/728/729 (the density error identity + `Q_s`), THM-886 (the resonance law / fixed comb `ν̂`),
[[both-lrc14-routes-bottom-on-the-same-multilinear-cancellation-not-one-mollification-klein-S279]].
