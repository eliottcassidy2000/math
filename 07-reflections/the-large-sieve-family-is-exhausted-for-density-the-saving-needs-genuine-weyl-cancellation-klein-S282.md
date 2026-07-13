# The width-weighting tames the clustering but not the decay — the density saving needs genuine Weyl cancellation, and the large-sieve family is exhausted

*klein-2026-07-13-S282. Owner: prove the width-weighted 2nd-moment bound. Doing so corrects S281 twice
and reaches an honest technical verdict: the width-weighted Montgomery–Vaughan inequality **does** tame
the endpoint clustering (a real finding), but it still yields only `Q_s=O(r²)` — the large-sieve family
cannot supply the `1/w²` decay that the saving needs. That decay is genuine **oscillatory (Weyl)
cancellation**, softer than covering's but of the same category.*

---

## Two corrections to S281

**(1) Thin arcs DOMINATE, they are not negligible.** For a spread cluster (`diam 199`, `w=997`),
`#thick=1` of `58` arcs, and the thin arcs carry the *entire* 2nd moment (`thin·w²=2.21 ≈ full 2.19`,
`thick·w²=0.05`). So `Q_s` is a **width-weighted exponential sum** `f̂(ℓw)=Σ_i w_i\,\mathrm{sinc}(ℓw w_i)
e(-ℓw c_i)` dominated by the many small-width arcs — each contributing `∝ w_i`. S281's "thin arcs are
weight-suppressed / negligible" was backwards: they are the main term.

**(2) The width-weighting DOES tame the clustering.** With midpoints `c_i` and weights `w_i`, the
Montgomery–Vaughan off-diagonal `Σ_i w_i²/δ_i` (`δ_i=` nearest-neighbour of `{w c_i}`) is **bounded and
constant in `r`** (`≈5\times10^{-3}` across `r=3\to81`), while the *unweighted* `Σ_i 1/δ_i` blows up
`∝ r²` (`38 \to 10^5`). So the clustered arcs genuinely carry small weight — the width-weighted large
sieve is not defeated by clustering. This is the correct form of the S281 intuition.

## But it still gives `O(r²)` — the decay is missing

Montgomery–Vaughan (even with the `\mathrm{sinc}` dyadic refinement) gives
`Σ_{ℓ≠0}|f̂(ℓw)|² ≤ O(Σ_i w_i²/δ_i · \mathrm{polylog}) = O(1)` — **bounded in `r`, but `w`-independent**.
Hence `Q_s=(2πw)²·O(1)=O(w²)`, and on the peel `w=d∼r` this is `O(r²)` — no saving. The reason is
structural: the large-sieve family bounds `Σ_ℓ|f̂(ℓw)|²` by an `\ell`-flat mean value, which cannot
see the `1/w²` decay that the sharp `Q_s=O(r)` needs. The `(2πw)²` factor then re-inflates any
`O(1)`-bounded mean value back to `O(w²)`.

## The verdict: the saving is genuine Weyl cancellation

Every soft tool now tried — crude Fourier (`O(r²)`), large sieve (`O(r³)`, worse), Montgomery–Vaughan,
width-weighted MV — caps at `O(r²)`. The `O(r^{2-ε})` saving (which, by S281, is *all* the density row
needs) requires the genuine **oscillatory cancellation** of the width-weighted Weyl sum
`Σ_i w_i e(-ℓw c_i)` under `×w` — i.e. the equidistribution of the arc midpoints, not merely their
spacing. This is softer than the covering side (any `ε>0` vs a sharp constant; 1-linear vs multi-linear
Gowers) but it is the *same category* of estimate: a real equidistribution/Weyl bound, not a mean-value
that the large-sieve toolbox delivers.

## Honest state of the density route (S273–S282)

- **Rigorous:** the whole reduction chain — the `Φ`-transfer (THM-710), the endpoint Fourier reduction
  (THM-727), the 1-D DFT of the derivative (THM-728), the autocorrelation-discrepancy identity
  (THM-729), the crude `Q_s≤4π²r²/3`, and the *any-power-saving-suffices* downgrade.
- **Confirmed:** the sharp `Q_s=O(r)` (S280); the width-weighting tames clustering (this session).
- **The single remaining piece:** an oscillatory (Weyl) cancellation `Σ_i w_i e(-ℓw c_i)=o(\text{trivial})`
  for the arc midpoints under `×w`. Soft (any `ε>0`), 1-linear, but a genuine equidistribution estimate
  — the large-sieve family provably does not reach it. This is a real analytic task, not a further
  one-session reduction.

Net: the density route is reduced as far as soft (mean-value) methods go; the last inch is genuine
equidistribution — the honest bottom of the density side, cleanly separated from and softer than
covering's sharp multi-linear crux (mac-mini-S78).

*Files: `04-computation/lrc14_thin_thick_klein_S282.py`, `lrc14_mv_weighted_klein_S282.py` (+ outs).
HYP-6440. Corrects
[[density-needs-any-power-saving-not-the-sharp-bound-and-the-large-sieve-is-worse-klein-S281]].*
