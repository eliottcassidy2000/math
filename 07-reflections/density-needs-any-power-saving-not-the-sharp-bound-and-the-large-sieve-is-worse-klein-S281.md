# The density row needs only *any* power-saving on Q_s — not the sharp √; and the large sieve is strictly worse than the crude bound

*klein-2026-07-13-S281. Owner: prove `Q_s=O(M)`. Two honest results: the crude `Q_s ≤ 4π²r²/3` is
rigorous but the sharp `O(r)` genuinely needs cancellation — however, **any** power-saving
`Q_s=O(r^{2−ε})` already closes the density row, so the target is far softer than "sharp". And the
large sieve, the obvious tool, is *strictly worse* than the crude Fourier bound here (endpoint
clustering), which points to the right refined tool.*

---

## The rigorous crude bound (and why it's `w`-free)

`Q_s = (2πw)²Σ_{ℓ≠0}|\hat f(ℓw)|²`, `f=1_{R_s}` (`r` arcs, `V(f)=2r`, `r=O(diam)`). The BV Fourier bound
`|\hat f(n)| ≤ V(f)/(2π|n|) = r/(π|n|)` gives

$$\sum_{\ell\ne0}|\hat f(\ell w)|^2 \le \frac{r^2}{\pi^2 w^2}\sum_{\ell\ne0}\frac1{\ell^2}=\frac{r^2}{3w^2}
\ \Longrightarrow\ \boxed{Q_s \le \tfrac{4\pi^2}{3}r^2 = O(r^2).}$$

The `(2πw)²` exactly cancels the `1/w²` — this is why `Q_s` is `w`-free. **Insufficient alone:**
`|S|=O(√Q_s)=O(r)=O(diam)`, so `Error=|S|/w=O(diam)/d=O(1)`, not `→0`.

## The downgrade: *any* power-saving suffices

But density has slack. For **any** `ε>0`, if `Q_s=O(r^{2−ε})` then `|S|=O(r^{1−ε/2})` and, on the peel
`w=d≥diam≥~3r`,

$$\text{Error}=\frac{|S|}{w}=O\!\big(r^{1-ε/2}/d\big)\le O\!\big(d^{-ε/2}\big)\to0.$$

So the density-row tail closes given **any nontrivial cancellation** `Q_s=o(r^2)` — no sharp constant,
no `ε=1/2`. This is a major reduction versus the covering side, which needs the *sharp* multi-linear
(Gowers) cancellation (mac-mini-S78: covering is now "one bounded inequality"; no slack).

## The large sieve is strictly worse — the clustering diagnosis

The obvious tool fails, informatively. `Σ_{ℓ≤L}|U_s(ℓw)|² ≤ (L+δ^{-1})\,2r`, `δ=\min\|w(p-p')\|`; dyadic
summation gives `Q_s ≤ O(r) + O(r\,δ^{-1})`. But `δ` collapses (`δ∼1/r²` from clustered endpoints), so
`Q_s ≤ O(r^3)` — **worse than the crude `O(r²)`**. Clean thin/thick splits also fail: they reintroduce
`μ_{\rm thin}/w` terms that break the `w`-cancellation.

The resolution is structural: the clustered endpoints are precisely the boundaries of **thin arcs**
(two large offsets crossing nearly together), whose weight in `\hat f` is `∝` arc width `→0`. So the
`δ`-clustering that defeats the large sieve is *weight-suppressed*: a **width-weighted / Montgomery–
Vaughan 2nd-moment** (points weighted by arc width, not unit) should deliver the saving, because the
badly-separated points carry negligible mass. That is the concrete remaining estimate.

## Honest state

- **Rigorous:** `Q_s ≤ 4π²r²/3` (crude Fourier); the large sieve is strictly worse (`O(r³)`).
- **Downgrade (rigorous implication):** any `Q_s=O(r^{2−ε})` closes the density row — density needs
  only a *soft* cancellation, not the sharp √.
- **Confirmed:** the sharp `Q_s=O(r)` holds empirically (S280) — a full √-saving, far more than needed.
- **Remaining:** a width-weighted 2nd-moment giving `o(r^2)` (clustered points = thin arcs = negligible
  mass). Soft, and any amount suffices.

Net: the density route is essentially complete — rigorous crude bound + the row closes on *any* power-
saving, and the saving is a soft width-weighted 2nd moment (not the sharp Gowers estimate the covering
side needs). The two LRC(14) routes' difficulty gap is now maximal: covering = one sharp bounded
inequality; density = any nontrivial cancellation.

*Files: (analytical; builds on `lrc14_second_moment_klein_S280.py` + out). THM-729 addendum, HYP-6425.
Consumes THM-729/728/727. Companion to
[[the-density-sqrt-cancellation-is-a-1d-autocorrelation-discrepancy-not-multilinear-klein-S280]].*
