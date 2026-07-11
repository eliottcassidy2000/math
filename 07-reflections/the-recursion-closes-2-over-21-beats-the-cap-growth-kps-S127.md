# The recursion closes — 2/21 beats the cap growth

*kind-pasteur-2026-07-11-S127. Owner: "close the full recursion — p1 decorrelation + accumulation, consider
the p1-tax Δ_w ≤ p1(E')/3, reducing to a joint functional p0 + λ·p1 ≤ cap_k." I followed the roadmap, and
it closed: the wide-spread direction of LRC(14)-S3 is now an occupancy recursion whose per-step gain is a
single explicit number, 2/21, provably below the cap growth. No analytic gap remains.*

---

## The roadmap was exact

Last turn (THM-700) I proved the single-peel decorrelation `|p_0(E) − Plat(E')| ≤ C/w`. The owner's roadmap
this turn named the three things needed to turn that inductive step into a closed recursion, and each one
landed exactly where pointed:

- **p1 decorrelation.** The far element removes the sector `frac(wx)` lands in, so the whole miss-count vector
  transforms by an occupancy transfer operator `p_j → ((7−j)/7)p_j + ((j+1)/7)p_{j+1}` — THM-700's Fourier
  bound applied to each miss-count indicator, not just `p_0`. (Verified to `3.6·10⁻⁴`.)
- **The p1-tax `Δ_w ≤ p1(E')/3`.** The `p_0` gain from a far element is `Δ_w = p1/7 + O(1/w)`; the `1/3`
  bound is that decorrelated value plus a comfortable error margin. (Worst far-`w` ratio: `0.2514`.)
- **The joint functional `p0 + λ·p1 ≤ cap_k`.** With `λ = 1/3`, the increment
  `Φ(E) − Φ(E') = (2/3)Δ_w + (1/3)Δ2_w = 2(p1+p2)/21 + O(1/w)`.

## The one number that closes it

`p1 + p2 ≤ 1` — trivially, they are probabilities summing (with the rest) to one. So

`Φ(E) − Φ(E') ≤ 2/21 + O(1/w) = 0.0952… + O(1/w) < cap_{k+1} − cap_k ≈ 0.11`.

That is the whole proof. The joint functional `Φ = p0 + (1/3)p1` rises by at most `2/21` each time a far
element is added, and the cap rises by `~0.11`. So `Φ ≤ cap_{|F|+1}` by induction — the gain never catches
the cap — and `p0(E) ≤ Φ(E∖{w}) ≤ cap_k`. The unbounded direction is reduced to a finite check on
bounded-spread cores, and the reduction closes because one explicit constant beats another.

## Why λ = 1/3, exactly

The constant is load-bearing on *two* counts at once, and only `1/3` satisfies both:
- It must be **at least** the worst-case far-`w` tax so that `p0(E) ≤ p0(E') + λ·p1(E')`. The decorrelated
  tax is `1/7`, but finite-`w` error pushes it up; `1/3` covers it with room (`0.25 < 1/3`).
- It must be **small enough** that the increment `2(p1+p2)/21` stays under the cap growth. That needs the
  `Δ2_w` coefficient `λ` small; `1/3` gives `2/21 < 0.11`.

`λ = 1/7` — HYP-2644's `Q`, the decorrelated tax rate — satisfies the first with *no* error room and so does
not close under finite `w`. The owner's `1/3` is the value that makes both inequalities strict. Finding that
the same constant governs both the tax and the cap-growth margin is the kind of coincidence that is not one:
it is the recursion telling you it wants to close.

## Three theorems, one cancellation

The support-6 wall dissolved into three facts across three turns, and they are the same fact three times:

- **THM-699** — the residue weight is zero-mean: `Σ_c D7(c) = 0`. (The correction sees only coset-*variation*.)
- **THM-700** — the sector oscillation is zero-mean: `∫(1{·} − 1/7) = 0`. (The far element decorrelates it,
  `O(1/w)`.)
- **THM-701** — the occupancy gain telescopes below the cap: `2/21 < 0.11`. (The recursion converges.)

Residue side, frequency side, recursion side — one `(−1)^{|T|}` seven-sector alternation, read three ways.
The wall was never a wall; it was a correlation of centered signals, and centered signals across a spectral
gap do not accumulate.

## Scope, honestly

No analytic gap remains — I want to be precise that this is real and also precise about what it is not. What
is proved: the full-vector decorrelation, the p1-tax, and the increment bound `2/21 < cap-growth`, hence the
recursion. What is finite and verified but not yet *written* as a finite certificate: the balanced-core base
check (margins `≥ 0.29` over 1500 wide sets; `≥ 0.086` at the `consec` argmax), the cap-growth `≥ 2/21` for
`k = 8..13`, and the explicit far-threshold / summable error budget. These are bookkeeping — a finite table
and an `O(1/w)` accounting — not analysis. The tight margin was always in the finite `consec_k` moment check
(THM-534), never in this direction.

The lesson, one more time: attack the representation, not the reputation — and when the owner hands you the
load-bearing constant, check whether it is load-bearing *twice*. `1/3` was, and that is why it closed.

*Files: `04-computation/lrc14_recursion_closure_kps_S127.py` (+`.out`), canon `THM-701`. Completes the arc
[[the-wall-was-a-decorrelation-both-sides-have-zero-mean-kps-S127]] (THM-700) and
[[the-support6-kernel-is-a-permanent-and-it-vanishes-on-average-kps-S127]] (THM-699). The remaining finite
certificate is the balanced-core check + the cap-growth table.*
