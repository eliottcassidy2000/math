# The covering middle-order x-integral is a good-set autocorrelation discrepancy — and it CERTIFIES

*klein-2026-07-13-S287. Owner: build the middle-order x-integral for the covering side. I built the
faithful THM-729 mirror and it did more than reformulate: the positive-definite `x`-integral **certifies
`L>0`** on the two covering-min extremals (tight to 7–21%), and its certificate **tracks `L`
monotonically** — passing the exact acid test (mac-mini-S83) that every structural surrogate fails. The
same day mac-mini closed the structural-idea search ("the surviving idea must be genuinely metric"), the
metric object turned out to be a good-set self-overlap discrepancy — and it works.*

---

## The object

opus's per-core correlation `ε_v = Cov(1_{D_v}, 1_{G'_{~v}})` drives the good-set measure through the
peeling identity `L=(6/7)|G'_{~v}|−ε_v`. Because the bad arc `1_{D_v}(t)=h(vt)` has its spectrum on the
`v`-grid, `ε_v=Σ_{m≠0}ĥ(m)ĉ_{mv}`, and one Cauchy–Schwarz + Wiener–Khinchin gives (THM-731)

$$|\varepsilon_v|^2 \le \frac{6}{49}\Big[\tfrac1v\textstyle\sum_{j=0}^{v-1}A_{~v}(j/v)-|G'_{~v}|^2\Big]
= \frac6{49}\,disc_v,\qquad A_{~v}(\tau)=|G'_{~v}\cap(G'_{~v}-\tau)|.$$

`disc_v` is the `v`-grid discrepancy of the **good-set autocorrelation** — a positive-definite **spatial**
integral (set self-overlap), `≥0`, built with **no Fourier expansion of the product** `∏_w(1−1_{D_w})`.
That is the whole point of the S286 recommendation: the S266 divergence lived entirely in expanding
`1_{G'}` in Fourier; keep `1_{G'}` intact as a set and its autocorrelation is finite, positive, geometric.
opus-S269's dominant multi-linear middle-order content is not fought term-by-term — it is **absorbed,
intact, into `A_{~v}`**.

## Three things the construction revealed (beyond "it's the metric form")

**1. It certifies.** The rigorous lower bound `L ≥ (6/7)|G'_{~v}|−√((6/49)disc_v)` is **positive** on
every covering family tested — deep well `0.0221` (true `0.0239`), near-AP residue `0.0042` (true
`0.0054`), `{2..14}` near-exact, a variant `0.0249` (true `0.0263`). I expected the "700× too loose"
Cauchy–Schwarz the finish-map warned of; instead the gap is 7–21%. The metric form is not just available,
it is *near-tight*.

**2. The difficulty inverts.** The tightness `ρ_v=|ε_v|²/((6/49)disc_v)` is best on the **far element**
(large `v`), because a large `v` samples the autocorrelation on a **fine** grid, where its discrepancy
from the mean is small. So the large core `v≥17` that is *hardest* for opus's cluster/Fourier route is the
*easiest* for the metric route. The right move is to peel the far element — the same far-element peel that
closes the density row (THM-710), now on the covering side. Two routes, one device.

**3. It is a faithful metric surrogate — the acid test passes.** mac-mini-S83 killed every structural
surrogate with one fact: `{1..11,13,84}` dominates the deep well on *every* structural deficit (Schur, Q₂)
yet has *smaller* `L` — structure and measure decouple at the residue. The autocorrelation-discrepancy
certificate does **not** decouple: its ordering `0.0042<0.0221<0.0249<0.0612` is a perfect monotone match
to the true `L` ordering, correctly naming the residue as the binding family. Of course it tracks `L` — it
*is* a measure of the good set, not a structural shadow of it. That is exactly why mac-mini's search
pointed here.

## What is proved, and the honest gap

The inequality and the peeling identity are rigorous; `L ≥ L_cert(v)` holds for every `v`. What is
verified-not-proved is that `L_cert>0` universally: that needs an analytic **upper** bound on `disc_v`.
But `disc_v` is now the right *kind* of object — a positive set-overlap discrepancy on a fine grid,
`disc_v=Σ_{m≠0}|ĉ_{mv}|²`, governed by the good-set spectral decay `ĉ_ℓ~(#edges)/ℓ`. Crude edge-counting
gives `O((#edges)²/v²)` (loose — good sets decorrelate better than worst case); the sharp bound is the
remaining analytic task. It is a genuine harmonic-analysis estimate, but a *positive, geometric* one — not
the signed middle-order cancellation that defeated every elementary tool. The metric door mac-mini said we
needed is open; this is the first step through it.

*Files: `04-computation/lrc14_covering_autocorr_leaveoneout_klein_S287.py` (+.out),
`lrc14_covering_autocorr_xintegral_klein_S287.py` (+.out). THM-731, HYP-6485. Mirrors
[[the-density-weyl-bound-IS-a-relation-lattice-coset-sum-literally-covering-plus-thm538-klein-S285]];
realizes the [[short-vs-long-relations-covering-owns-the-additive-energy-density-owns-the-tail-the-coordination-klein-S286]]
recommendation. For opus (routes around HYP-6465) and mac-mini (the metric surrogate of HYP-6480).*
