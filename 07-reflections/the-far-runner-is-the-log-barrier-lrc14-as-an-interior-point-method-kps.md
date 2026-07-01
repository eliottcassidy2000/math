# The far runner is the log-barrier: LRC(14) as an interior-point method, with the Riesz comb as the central path

*kind-pasteur-2026-06-30. The owner asked to "use w as a log-barrier smoothing inside a Beurling–Selberg / Fourier-positivity signed certificate," and to look back at where the project has needed signed vs. unsigned series. This reflection assembles a signed-vs-unsigned survey, three fresh computations, and the existing canon (THM-515 Riesz product, THM-563 single-far periodicity, HYP-3130 Beurling–Selberg minorant, HYP-3132 multi-far floor) into one organizing frame: **the far speed is a self-concordant barrier, and peeling far runners follows a central path to a bounded-core analytic center.***

> **Concurrent convergence (added at close-out, 2026-06-30).** opus **HYP-3769** independently reached the *same* interior-point frame — on the **covering-min `M`**, not the sector `p0`: the exact split `1/M = (n−1) + 1/rung` reads `ν = n−1` (the danger-arc count) as the self-concordant barrier parameter and the **tight AP as the analytic center** (rung 1), coverings pushed off it (rung ≥ 2). opus pushed first and **owns the covering-min framing**; this is the independent sector-`p0` version (`Plat(C)` = analytic center, `1/w` = barrier coordinate). Two functionals, one frame — mutual corroboration, not priority. opus also works with the log-barrier `F=−∫log m_w` directly (self-concordant), consistent with the correction (below) that the *tent* form is not.

## The one law the survey found

Across every signed/unsigned instance in the repo — the `(−1)^{|T|}` singular series (THM-504, `A_3=∞`), the signed SPEC sum (HYP-3129), the Dedekind sawtooth (THM-563), the Legendre/Gauss `i√7` (the three modes), Pfaffian-vs-permanent (THM-515 §D), Euler's pentagonal signs — the same law holds:

> **The hard, open content is always the signed / odd / imaginary part; positivity methods (Mayer gas, Schur envelope, absolute lattice bounds, SOS) reach only the even / real / positive part, and they diverge or overcount ~100× on the signed part.**

The real side (`ℚ(cos 2π/7)`, equioscillation, cap, binding pairs) is closed. The imaginary side (`ℚ(√−7)`, the Gauss sum `i√7`, the decorrelation floor `R′`, the signed SPEC) is the open frontier. Any certificate must **keep positivity while extracting a signed cancellation** — which is exactly what a Riesz product does (THM-515), and exactly where the uniform Beurling–Selberg minorant fails: HYP-3130's `C^∞` minorant *closes the apex block* `r≤6` but **degrades the coupling ratio**, so it cannot be used uniformly to bridge the far–core coupling.

## The frame: `1/w` is a central-path coordinate

The fix is to stop applying the minorant to the coupling and instead treat the far speed as the barrier. Peel one far runner `w` from a bounded near-cover core `C`. The exact identity (validated, `lrc14_w_as_logbarrier_centralpath_kps.py`, `p0` sector proxy):

```
p0(C ∪ {w})  =  Plat(C)  +  Δ_w ,     Δ_w = (w·Δ_w)/w ,   |Δ_w| ≤ (period-max)/w
```

with `p0(C∪{w}) == Plat(C) + wD(cells_C, w)/w` holding **exactly** for every `w` (match verified for `w = 15…800`, two cores), and the residual `|p0 − Plat|` decaying like `1/w` (consec-8 core: `0.021 → 0.00013 → 0.000035` at `w = 15, 100, 800`) inside the `period-max/w` envelope.

Read this as an **interior-point method**:

- `Plat(C)` is the **analytic center** — the mean-field core measure with the far runner replaced by its mean `6/7`. It is a *bounded, finite* object, positive by LRC≤13 (induction) plus a Beurling–Selberg minorant applied **only to the core** (never the coupling — so HYP-3130's failure mode is sidestepped).
- `1/w` is the **central-path coordinate**. `w = ∞` is the center; finite `w` sits off-center by the barrier residual `Δ_w`.
- `Δ_w` is the **barrier residual**: signed, but **bounded and periodic** (a Dedekind sum, THM-563), decaying `O(1/w)`. It is the signed content — now confined to a *finite periodic object*, not the divergent `|T|≥3` absolute tail. The approach to the center is oscillatory (the sawtooth), enveloped by `1/w`.

  > **CORRECTION (same-day deep study, `lrc14_barrier_residual_dedekind_selfconcordance_kps.py`):** calling `Δ_w` a *self-concordant* residual was aspirational. Rigorously, the sector/tent functional is **piecewise-linear (`F″=0`) → NOT self-concordant**. Self-concordance holds *only for the log-barrier form* (`−log(x)` ratio `=2`; the arc barrier `−log sin πx` ratio `=2|cos| ≤ 2`). The tent and the log barrier compute the **same signed number** `s(h,k)` (Rademacher/`B₂` vs cotangent form). So "interior-point" becomes *literal* only if the functional is recast with the log barrier — which also hands us Dedekind **reciprocity** (a Euclidean `O(log)` descent), the lever the sector period-max lacks. See the companion reflection `the-barrier-residual-is-a-dedekind-sum-...`.

Peeling far runners one at a time = **following the central path inward**; the multi-far floor (HYP-3132, the open `r=2…6` core) is the several-step central path, whose joint residual is the multi-fold Dedekind sum (the Dedekind ladder).

## The far runner *is* the Riesz test measure

The synthesis that ties it to the owner's "Fourier-positivity" ask, and to the attacker=defender reframing: the far runner's safe-factor `1_safe(wτ)` has Fourier support **exactly on the lacunary comb `{k·w : k∈ℤ}`**. A Riesz product on a lacunary/dissociated set is automatically a probability density (`≥0`, Bonami trivial). So:

> The far runner, the log-barrier, and the Riesz test measure are **the same object** on the comb `{kw}`. Pairing the core certificate against this comb reproduces `Plat(C) + Δ_w` exactly; positivity is automatic (`Riesz ≥ 0`), and `w → ∞` makes the comb infinitely lacunary, the Riesz product exact, and the residual `→ 0`.

This lifts "attacker = defender" from the *extremizer* level (the hard sets are Riesz products — confirmed: `{0,1,2,30,31,32,60,61}` is an exact defected Riesz product `\hat 1_B(ξ)\hat 1_B(30ξ) − e(62ξ)`, and its power spectrum peaks at `k/30`) to the *certificate* level (the test measure is a Riesz product on the far comb). The `w·Δ_w` object we've mapped is a two-scale resonance — a fast period-7 carrier (the danger arc) under a slow scale-30 envelope (the Riesz tower) — so the barrier residual's own structure is dictated by the comb.

## The signed step made real: `ℚ(√−7)`

The residual `Δ_w` is where positivity fails, and by the one law it lives in the odd/imaginary part — the Gauss sum `i√7`, because `7 ≡ 3 (mod 4)`. The concrete move (executing the flagged-but-untried HYP-3252/HYP-3254): reorganize `Δ_w` over `ℚ(√−7)`, where `7` ramifies and the imaginary Gauss sum becomes **real**. Then the signed Dedekind residual is a real quantity, sign-analyzable, and the Riesz-on-`{kw}` positive measure carries it — turning the signed barrier step into a positivity statement in `ℚ(√−7)`.

## What is new here vs. assembled

- **Assembled (canon):** the mean-field split `L = (6/7)L_∞(C)+Δ_w` (HYP-3132); `Δ_w` bounded/periodic `O(1/w)` (THM-563); Riesz product as the signed-to-positive tool (THM-515); Beurling–Selberg minorant closes the apex block (HYP-3130); `w`-as-regularizer noted in passing.
- **New (this reflection):** (1) the **interior-point/central-path organizing frame** — `1/w` is the barrier coordinate, `Plat(C)` the analytic center, far-runner peeling = following the central path — which gives a *principled induction structure* for the open multi-far floor; (2) the identification **far runner = log-barrier = Riesz test measure** on the lacunary comb `{kw}`, lifting attacker=defender to the certificate; (3) the prescription to apply the Beurling–Selberg minorant **only to the bounded center**, dodging HYP-3130's coupling degradation; (4) the `ℚ(√−7)` real-ification of the barrier residual.
- **Honest status.** Not a proof. Validated: the exact central-path decomposition and its `O(1/w)` envelope (on the `p0` proxy). Assumed: `Plat(C) > 0` (LRC≤13 induction + core minorant). Open: the multi-step (multi-far) central path — the same `r=2…6` floor the repo already isolates — and whether the `ℚ(√−7)` reorganization actually makes the joint Dedekind residual sign-definite.

## The concrete next probe

Build the truncated Riesz product on `{w, 2w, …, Kw}`, pair it against a bounded-core certificate, and check it (a) reproduces `Plat(C)+Δ_w`, (b) gives a positive lower bound — the far-comb analogue of mac-mini's `1.857 → 1.41` Riesz probe (HYP-2540). If the far-comb Riesz product clears `> 0` where the uniform minorant degraded, the interior-point route is live.

— Related: [[lrc14-thread]], THM-515 / HYP-2540 (Riesz product), THM-563 (single-far Dedekind periodicity), HYP-3130 (Beurling–Selberg minorant, apex block), HYP-3132 (multi-far floor), HYP-3129 (signed SPEC / L2-CS), HYP-3252/3254 (`ℚ(√−7)`), OPEN-Q-108. Companion: `lrc14-the-far-discrepancy-is-a-riesz-dedekind-ladder-dual-to-the-plateau-kps.md`, `lrc14-is-the-lonely-measure-and-the-key-is-a-riesz-product.md`. Artifacts: `04-computation/lrc14_w_as_logbarrier_centralpath_kps.py`, `lrc14_spectral_predictor_and_dft_peaks_kps.py`, `05-knowledge/results/lrc14_{w_as_logbarrier_centralpath,spectral_predictor_and_dft_peaks}_kps.out`.
