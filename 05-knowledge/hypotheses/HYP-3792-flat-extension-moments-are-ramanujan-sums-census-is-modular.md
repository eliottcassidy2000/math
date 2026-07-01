# HYP-3792 — The flat-extension moments are Ramanujan sums; the 11-census is organized by (Z/N)*; ζ is the unit-restricted Dirichlet denominator

- **Status:** CONFIRMED (all identities verified exactly / high-precision numerically)
- **Source:** kind-pasteur-2026-07-01-S7
- **Prompt:** improve understanding of the census and the Carathéodory–Toeplitz flat-extension; search prior work; consider the Riemann zeta function and modular multiplication specifically.
- **Depends on / relates:** HYP-3789 (flat-extension atoms = units), HYP-3773 (margin = Dedekind sum → −1/12 = ζ(−1)), THM-503 (no Euler product for L — refined here), THM-501/515 (singular series = lonely measure), THM-523 ({AP, Goddyn–Wong} dichotomy), THM-560/522 (census template, scale-invariance), HYP-3762 (three-gap/CF geometry), OPEN-Q-108.

## Claim (six verified facts)

1. **AP maximizers = units.** The covering-min maximizers of `{1..N−1}` are exactly `{a/N : a ∈ (Z/N)*}` (N=12,13,14 → φ=4,12,6 atoms). HYP-3789's "flat-extension atoms" are literally the units.
2. **Moment matrix = Ramanujan-sum Toeplitz.** For the atomic measure `μ = Σ_{a∈(Z/N)*} δ_{a/N}`, `μ̂(k) = c_N(k)` (Ramanujan sum). The Toeplitz matrix `[c_N(j−ℓ)]` has rank `φ(N)` exactly (flat at window `φ, φ+2, φ+4`); moments match the multiplicative formula `c_N(k)=μ(N/g)φ(N)/φ(N/g)`, `g=gcd(N,k)`.
3. **ζ enters as the denominator.** `Σ_{N≥1} c_N(k) N^{−s} = σ_{1−s}(k)/ζ(s)` (verified s=2,3; k=1,6,12). This is the unit-restricted local density THM-503(4) flagged as untested — it refines (not contradicts) THM-503: the archimedean measure value L has no Euler product, but the atomic-skeleton moments are multiplicative Ramanujan sums carrying 1/ζ.
4. **Census organized by (Z/N)*.** The 78 cores `{1..13}\{i,j}` sort by `meas(L_C@1/14)`, each labelled by its maximizers' `(Z/N)*`. Minimum = full-orbit `(Z/10)*` core `{1..13}\{6,10}` (`313/9702 = 0.03226`); runner-up = partial-orbit `(Z/19)` two-speed clash (`0.03238`) — the 11-speed echo of THM-523's `{AP, Goddyn–Wong}`.
5. **Collar width = modular inverse.** Binding runner at atom `a/N` has speed `a^{−1} mod N` (and `N−a^{−1}`); `meas(L_C)=Σ_{a∈(Z/N)*} width(a)`, a Dedekind-shaped units-sum (kin to HYP-3773's `−1/12 = ζ(−1)`).
6. **OPUC/Verblunsky thermometer** (fills a repo gap): all `|α_k|<1` (PD, valid measure); the extremizer is near-atomic (`|α_{5,9,11}|≈0.92–0.98`, Szegő product `1e−5`) vs the fatter `{1..11}` (`max 0.65`, product `0.19`). Minimize `meas` ⟺ maximize atomicity ⟺ push Verblunsky coefficients to the unit circle.

## Honesty

Verified identities + an organizing principle (WHY the extremizer is `(Z/10)*`), NOT a proof of `inf meas ≥ 1/36`. The census is a search over one family (2-drops of `{1..13}`); modular labels explain the minimizer and split the near-tight locus into `{full unit orbit, partial two-clash}` but do not yet bound all primitive 11-cores. Next: prove the two-family dichotomy exhaustive (à la THM-523 finiteness, band 1/14, 11 speeds; THM-522 bounds the search to scale-1), then evaluate both families in closed form and check `≥ 1/36`.

## Artifacts

- Reflection: `07-reflections/the-flat-extension-moments-are-ramanujan-sums-the-census-is-modular-kps.md`
- Scripts: `04-computation/lrc14_census_toeplitz_ramanujan_kps.py`, `04-computation/lrc14_census_modular_verblunsky_kps.py`
- Outputs: `05-knowledge/results/lrc14_census_toeplitz_ramanujan_kps.out`, `05-knowledge/results/lrc14_census_modular_verblunsky_kps.out`
