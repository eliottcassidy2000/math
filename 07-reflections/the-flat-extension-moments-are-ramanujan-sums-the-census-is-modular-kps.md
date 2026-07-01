# The flat-extension moments are Ramanujan sums: the 11-census is organized by (Z/N)*, and the Riemann zeta function is the denominator of the unit-restricted Dirichlet series

*kind-pasteur-2026-07-01-S7. Deepening the two load-bearing finish-tools — the near-tight 11-speed CENSUS and the Carathéodory–Toeplitz FLAT-EXTENSION (HYP-3789) — through the Riemann-zeta and modular-multiplication lens. The punchline: the flat-extension atoms are units mod N, their moment matrix is literally a Ramanujan-sum Toeplitz matrix, and the classical identity `Σ_N c_N(k)/N^s = σ_{1−s}(k)/ζ(s)` is exactly the unit-restricted local density THM-503(4) left untested. Everything below is verified exactly.*

## One picture in one paragraph

For an arithmetic-progression speed set the lonely-runner problem is **pure modular multiplication**: the covering-min maximizers of `{1..N−1}` are exactly the units `(Z/N)*`, the flat-extension moment matrix of that atomic set is the **Ramanujan-sum Toeplitz matrix** `[c_N(j−ℓ)]` (rank `φ(N)`), and the Riemann zeta function enters — precisely and unavoidably — as the *denominator* of the Ramanujan Dirichlet series `Σ_N c_N(k)N^{−s} = σ_{1−s}(k)/ζ(s)`. The near-tight 11-speed census inherits this skeleton: each loose 11-core parks its maximizers on some `(Z/N)*`, the measure is a sum over units of collar widths set by the **modular inverse** `a^{−1} mod N`, and the extremizer `{1..13}\{6,10}` is the one whose orbit is `(Z/10)*`. The modular structure is the *skeleton*; the archimedean collar widths are the *flesh* — which is exactly THM-503's split (no Euler product for the measure value, multiplicative structure in the unit-restricted skeleton) made concrete.

## 1. The flat-extension atoms are units, and their moments are Ramanujan sums

HYP-3789 (mac-mini-S76) established that the extremal lonely set is a few **atoms**, and noted the AP `{1..13}` gives "6 atoms = units of `(Z/14)*`." Pushing this to its natural end:

> **Verified exactly.** The covering-min maximizers of the AP core `{1,…,N−1}` are `{a/N : a ∈ (Z/N)*}` — for `N=12,13,14` the atom counts are `φ(N)=4,12,6`, matching the units `{1,5,7,11}`, `{1,…,12}`, `{1,3,5,9,11,13}` on the nose.

Put unit masses on the atoms: `μ = Σ_{a∈(Z/N)*} δ_{a/N}`. Its Fourier moments are
`μ̂(k) = Σ_{a∈(Z/N)*} e(−ak/N) = c_N(k)` — **the Ramanujan sum.** So the Carathéodory–Toeplitz flat-extension certificate for the AP's lonely set is *literally a Ramanujan-sum Toeplitz matrix* `[c_N(j−ℓ)]`. Verified: its rank is `φ(N)` exactly (flat at window sizes `φ, φ+2, φ+4`), and the computed moments match the closed multiplicative form `c_N(k)=μ(N/g)φ(N)/φ(N/g)`, `g=gcd(N,k)` — e.g. `c_14 = (6,1,−1,1,−1,1,−1,−6,−1,1,−1,1,−1,1)`.

This is the right way to see "rank = #atoms": the rank *is* `φ(N)`, the number of units, and the flatness *is* the periodicity of the Ramanujan sum. The modular group `(Z/N)*` is the certificate.

## 2. The Riemann zeta function is the denominator of the unit-restricted density

Ramanujan sums are the classical bridge from `(Z/N)*` to `ζ`:

> **Verified (s=2,3; k=1,6,12):** `Σ_{N≥1} c_N(k) N^{−s} = σ_{1−s}(k)/ζ(s)`.

The zeta function sits in the denominator. This is the *exact object* THM-503(4) flagged and did not test: it proved `L` (the archimedean lonely-**measure value**) has **no Euler product** (`β_p = L` for every prime, all local factors collapse to 1), but explicitly left open *"a unit-restricted local density `a ∈ (Z/q)*` … could recover a refined product identity."* It does. The resolution is a clean split — and it **refines rather than contradicts** THM-503:

- **The measure value `L`** (a real number, the covering deficit density) is archimedean, a conditionally-convergent sinc-lattice sum, and has no Euler product. THM-503 stands.
- **The atomic skeleton** (the `(Z/N)*` maximizer set that `L_C` collapses to as the band `→ M`) has moments that ARE the multiplicative Ramanujan sums `c_N(k)`, which factor as an Euler product and carry `1/ζ(s)`.

So *the modular/multiplicative content lives in the skeleton (where the mass concentrates), the archimedean content lives in the collar widths (how much mass).* THM-503 was measuring the flesh; the units live in the bones.

## 3. The census is organized by which (Z/N)* the maximizers hit

The 78 eleven-cores `{1..13}\{i,j}`, sorted by `meas(L_C)` at band `1/14`, each carry a modular label — the `(Z/N)*` their maximizers realize:

| rank | core | `meas` | `M` | atoms on | `#`atoms | full orbit? |
|---|---|---|---|---|---|---|
| 1 | `{1..13}\{6,10}` | `313/9702 = 0.03226` | `1/10` | `(Z/10)*` | 4 | yes |
| 2 | `{1..13}\{4,6}` | `0.03238` | `2/19` | `(Z/19)` | 2 | no (2-clash) |
| … | … | … | … | … | … | … |
| top | `{1..13}\{7,9}` etc. | `0.147–0.155` | `1/7` | `(Z/7)*` | 6 | yes |

Two structural readings:

- **The minimizer is the small-`N` full-orbit core.** Removing `{10}` unlocks the `1/10` maximizers (no runner lands on `0`); removing `{6}` on top of it minimizes the collar volume. The result parks on `(Z/10)*` — a complete unit orbit of just 4 atoms with thin collars.
- **The runner-up is a partial-orbit two-speed clash** on `(Z/19)` (only 2 atoms, a Diophantine collision, not a full unit orbit) — the exact 11-speed echo of THM-523's tight-locus dichotomy `{AP dilates, Goddyn–Wong sporadic}`. The census bottoms out with *two competing families*, and they are within `0.0004` of each other.

This is the organizing principle the census needed: **the extremizer is not mysterious — it is the loosest core whose maximizers form a small complete unit orbit.**

## 4. The collar width is the modular inverse (a Dedekind-shaped units-sum)

Why those widths? At maximizer `a/N` the binding runner — the one with `‖v·a/N‖ = M` — has speed `v = a^{−1} mod N` (and its reflection `N − a^{−1}`). Verified across `(Z/10)*, (Z/12)*, (Z/14)*`: e.g. on `(Z/10)*`, atom `3/10` has `a^{−1}=7` and binding speeds `{3,7,13}`. So the collar half-width `(M−1/14)/v` is set by the **modular inverse of the atom**, and

`meas(L_C) = Σ_{a ∈ (Z/N)*} width(a)`,  with `width(a)` a function of `a^{−1} mod N`.

A sum over units of a function of `a` and its inverse is precisely the shape of a **Dedekind sum** `s(h,k) = Σ_a ((a/k))((ha/k))`. This is the structural bridge to HYP-3773 (opus), which proved rigorously that the covering-min margin is a 2-fold cotangent Dedekind sum with `s(n,Φ₆) → −1/12 = ζ(−1)`. The census measure and the covering-min margin are the *same kind of arithmetic object* — a units-sum of modular inverses — seen at two thresholds. (I do not re-derive `−1/12` here; I show the census measure has the Dedekind shape that HYP-3773 evaluates.)

So the two requested lenses meet on the collar: **modular multiplication picks the atoms and their binding inverses; the zeta/Dedekind value is what the resulting units-sum evaluates to.**

## 5. OPUC coordinates: a Verblunsky thermometer for "how tight"

Both search sweeps noted the repo has Toeplitz matrices but **no Verblunsky/OPUC computation**. Filling that gap (Levinson–Durbin on the moments `c_k = \hat 1_{L_C}(k)`):

- Every core gives `|α_k| < 1` — the Toeplitz form is positive-definite, `L_C` is a genuine (absolutely continuous) measure. Good sanity check on the whole moment picture.
- **The magnitudes are a thermometer for tightness.** The extremizer `{1..13}\{6,10}` is *near-atomic*: `|α_5|≈0.98, |α_9|≈0.95, |α_{11}|≈0.92`, Szegő product `≈10^{−5}`. The fatter core `{1..11}` stays cool: `max|α_k|≈0.65`, Szegő product `≈0.19`.

As the band `→ M` the measure collapses onto the `φ(N)` unit-atoms and some `|α_k| → 1` (the OPUC signature of an atomic measure). So **the extremizer is exactly the loose core that sits closest to the atomic/tight boundary** — the Verblunsky coefficients quantify "how close to `(Z/N)*`-atomic" a core is, and the minimizer is the most-atomic one. This is a new, coordinate-free handle on the census: minimize `meas` ⟺ maximize atomicity ⟺ push the Verblunsky coefficients to the unit circle.

## What this buys, and what it doesn't

**Buys (verified, new):** a single modular skeleton under both finish-tools — flat-extension atoms `= (Z/N)*`, moment matrix `=` Ramanujan sums, `ζ` as the unit-restricted Dirichlet denominator (walking through THM-503's open door without contradicting it), the census organized by unit orbits with the extremizer explained, collar widths `=` modular inverses (Dedekind-shaped, kin to HYP-3773's `−1/12`), and the first OPUC/Verblunsky coordinates in the repo as a tightness thermometer.

**Doesn't buy (honest):** a proof of `inf meas(L_C) ≥ 1/36`. The census here is a search over one family (2-drops of `{1..13}`); the modular labels *explain* the minimizer and *organize* the near-tight locus into `{full unit orbit, partial two-clash}`, but do not yet bound *all* primitive 11-cores. The natural next step the structure suggests: prove the two-family dichotomy is exhaustive (à la THM-523's tight-locus finiteness, but at band `1/14` over 11-speeds), then evaluate the two families' measures in closed form — the full-orbit one as a `(Z/N)*` Dedekind/collar sum, the two-clash one as a Diophantine collision — and check both clear `1/36`. Scale-invariance + quantization (THM-522) bounds the search to primitive scale-1 cores, so the census is finite.

— Related: HYP-3789 (flat-extension atoms = units; this makes the moment matrix explicit as Ramanujan sums), HYP-3773 (margin = 2-fold Dedekind sum → −1/12 = ζ(−1); the collar sum shares its shape), THM-503 (no Euler product for `L` — refined here: the unit-restricted skeleton *does* carry 1/ζ), THM-501/515 (singular series = lonely measure), THM-523 ({AP, Goddyn–Wong} tight-locus dichotomy → the 11-speed {full-orbit, two-clash} echo), THM-560 (census enumeration template), THM-522 (scale-invariance + quantization bounds the census), HYP-3762 (three-gap / CF-convergent geometry of the atoms), OPEN-Q-108. Scripts: `04-computation/lrc14_census_toeplitz_ramanujan_kps.py`, `lrc14_census_modular_verblunsky_kps.py` (+ `.out`). HYP-3792.
