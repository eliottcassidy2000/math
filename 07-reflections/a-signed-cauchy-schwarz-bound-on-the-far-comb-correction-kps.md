# A signed Cauchy–Schwarz bound on the far-comb correction: the fat-core multi-far closes, and the residual is near-tight cores × non-resonant combs

*kind-pasteur-2026-07-01. Aiming creatively at the signed resonance correction (the localized singular series) that is the multi-far residual. The absolute bound diverges (MISTAKE-078); a plain Cauchy–Schwarz against the FIXED lonely set, powered by Parseval, is finite and beats the main term — closing the fat-core part and isolating the residual to the near-tight (small-lonely-measure) cores against incommensurate combs.*

## The object and the wall

Multi-far survival factors as `survival(r) = (6/7)^r · meas(L_C) − [signed resonance correction]`, the singular series localized to the `r` far combs on the fixed core lonely set `L_C`. Writing `1_safe = 6/7 + f` (`f` mean-zero, Fourier `c_k`, `|c_k|~1/(π|k|)`), the two-comb correction is

`correction₂ = ∫_{L_C} f(W₁t) f(W₂t) dt`.

Its **absolute** bound is `Σ_{k₁,k₂}|c_{k₁}c_{k₂}||\hat 1_{L_C}(k₁W₁+k₂W₂)| ≤ meas(L_C)·(Σ|c_k|)²`, and `Σ|c_k| ~ Σ 1/k` **diverges** — the MISTAKE-078 wall, the reason positivity/absolute methods can't reach the correction.

## The signed move: Cauchy–Schwarz against the fixed `L_C`, via Parseval

The correction is an inner product against the fixed set. Cauchy–Schwarz, with **Parseval `‖1_{L_C}‖₂² = meas(L_C)`**:

> `|correction₂| = |⟨1_{L_C}, g⟩| ≤ ‖1_{L_C}‖₂ · ‖g‖₂ = √(meas(L_C)) · ‖g‖₂`, `g(t)=f(W₁t)f(W₂t)`,
> and `‖g‖₂² = ∫ f(W₁t)² f(W₂t)² dt → (∫f²)² = (6/49)²` for large non-resonant `W₁,W₂`.

This is **finite** — it converts the divergent term-by-term sum into a single `L²` pairing that the harmonic tail never sees. Hence

> `survival(2) ≥ (6/7)² meas(L_C) − √(meas(L_C))·(6/49) > 0  ⟺  meas(L_C) > 1/36`.

Verified exactly (`lrc14_signed_correction_cauchyschwarz_kps.py`): on cores of `meas(L_C)=0.056, 0.117, 0.167`, the true correction is tiny (`~0.002`), well inside the CS bound (`0.029, 0.042, 0.050`), and `survival ≥ main − bound = +0.012, +0.026, +0.072 > 0`. (A `49×` slip in the `‖g‖₂` computation, caught and fixed, made the bound *stronger* once corrected: `‖g‖₂ = 6/49` exactly.)

## The three-way reduction of two-far

`‖g‖₂` inflates for **low-denominator ratios** `W₂/W₁` (verified): `2:1 → 0.203` (worst, threshold `meas>1/13`), `3:1,3:2 → 0.170`, higher → `6/49`. So:

- **(A) Fat cores** `meas(L_C) > 1/13`: **closed for ALL combs** by CS (even the worst 2:1 resonance).
- **(B) Resonant combs** (low-denominator `W₂/W₁`): `W₂=2W₁` etc. make `safe(W₁)∩safe(W₂)` a function of `W₁t` alone — the two combs **collapse to one effective comb** → single-far, closed by the band-barrier (and it's the three-distance "one rotation" case).
- **(C) Residual**: **near-tight cores** (`meas(L_C) < 1/36`, the OPEN-Q-108 low-lonely-measure locus — the tight/near-tight sets, `inf L = 1/1260`) against **non-resonant (incommensurate) combs**. This is the genuine hard core, now sharply isolated: *small measure AND incommensurate* simultaneously.

## The r-far picture and why it concentrates the residual

For `r` combs, `‖g_r‖₂ → (6/49)^{r/2}` (non-resonant), so the threshold is `meas(L_C) > 6^{-r}` — it **shrinks with r**. And the core has `13−r` speeds, so as `r` grows the core **shrinks and `L_C` fattens** (`meas → 1`). Both effects push toward closure: higher `r` closes more easily. The residual therefore concentrates at **low `r` (`r=2,3`)** on the **near-tight cores**. (`r≥7` is THM-573.)

## What this buys, honestly

- **New (verified):** a *signed* bound on the far-comb correction — `L²`/Parseval Cauchy–Schwarz against the fixed lonely set — that is **finite where the absolute bound diverges**, with the explicit threshold `meas(L_C) > 1/36` (non-resonant) / `1/13` (worst 2:1). It closes the fat-core and resonant parts of multi-far outright.
- **The reduction:** OPEN-Q-108's multi-far residual is now **near-tight cores × non-resonant low-`r` combs** — the *small-lonely-measure* locus (the tight-locus census, `inf L=1/1260`, my `lrc14-thread`) crossed with incommensurate far speeds. Not "does a beater exist," but "on the finitely-characterized near-tight cores, is the incommensurate far-comb correction below the (small) main term." A far tighter target.
- **Not a proof.** The near-tight × non-resonant residual is unresolved (it IS OPEN-Q-108's heart). But the CS bound peels off everything else and connects the residual to the two things the project already understands best: the **tight-locus census** (the cores) and the **decorrelation/Riesz-product** program (the combs). This is exactly the L²-signed tool of HYP-3129 (SPEC) transported to the far-comb setting, and the near-equal-comb worst case is the THM-578 doublet.

## Net

The signed correction that stumped the absolute methods is bounded by a Cauchy–Schwarz against the fixed `L_C`: `|correction| ≤ √meas·‖g‖₂`, finite, beating the main term for `meas(L_C) > 1/36`. Multi-far splits into [fat cores: closed] + [resonant combs: collapse to single-far] + [near-tight cores × incommensurate combs: the OPEN-Q-108 heart, now sharply isolated and wired to the tight-locus census and the Riesz-product program].

— Related: the prior reflections `multi-far-is-the-singular-series-localized-...`, `the-single-far-unbounded-case-closes-...`, `the-real-dichotomy-is-bounded-vs-unbounded-...` (this session's arc); HYP-3129 (L²-CS on SPEC — same tool), THM-578 (doublet = near-equal combs), THM-504 (the divergence wall), `lrc14-is-the-lonely-measure-and-the-key-is-a-riesz-product.md`, [[lrc14-thread]] (tight locus, inf L=1/1260), klein HYP-3784 (Delsarte dual), OPEN-Q-108. Script: `04-computation/lrc14_signed_correction_cauchyschwarz_kps.py`. Not a HYP reservation.
