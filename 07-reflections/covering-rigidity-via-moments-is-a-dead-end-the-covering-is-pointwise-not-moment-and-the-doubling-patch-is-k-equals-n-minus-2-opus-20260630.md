# Chasing covering-rigidity vigorously — an HONEST NEGATIVE: the conserved overlap ∫(m−1)=1−2/n is satisfied by EVERY tight set (not a distinguishing rigidity), and the second moment ∫m² does NOT separate tight from non-tight (the AP is NOT the smoothest covering — n=8 sporadic has ∫m²=4.60<4.78=AP; the NON-tight 12→26 swap has the LOWEST 6.11<6.13=GW<6.33=AP at n=14), so MOMENT-based covering-rigidity cannot force the AP-skeleton; the reason is structural — tightness is the POINTWISE condition m(t)≥1 ∀t, which moments average away — so the right tool is the pointwise patch-tuning (widest-hole ⇒ bounded patch), not moments; the salvage is a clean SHAPE result: the doubling patch is specifically k=n−2 → 2(n−2) (necessary: 2k∉AP i.e. k≥n/2; the tight one is k=n−2, verified n=8,14), Jacobsthal-gated

*opus-2026-06-30. Owner: chase covering-rigidity vigorously. I did — and the moment form is a dead end for the
lowness. Recording the negative + why + the shape salvage, so no one re-runs it.*

## The angle, and why it fails
"Covering-rigidity" (patch-tuning Angle 1): tightness `M(S)≤1/n` ⟺ danger sets `D_v` cover `[0,1)`, and
`∫m(t)dt = Σmeas(D_v) = 2−2/n` gives a **conserved overlap `∫(m−1)=1−2/n`** (`m(t)=#{v:‖vt‖≤1/n}`). The hope:
this budget forces the AP-skeleton. **It does not.**
- **The conserved overlap is not a rigidity.** `∫(m−1)=1−2/n` holds for *every* set whose danger sets cover
  `[0,1)` — it's the definition of covering with fixed per-set measure `2/n`, not a property of the AP.
- **The second moment doesn't separate tight from non-tight.** `∫m²` (covering "smoothness"): at n=14,
  `AP=6.33`, `GW(tight)=6.13`, and the **non-tight** `12→26` swap `=6.11` — the *non-tight* set is the
  smoothest. At n=8, `AP=4.78 > sporadic=4.60`. So the AP is **not** the min-`∫m²`, and low `∫m²` does not
  imply tight. Moments are blind to the distinction.

## Why moments are the wrong tool (the real lesson)
Tightness is `min_t m(t) ≥ 1` — a **pointwise** (`L^∞`) condition: *one* uncovered point (the widest hole)
breaks it. `∫m`, `∫m²` are **averages** (`L^1, L^2`) that wash out exactly that pointwise structure. A set can
have small moments yet leave a hole (non-tight), or large moments yet cover (tight). **So no moment functional
of `m` can characterize the tight locus.** The right object is the *pointwise* covering — which is exactly the
**patch-tuning** analysis (the widest hole-interval must fit in one `D_g`-interval ⇒ the bounded patch). That
stands; the moment detour does not add to it.

## The salvage: the doubling patch is k = n−2
Chasing the *shape* (not the moments) did clarify the doubling operad. A single-swap doubling `k→2k` needs
`2k ∉ AP` — i.e. **`k ≥ n/2`** (else `2k` is already a runner, not a new patch). Among those, the **tight** one
is `k = n−2`:
| n | 8 | 10 | 12 | 14 | 16 |
|---|---|---|---|---|---|
| `2k∉AP` candidates (`k≥n/2`) | 4–7 | 5–9 | 6–11 | 7–13 | 8–15 |
| doubling-TIGHT `k` | **6** | — | — | **12** | — |
So the doubling operad is precisely the **`n−2 → 2(n−2)`** swap (`n=8: 6→12`, `n=14: 12→24`), and it is tight
only when the **Jacobsthal gate** opens for `v=n−2` (rare: `n=8,14` among `n≤16`). `k=n−2` is the largest `k`
whose double `2(n−2)` still lands where the surviving AP runners can cover the outer residual — the geometric
content of the Jacobsthal window. This confirms & sharpens the repo's picture (`v=12` at n=14) and gives the
necessary condition `k≥n/2` cleanly.

## Consequences / redirection
- **Angle 1 (covering-rigidity via moments) is REFUTED as a lowness lever.** Do not pursue `∫m^p` functionals
  to force the tight-skeleton — they average away the pointwise covering.
- **The pointwise patch-tuning is the surviving tool** for finiteness (bounded patch) and for the lowness
  (which holes are unpatchable). The remaining hard link ("tightness ⇒ AP-skeleton") must be attacked
  pointwise, not by moments — likely via the **hole geometry** (which `k` open unpatchable holes) rather than
  any averaged quantity.
- **Better creative angles remain**: patch-tuning Angles 2–4 (Fourier-*positivity* — note this is `L^∞`/pointwise,
  NOT a moment, so still live; Stern-Brocot resonance; inhomogeneous patch). The Fourier-positivity angle is
  the natural pointwise successor: `Σ_v f(vt) ≥ 1` with `f`'s spectrum `\hat f(m)=sin(2πm/n)/(πm)` (vanishing on
  `m≡0 mod n`) — a genuine min-of-exponential-sum condition, not an average.

## Status
- **Negative (opus, verified):** conserved overlap `1−2/n` and `∫m²` do not distinguish tight sets; moment
  covering-rigidity cannot force the AP-skeleton (covering is pointwise `L^∞`, moments are `L^1/L^2`).
- **Salvage (opus, verified n=8,14):** doubling patch = `n−2 → 2(n−2)`; necessary `k≥n/2`; Jacobsthal-gated.
- **Redirection:** attack the lowness pointwise (hole geometry / Fourier-positivity), not via moments.

Related: PATCH-TUNING-… (Angle 1 refuted here, Angles 2–4 live; the pointwise patch bound stands),
the-avoided-arc-argument-… (difference-closed = AP), THM-631/HYP-2917 (doubling operad, v=n−2), OPEN-Q-108.
