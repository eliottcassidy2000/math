# The modular-forms bound on the covering-min rung a(n) FAILS — beater EXISTENCE is combinatorial, not cuspidal: the genus of X₀(2n) does not track beaters (they occur at genus 0, n=8,9; the construction wins at genus 1,2, n=12,13,14), the cusp form f₁₄=(14a) coefficients [1,−1,1,0,0] do NOT match a(n)=[2,2,4,4,3], and the construction-vs-covering-min margin gaps are irregular rationals with no L-value/norm pattern; integrating klein-S60's ILP-verified transition (HYP-3778: a(n)=n for n≥12 via the radius-1 band), the rung a(n) is a FINITE combinatorial phenomenon (n=7..11), while the modular structure (−1/12=ζ(−1), the Eisenstein/E₂ anomaly) governs ONLY the construction's margin VALUE — so the clean law is **VALUE is modular, EXISTENCE is combinatorial**; the cusp form's real role is the OPEN-Q-108 proof obstruction (klein-S56), not a(n), and the rigorous transition-proof tool is character-sum/Weil (the repo's √V cancellation, THM-546), not modular forms

*opus-2026-07-01. Owner: work on the modular-forms bound (a cusp-form-norm bound on a(n)). I tested it three
ways — it fails cleanly — and the failure sharpens the framework into a value/existence separation, integrating
klein-S60's freshly ILP-verified transition.*

## The test: three independent NEGATIVES
The HYP-3775 framework suggested the "residual" (beaters, the covering-min rung a(n)) might be the genus-1 cusp
form of X₀(2n), giving a modular-forms bound on a(n). Tested against klein's data (a(n)=2,2,4,4,3 for n=7..11;
a(n)=n for n≥12, HYP-3778):
1. **Genus mismatch.** `genus X₀(2n)` for n=7..14 = `1,0,0,1,2, 1,2,2`. Beaters exist at n=7..11 (genus
   `1,0,0,1,2` — **including genus 0** at n=8,9, where there is NO cusp form); the construction wins at n=12..14
   (genus `1,2,2` — cusp form EXISTS, yet no beater). The genus does **not** separate beater/no-beater.
2. **Coefficient mismatch.** f₁₄ = η(q)η(q²)η(q⁷)η(q¹⁴) (curve 14a) has `aₙ = 1,−1,1,0,0` for n=7..11 — nothing
   like `a(n)=2,2,4,4,3` (e.g. f₁₄'s a₁₀=a₁₁=0 but a(10)=4, a(11)=3).
3. **Gap irregularity.** The margin gap `construction − covering-min` for n=7..11 is `5/559, 2/285, 5/2409,
   6/3367, 8/3441` — irregular rationals, all Θ(1/n²), no L-value / Petersson-norm pattern; and it drops
   **abruptly to 0** at n=12 (a sharp threshold, not a modular decay).
**Conclusion: there is no modular-forms bound on a(n).** The covering-min rung is not modular.

## Why: value vs existence (the clean law)
klein-S60 (HYP-3778) has now **ILP-verified the transition is REAL** — a(n)=n for n≥12 (up to speeds 4n),
correcting the earlier H3 pessimism — and its mechanism is the **radius-1 band over-constraint** (HYP-3737): a
spread (disordered) base cannot cover `Z/D` at radius 1 for the band `D∈(n,2n−2]`, forcing the consecutive
core + outlier. That is a **covering/combinatorial** fact. So:
> **EXISTENCE is combinatorial.** Whether a disordered set beats the construction (a(n)<n) is governed by the
> radius-1 band, a finite covering condition. Beaters are a FINITE phenomenon (n=7..11); their rungs 2,2,4,4,3
> are irregular/non-modular (no closed form — as a genuinely arithmetic sequence should be, HYP-3732).
> **VALUE is modular.** The CONSTRUCTION's margin `= −12 s(n,Φ₆)/n² → −1/12 = ζ(−1)` is the Eisenstein/E₂
> anomaly (HYP-3768/3774/3775) — the Dedekind/eta cocycle of the *ordered* (scaled-interval/AP) set. This is
> the only genuinely modular part, and it is about the VALUE, not the existence of a competitor.
The klein/mac-mini "Eisenstein bulk + cusp residual" framework conflated these; the test separates them.
HYP-3775's criterion (margin is a Dedekind sum ⟺ speeds are a scaled interval) stands — it is a statement about
the VALUE. The "beaters = cusp form" reading (about existence) is **withdrawn**.

## The cusp form's real role: the PROOF obstruction, not a(n)
So what IS the genus-1 cusp form f₁₄? klein-S56: the RESIDUAL = the **ι-odd genus cusp form = OPEN-Q-108's
Borsuk–Ulam odd index** — i.e. the obstruction to *proving* the covering-min structure (the ι-odd degree that
must vanish), NOT the value a(n). This is consistent: a(14)=14 (no beater, klein-S60), so at n=14 the
covering-min IS the pure Eisenstein construction — there is no beater "residual" to be a cusp form; the cusp
form is the residual of the PROOF (why no beater exists), which is genuinely hard (genus 1 ⇒ nontrivial).

## The rigorous path (the redirect): character sums / Weil, not modular forms
klein's transition is ILP-verified only up to speeds 4n; the open part is "no beater with speeds > 4n." That is
a **covering-realizability** statement, and its natural rigorous tool is a **character-sum / Weil bound** — the
*elementary* Ramanujan–Petersson (square-root cancellation), which the repo already carries as the peel-deviation
`|Δ_w·w| ~ C√V` bound (THM-546 / HYP-2852). "Cusp form = square-root cancellation = Weil bound" is the correct
shadow, but it lives at the level of the **danger-cover fluctuation**, not the modular curve X₀(2n). So the
modular-forms *language* (Eisenstein/cusp) is a useful picture for the VALUE (−1/12), but the transition PROOF
must be done with character sums on the band, not cusp-form norms.

## Status
- **Refuted (opus, three ways):** the modular-forms bound on a(n) — genus mismatch (beaters at genus 0; no
  beaters at genus 1,2), f₁₄-coefficient mismatch, irregular gaps with abrupt threshold. No modular bound on a(n).
- **Integrated (klein-S60/HYP-3778):** the transition a(n)=n for n≥12 is ILP-verified (radius-1 band,
  combinatorial); beaters are finite (n=7..11), non-modular.
- **Sharpened:** the clean law is **VALUE is modular (−1/12/Eisenstein, the construction), EXISTENCE is
  combinatorial (the band, a(n))**; HYP-3775's Dedekind⟺interval criterion is a VALUE statement and stands; the
  "beaters = cusp" existence reading is withdrawn.
- **Redirected:** the cusp form f₁₄ is the OPEN-Q-108 proof obstruction, not a(n); the rigorous transition proof
  (klein's >4n) is a character-sum/Weil problem (the repo's √V cancellation, THM-546), not modular forms.
- **Open:** the transition proof for speeds >4n (klein's HYP-3778/3737 open part) via character sums; OPEN-Q-108.

Related: HYP-3778/klein-S60 (transition real, ILP), HYP-3737 (radius-1 band mechanism), HYP-3775/opus-S7
(Dedekind⟺interval — the VALUE criterion, stands), HYP-3774/3768 (−1/12/Eisenstein, the VALUE), HYP-3586
(genus jump), THM-546/HYP-2852 (√V peel-deviation = the elementary Ramanujan–Petersson shadow), OPEN-Q-108.
HYP-3779 (this). Scripts: 04-computation/lrc_modular_bound_test_genus_f14_opus_20260701.py.
