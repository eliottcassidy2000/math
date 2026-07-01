# Deep reasoning + Fourier-positivity chase: the LRC tight-skeleton is COMBINATORIAL COVERING (a multiple of every q∈{2,…,n−1}), NOT an averaged/energy condition — Fourier-positivity Σf(vt)≥1 LOCALIZES to the q-witness at t=1/q (it IS THM-523, not a genuine analytic min: the covering-min sits at low-denominator rationals like t=1/7 where only runner 7 covers), so the whole "analytic" route reduces to the combinatorial covering; and BOTH averaged lenses are DEAD ENDS with clean counterexamples — ∫m² doesn't separate tight from non-tight (the non-tight 12→26 swap has the LOWEST ∫m²), and additive energy REFUTED (AP+1={2..14} has MAXIMAL additive energy yet is non-tight, contains 14≡0 mod 14) — because tightness is POINTWISE (L∞: min_t m(t)≥1) and averaging washes it out; the repo's developed lens is exactly this pointwise view (Fourier-Toeplitz PSD certificates HYP-2974; singular-series THM-523; the PROVED+Lean binding-pair HYP-+2909: tight ⟹ a pair sums to 0 mod n), and the doubling operad is a COVERING-PRESERVING multiple-swap (replace runner k, which covers q=k, by 2k, which still covers q=k), Jacobsthal-gated

*opus-2026-06-30. Owner: chase Fourier-positivity, reason deeply about other possibilities, search the repo
for distinct lenses. The chase + the search converge: the tight-skeleton is combinatorial covering; averaged
lenses can't see it; the analytic route localizes to the q-witness.*

## The deep question and the answer
"What IS the tight-skeleton (which S have `M(S)=1/n`)?" Answer, from the chase: **a COVERING condition.**
`M(S)=1/n < 1/q` for every `q<n` forces (q-witness, THM-523) **a multiple of every `q∈{2,…,n−1}` in S**
(else `t=1/q` gives `M≥1/q`). Verified: AP, GW cover all `q`; the non-tight `12→26` misses `q=12`, `drop-7`
misses `q=7`. The covering-min `min_t D(t)` (`D=#{v:‖vt‖≤1/n}`) is achieved at **low-denominator rationals**
(`t=1/7` for AP & GW, where only runner `7` sits in the danger zone, `D=1`) — so tightness is a
**finite/combinatorial** condition at the Farey resonances, not a genuine analytic minimum.

## Fourier-positivity chased ⇒ it localizes to the q-witness (= HYP-2974 / THM-523)
`M(S)≤1/n ⟺ Σ_{v∈S} f(vt) ≥ 1 ∀t` (`f=1_{‖·‖≤1/n}`, `f̂(m)=sin(2πm/n)/(πm)`, vanishing on `m≡0 mod n`).
Chasing this **rediscovered the repo's Fourier-Toeplitz PSD lens** (HYP-2974): `D_S−1≥0` a.e. ⟺ the Toeplitz
moment matrices are PSD; a negative eigenvalue is a Farkas certificate of a hole. **AP/GW are the
borderline-PSD (extremal) cases** — exactly why the Beurling–Selberg minorant is delicate there (the AP is
tight to zero margin, so any minorant loss breaks it). The spectrum `D̂(k)=Σ_{v|k,v∈S} f̂(k/v)` is a
**divisor sum** — the multiplicative/covering structure, again pointing at combinatorial covering, not an
average. **Fourier-positivity is the right framing but reduces to the combinatorial covering (the q-witness).**

## The averaged lenses are DEAD ENDS (clean counterexamples — the real lesson)
Tightness is `min_t m(t) ≥ 1` — **pointwise (L∞)**. Averaged functionals wash out the one point (the widest
hole) that matters:
- **Moments (`∫m²`):** does NOT separate tight from non-tight — the non-tight `12→26` swap has the *lowest*
  `∫m²` (6.11 < 6.13 GW < 6.33 AP); the AP is not even the smoothest.
- **Additive energy (`A(S)=∫|T_S|⁴`):** **REFUTED as a tight-distinguisher** — `AP+1={2,…,14}` has *maximal*
  additive energy (it's a shifted AP) yet is non-tight (`M=1/8`, contains `14≡0 mod 14`); and GW (tight) has
  the *same* `A` as its non-tight neighbor `12→26`. (The repo's HYP-+2873 correctly uses `A` for the AP's
  *extremality*/Fejér-optimality — NOT as a tightness criterion; consistent with this refutation.)
> **META-INSIGHT: no averaged/energy functional of the runner set can characterize the tight locus.** Two
> clean counterexamples prove it. The tight-skeleton is irreducibly pointwise. Do not chase `∫m^p`, additive
> energy, or spectral-moment characterizations of tightness — they are provably blind.

## The repo's lens-map, re-read through this (search integration)
The developed lenses are ALL the pointwise/combinatorial view:
- **THM-523 singular-series** (most developed): the covering / q-witness — exactly my "tight ⟹ covering."
- **HYP-2974 Fourier-Toeplitz PSD** (degree-160 Fejér certificates): the pointwise Fourier-positivity;
  AP/GW the borderline cases — what my chase landed on.
- **HYP-+2909 binding pair (PROVED + Lean, sorry-free):** `M(S)=1/n ⟹` some pair `s_i+s_j≡0 mod n` at `t*`;
  AP & GW both bind `{1,13}`. The rigorous *forward* half of the skeleton.
- **HYP-3750 near-AP classification / HYP-2893 Jacobsthal accelerations:** sporadics = "duplication+drop"
  (drop `v`, add `2v`), Jacobsthal-gated (`n≡1 mod 6` for the GW law), plus cross-type — my census.
- **HYP-3740 lowness lemma:** the covering-min threshold analog; the isolated hard core.
> The **doubling operad is covering-preserving**: replacing runner `k` (which covers `q=k`) by `2k` (which
> STILL covers `q=k`, being a multiple) keeps S covering — *that* is why it's a legal retiling; Jacobsthal
> gates which `q` may substitute `2q` for `q` (the finer rational conditions). This unifies THM-523 (covering)
> with HYP-2893 (Jacobsthal) via one idea: **tight sets are minimal-ish coverings, deviating only by
> covering-preserving multiple-swaps.**

## Verdict and the crux
- **The tight-skeleton = combinatorial covering** (mult of every `q`) refined by finer Farey/Jacobsthal
  conditions. Analytic (Fourier) ⇒ localizes to it; averaged (moments/energy) ⇒ blind to it.
- **The lowness ("tight ⟹ AP-skeleton")** = "covering forces a multiple of every `q`; the minimal covering is
  the AP (`q` itself); deviations are bounded (patch-tuning) & gated (Jacobsthal) covering-preserving swaps."
  The remaining hard link is proving the deviations are *finite* — which is HYP-3740's lowness (the repo's
  most-promising-underexplored lens) and OPEN-Q-108's fattening.
- **Right next attack:** the pointwise Fourier-Toeplitz PSD certificate (HYP-2974) or the lowness lemma
  (HYP-3740) — NOT any averaged quantity. My patch-tuning bounded-patch is the geometric companion (finiteness
  once the AP-skeleton is forced).

## Status
- **Refuted (opus, clean):** moment `∫m²` and additive-energy characterizations of tightness (counterexamples:
  non-tight lower `∫m²`; `AP+1` max energy non-tight). Meta: averaged lenses are provably blind (pointwise).
- **Rediscovered/confirmed:** Fourier-positivity = the q-witness / covering (THM-523) = Fourier-Toeplitz PSD
  (HYP-2974); AP/GW borderline. Divisor-sum spectrum `D̂(k)=Σ_{v|k}f̂(k/v)`.
- **Unified:** doubling operad = covering-preserving multiple-swap (`k→2k` still covers `q=k`); Jacobsthal =
  the finer gate. Ties THM-523 + HYP-2893 + HYP-+2909.
- **Crux:** the pointwise lowness (HYP-3740) / fattening (OPEN-Q-108); attack pointwise, not averaged.

Related: covering-rigidity-…-dead-end (the moment negative), PATCH-TUNING-… (the geometric companion),
THM-523 (q-witness/singular-series), HYP-2974 (Fourier-Toeplitz PSD), HYP-+2909 (binding pair, proved),
HYP-2893/3750 (Jacobsthal/near-AP), HYP-3740 (lowness), OPEN-Q-108; HYP-3753 (patch tuning).
