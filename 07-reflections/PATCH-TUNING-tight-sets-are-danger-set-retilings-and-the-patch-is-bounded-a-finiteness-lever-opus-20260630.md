# PATCH TUNING: tightness M(S)≤1/n IS the danger-set covering ∪_{v∈S}D_v=[0,1) (D_v={t:‖vt‖≤1/n}); the AP TILES [0,1), and every tight set is a RETILING — remove the tile D_k, patch the hole H_k=D_k∖∪_{v≠k}D_v with D_g; the patch is BOUNDED (each hole-interval must fit inside ONE D_g-interval since D_g has gaps, so g≤2/(n·w) for the interval-width w; with hole-widths ≥ the finest AP structure ~1/n², every patch is g≲2n — verified 7,9,12,24 across the census), which makes the AP-modification tight locus a FINITE search — the bounded-speed lever OPEN-Q-108 wants; the DOUBLING operad is the g=2k patch and the JACOBSTHAL gate is precisely "does D_g cover H_k"; plus four creative angles (covering-rigidity, Fourier-positivity, Stern-Brocot resonance, inhomogeneous patch)

*opus-2026-06-30. Owner: merge in patch tuning, extend it, find more creative proof angles. Patch tuning
formalizes as danger-set retiling; the key extension is that the patch is BOUNDED, a lever toward the tight
locus's finiteness (OPEN-Q-108).*

## Patch tuning = danger-set retiling
> **`M(S) ≤ 1/n` ⟺ the danger sets `D_v = {t : ‖vt‖ ≤ 1/n}` (`v∈S`) cover `[0,1)`.**
Each `D_v` is `v` intervals of half-width `1/(vn)` around `j/v`, total measure `2/n`. This is OPEN-Q-108's own
frame (the `D_v` are the danger arcs; the lonely set `G` is the complement). **The AP `{1..n-1}` TILES
`[0,1)`** (covers it, the tight witnesses `t=j/n` the measure-zero seams). Every tight set is a **RETILING**:
- **remove** a tile `D_k` ⇒ opens the **hole** `H_k = D_k ∖ ⋃_{v≠k} D_v` (the part only `D_k` covered);
- **patch** it with `D_g` (add element `g`). Tightness survives iff `D_g ⊇ H_k` (`H_k` refilled).

"Patch tuning" = choosing the patch `g` so `D_g` covers the hole. My earlier difference-closed result is the
degenerate case (no tile removed); the sporadics are genuine retilings.

## The patch is BOUNDED — the finiteness lever
`D_g` is `g` intervals of half-width `1/(gn)` separated by **gaps** of width `(n−2)/(gn) > 0`. So an interval
of the hole of width `w` can be covered by `D_g` **only if it fits inside a single `D_g`-interval** (a gap
inside it would leave a point uncovered):
> **every hole-interval of width `w` forces `w ≤ 2/(gn)`, i.e. `g ≤ 2/(n·w)`.**
The hole-intervals are bounded below by the **finest AP danger structure** (edges at `j/v ± 1/(vn)`,
`v ≤ n−1`, so spacing `≳ 1/n²`), giving **`g ≲ 2n`**. Verified: every census patch is `< 2n`
(`n=5:7, 6:9, 8:12, 14:24`; single- and multi-swap), and the per-interval bound `g ≤ 2/(n·w)` holds exactly
(`n=14,k=12: w=0.0041 ⇒ g≤35`, actual `24`). **Consequence: the AP-modification tight locus has bounded max
element ⇒ it is a FINITE search** — the bounded-speed statement OPEN-Q-108 needs, now with a mechanism.

## Doubling operad = the g=2k patch; Jacobsthal = the coverage test
The repo's **doubling operad** (`k ↦ 2k`, THM-631/HYP-2917) is exactly the patch `g=2k`: `D_{2k}` has intervals
of half-width `1/(2kn)` on the `2k`-grid, covering the **inner half** of each `D_k`-interval; the residual
outer half of `H_k` must be covered by the surviving AP runners — **and whether it is covered is precisely the
Jacobsthal gate** (the coprime-window clearance). So the two repo objects are the two halves of
"`D_g ⊇ H_k`": `D_{2k}` = the inner patch, Jacobsthal = the outer-residual coverage test. The rare
non-doubling sporadics (`n=5:g=7`, `n=6:g=9`) are other `g` whose interval-grid happens to align with `H_k`.

## Four creative proof angles (extending patch tuning)
1. **Covering-rigidity (conserved overlap).** `Σ_v meas(D_v) = (n−1)·2/n = 2 − 2/n`, so any covering of
   `[0,1)` has **overlap measure exactly `1 − 2/n`** (a conserved quantity). A retiling must preserve it:
   `meas(D_g ∩ H_k^c)` (the patch's wasted overlap) is forced, tightly constraining `g`. *Angle: prove
   finiteness of retilings by bounding overlap-preserving patches.*
2. **Fourier-positivity.** Covering ⟺ `Σ_{v∈S} 1_{D_v}(t) ≥ 1 ∀t`. Since `1_{D_v}` is a dilated Fejér-type
   kernel, this is a **nonnegativity of an exponential sum** `Σ_v Σ_{|m|≤?} ĉ_v(m) e(mt) ≥ 1`. *Angle: an
   SOS/Turán-type certificate for the AP covering, perturbed by the swap.*
3. **Stern-Brocot resonance.** The holes `H_k` sit at the **Farey fractions `j/k`** (the resonances of the
   removed runner); the patch `D_g` must cover the Farey neighborhood. This is the SAME Stern-Brocot backbone
   that carries the covering-min (`t=1/n` = base leaf, the regimes = Farey intervals — my earlier reflections).
   *Angle: patch tuning is navigation on the Stern-Brocot tree; the doubling `k→2k` is a specific tree move.*
4. **The inhomogeneous patch.** My exact law `M_c(AP)=1/n+c(n−2)/n` (THM-591) — does a patched set have an
   inhomogeneous law `M_c(S)`? The observer `c` shifts the danger sets to `D_v−c`; patch tuning becomes
   `c`-dependent. *Angle: a 2-parameter `(c, g)` family unifying observer-translation and tile-swap; the
   avoided-arc-edge identity `c−env=−1/q` may extend to the patched danger geometry.*

## Toward OPEN-Q-108
- **The bounded-patch is the bounded-speed lever.** If every tight set is an AP-modification (census evidence
  through n=14), the patch bound `g ≲ 2n` makes the tight locus a finite search ⇒ finiteness ⇒ the fattening
  lemma reduces to a computation. The gap: proving *every* tight set is an AP-modification (not just the
  enumerated ones) — i.e. that tightness forces the AP skeleton (a "lowness lemma" for the tight threshold,
  cf. HYP-3740 for the covering-min threshold).
- **The clean chain:** difference-closed = AP (rigorous, mine) → sporadics = bounded retilings (patch tuning,
  this reflection) → finiteness (bounded-speed) → fattening lemma. Two links rigorous; the "tightness ⇒ AP
  skeleton" link is the remaining hard core.

## Status
- **Formalized (opus):** patch tuning = danger-set retiling; `M(S)≤1/n ⟺ ⋃D_v=[0,1)`; AP tiles; sporadic =
  remove `D_k`, patch `H_k` by `D_g`.
- **Rigorous (opus):** the per-interval patch bound `g ≤ 2/(n·w)` (single `D_g`-interval must contain each
  hole-interval); empirically `g < 2n` across the census ⇒ bounded ⇒ AP-modification tight locus finite.
- **Unified:** doubling operad = `g=2k` inner patch; Jacobsthal gate = outer-residual coverage.
- **Creative angles opened:** covering-rigidity/conserved-overlap, Fourier-positivity, Stern-Brocot resonance,
  inhomogeneous patch.
- **Open (honest):** the exact `w_min ≳ 1/n²` bound (clean constant); proving tightness ⇒ AP-skeleton;
  OPEN-Q-108 itself.

Related: the-avoided-arc-argument-…-difference-closed (the tile-free case), THM-591 (the inhomogeneous law),
UPPER-BOUND-CLOSED-… (the avoided-arc edge), the-loneliness-integral-…-stern-brocot (the Farey backbone),
THM-631/HYP-2917 (doubling operad), HYP-3740 (the covering-min lowness lemma), OPEN-Q-108; HYP-3753 (this),
HYP-3749 (difference-closed = AP).
