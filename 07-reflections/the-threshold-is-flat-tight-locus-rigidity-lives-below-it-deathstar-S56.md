# The threshold is flat: tight-locus rigidity lives below it

**death-star-2026-07-17-S56.** A structural finding from generalizing the live law (THM-991 → THM-996):
the loneliness threshold `1/14` is a *flat* invariant — it cannot see the difference between the two
primitive tight families. The information that distinguishes them, and therefore the information any
rigidity proof must consume, lives strictly *below* the threshold, in the sub-threshold coverage
spectrum. This sharpens "what we are missing" in our view of the LRC(14) object and prunes a whole
class of attacks on the hard core.

Companions: THM-996 (the census law), THM-991 (the n=14 live law), LEM-020 (moment-order-3 necessity),
`the-coverage-spectrum-one-grammar-four-instruments.md`, `the-lrc14-proof-object-what-every-lens-pins-death-star-S17.md`.

---

## 1. Two facts that only make sense together

**Fact A (resonance confinement, general).** Every tight family (`M = 1/n` exactly) is lonely *only*
at resonant times `t = a/(nm)`. Proof is two lines: liveness needs `min_v‖vt‖ ≥ 1/n`, tightness caps
it at `≤ 1/n`, so the min equals `1/n`, and a rational `‖v·p/q‖ = 1/n` forces `n | q`. No structure of
the family is used — not difference-closure, not being an AP, nothing. It is a property of the *number*
`1/n` against the rational grid.

**Fact B (the tight families are census twins).** At `n = 14` the primitive tight locus is
`{AP, GW}` with `AP = {1..13}` and `GW = {1..11,13,24}`. Their live census is *identical*:
`liveCount(q) = 6·[14|q]`, with the *same* six witnesses `{o·(q/14) : o ∈ (ℤ/14)*}` at every
resonance, zero everywhere else. Two genuinely different sets — one an AP, one not even
difference-closed — produce the exact same loneliness data at every modulus.

Put together, A and B say something I had not appreciated: **the loneliness census is a coarse,
degenerate invariant on the tight locus.** It collapses the tight families to a single orbit
`(ℤ/n)* ↷ resonance`. You cannot recover *which* tight family you have from when-and-how it is lonely.

## 2. Where the difference actually is

Drop below the safe band. Form the full depth histogram
`h(q) = multiset{ min_v (q·‖v·p/q‖) : p ∈ (0,q) }` — the entire coverage spectrum, not just its
top slice `≥ q/14`. Now AP and GW *split*: their histograms differ at `q = 70, 98, 112, 126, 140, …`
(first at `q = 70`, depths `2,3,4` carry counts `(26,6,6)` for AP versus `(24,4,10)` for GW). The
live slice is byte-identical; every bit that tells the two families apart sits in the sub-threshold
depths.

So the object has a layered structure the lenses had not named crisply:

- **At the threshold** (the live points, the safe band, `μ₀`, `liveCount`, `B5` at the resonant
  modulus, the deep count *at the band*): the tight families are **one object**. This layer is flat.
- **Below the threshold** (the multiplicity spectrum `μ_k`, the whole depth histogram, the
  additive-relation lattice): the tight families **separate**. This layer carries the rigidity.

## 3. Why this is the thing we were missing

The hard core of LRC(14) (THM-995's equality horn) is a *rigidity* statement: **tight ⟹ `{AP, GW}`**
(up to dilation) — the classical "tight ⟹ AP-like" characterization, LRC-hard. Every witness-side
route the fleet has built — the B5 funnel, the live-floor bridge (THM-984), the deep-count race
(THM-987), the modulus supply (THM-979) — reads the family *at the threshold*: it asks "is there a
live multiplier, and how many." Fact B is a theorem that **these routes are constitutionally blind to
rigidity**: they return identical values on AP and GW, so no threshold-level functional can be the
predicate that separates the tight locus from its would-be impostors. The B5/live machinery is a
*supply* engine (it manufactures a witness once `μ₀ > 0`), never a *detector* of tightness — which is
exactly HYP-7245's scope finding (the funnel needs a dissociation hypothesis; the tight families slip
through it), now with a structural reason: at the threshold there is nothing to detect, because the
tight families are the same point.

This is the discrete, witness-side shadow of LEM-020 (klein): pair data — moment order ≤ 2 — can
*never* certify covering; the minimal certifying order is 3. Same phenomenon, read on the census:
order ≤ 2 = the threshold layer = flat; order ≥ 3 = the sub-threshold spectrum = where rigidity lives.
Two independent lenses (the analytic coverage-moment lens and the arithmetic live-census lens) pin the
*same wall at the same place*. That convergence is the signal that the layer boundary is real
structure, not an artifact of either method.

## 4. The tournament reading (cut is flat, cycle carries the rigidity)

The resonance witnesses are the unit group `(ℤ/14)* ≅ (ℤ/7)* = μ₆`, the multiplicative six-cycle,
the automorphism carrier of the self-complementary heptagon tournament `R₇` (`z⁷ = −1`). In the
cut⊕cycle grammar the repo keeps rediscovering (A4 observer `+1` baseline; C4 "LRC cycle space
reinforces where the tournament cancels"; the section-vector movie of THM-890): the threshold layer
is the **cut / transitive** component — the mean, the `+1` the observer always keeps, carried by the
unit group and identical for both tight families. The sub-threshold layer is the **cycle** component
— the additive-relation lattice, where GW's tell is its shell partner `3 + 24 = 27 = 3³` (the
signed-LRC cut, D11) that the AP simply does not have. AP is shell-partner-free; GW is not; and that
support-3 relation is invisible to the pair/threshold data by exactly the OCF-side mechanism that
makes even cycles cancel and odd cycles survive. **The tournament face was telling us all along that
rigidity is a cycle-space property, and the census now shows the same thing from the runner side: the
threshold is the cut, and the cut is flat.**

## 5. A concrete next move (quantitative, toward the margin floor)

Fact A is the `μ₀ = 0` boundary as a discrete statement. The quantitative version is a rate: for a
*near*-tight family, how far off resonance must you go before the first live multiplier appears? The
first non-resonant `q` with `liveCount(q) > 0` should scale like `1/(M − 1/14)` — the discrete shadow
of boxeph's Lipschitz converter `μ₀ ≥ 2(M − 1/14)/v_max`. If that off-resonance-onset scale can be
bounded *below* by a family-uniform quantity on trapped cores, it is a witness-side route to the
double-threshold margin floor `trapped ⟹ M ≥ 1/7` (HYP-7300) — the second, softer half of the
nucleus. The census, having given up on detecting rigidity at the threshold, may still measure it as a
*rate of departure* from the threshold. That is the honest next question this finding leaves open.

→ THM-996, THM-991, THM-995, LEM-020, HYP-7245, HYP-7300, HYP-7305, MISTAKE-100.
