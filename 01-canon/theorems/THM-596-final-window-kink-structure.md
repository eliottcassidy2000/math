# THM-596: Final-window kink structure — band arithmetic, the absorption mechanism, and the FALSITY of the universal K=0 lemma

**Status:** PROVED (parts i–iii) + explicit COUNTEREXAMPLE (part iv) + refined reduction (part v)
**Author:** mac-mini-2026-07-01-S95 (HYP-3851)
**Verification:** exact `Fraction` computations via `lonely_profile.py` (this session's inline runs; counterexample re-checkable in one call).
**Context:** klein-S88 (HYP-3841) observed the final window `[1/15, 1/14]` overtaking-free in all tests and proposed "K₂ = 0" as a candidate lemma making the tangent-ladder's last rung defect-free universally; the owner asked for a proof. The truth is sharper and two-sided.

## (i) Band arithmetic (proved)

Every overtaking kink (THM-592(iii)) at radius `ρ ∈ (1/15, 1/14)` (open window) comes from a pair `v < w` with crossing points `t* = m/(w−v)`; writing `g = gcd(m, w−v)`, `q* = (w−v)/g`, `d' = (v·m mod (w−v))/g`, the kink radius is `ρ = d'/q*` with
```
q* ∈ (14·d', 15·d'),   d' ≥ 2.
```
`d' = 1` is impossible in the open window (no integer in `(14, 15)`): first admissible bands `q* = 29 (d'=2)`, `q* ∈ {43,44} (d'=3)`, `{57,58,59} (d'=4)`, …  Kinks with `ρ = 1/14` exactly (the `d'=1` boundary) do not affect the ladder integral over the open window.

*Proof.* At a same-side endpoint crossing, `(w−v)·t* = j'−k ∈ ℤ` (subtract the two endpoint equations), so `t* = m/(w−v)`; `ρ` is the fractional part `v·m/(w−v) mod 1 = d/(w−v)`, and `g | d` since `g | m` and `g | (w−v)`; reduce. The window condition `1/15 < d'/q* < 1/14` is the band. ∎

## (ii) Exposure characterization (proved)

The kink at `t* = m'/q*` (reduced) is **exposed** (contributes to `K₂`) iff `t*` is `ρ`-lonely: the residue system `{u·m' mod q* : u ∈ S}` avoids the open band `(−d', d')`, with the crossing pair at `+d'` exactly. The mirror point `1 − t*` carries the mirrored kink (residues negate), so exposed defects come in pairs — total per crossing `2(1/v − 1/w)`.

## (iii) The absorption mechanism (proved) — why natural covering sets test clean

At the minimal band `q* = 29`: if the covering set's multiple of 14 is `28 ≡ −1 (mod 29)`, its residue at any crossing involving runner `1` (residue `m' = d'`) is exactly `−d'` — the closed-boundary case. Then THREE endpoints coincide at `(t*, ρ)` (the two crossing right-endpoints and 28's left endpoint), the lonely component **dies** at the event, and the would-be concave kink is absorbed into a convex death: net `K` contribution zero. More generally a runner at residue exactly `−d'` converts the overtaking into a component death. This is why every naturally-sampled covering set (klein-S88's 11/11; S94's 22/25 window census) shows `K₂ = 0`: the standard covering completions carry a `∓d'`-residue absorber.

## (iv) The universal K₂ = 0 lemma is FALSE (explicit counterexample)

```
S = {1, 2, 3, 5, 6, 7, 8, 9, 11, 12, 26, 30, 42}
```
is a covering 13-set (multiples of every `q ∈ {2..14}`: q=4→8,12; q=10→30; q=13→26; q=14→42; …) with **exposed concave kinks at `ρ = 2/29 ∈ (1/15, 1/14)`, total defect `K₂ = 29/15 ≈ 1.933`** (the `(1,30)` crossing at `t* = 2/29` plus its mirror; residues `·2 mod 29`: `{2,4,6,10,12,14,−13,−11,−7,−5,−6,−3,2}` — all clear of `(−2,2)`). Construction recipe: place the crossing (pair `(1, 30)`, `w−v = 29`), then choose the 14-multiple as `42 ≡ 13 (mod 29)` (residue `−3·2·…` clear) instead of `28` (which would absorb, by (iii)); check the remaining residues clear the band. Note `M(S) > 1/14` — the set is LRC-safe; it breaks only the *universality* of the K₂=0 reduction.

## (v) The refined reduction (what replaces "three numbers")

klein-S88's reduction (covering floor ⟺ `a − b/210 ≥ 0` with `a` = anchor mass at 1/15, `b` = slope cap, K₂ = 0) survives in the refined form: covering floor ⟺
```
a − (b + K₂^danger)/210 ≥ 0,
```
where `K₂^danger` = the exposed final-window defect **over the danger zone only** (sets with `M(S) ≤ 1/14 + δ`). Structure available to bound it: an exposed kink forces (a) a pair with `w − v ∈ {29g, 43g, 44g, 57g, …}`, (b) a `q*`-denominator `ρ`-lonely point with `ρ > 1/15` — a deep-witness event at a modulus in the sparse band list, which the danger-zone classification (forced-cover obstruction, HYP-3792/3789; witness CF-signature, definitions §witness) constrains, and (c) no `∓d'`-residue absorber among all 13 runners — 13 simultaneous non-membership conditions mod `q*`. Alternatively: anchoring the last rung at `2/29` instead of `1/15` empties the `q* = 29` band (the window `(2/29, 1/14)` admits only `q* ∈ {43, 57, 71, …}`-band kinks), pushing the counterexample family to higher moduli at the cost of a shorter rung — an anchor-vs-band tradeoff ladder.

## Cross-references

Discrete/continuous rigidity twins: THM-594(C). Ladder: THM-592(v) + S94 weighted sharpening. klein HYP-3841 (instrument confirmed; the candidate lemma refuted as universal, alive as danger-zone-restricted). Absorption = the "boundary case protects natural sets" phenomenon — same shape as the tight-set covering-patch domino (THM-593 interpretation).

-> THM-592, THM-593, THM-594, HYP-3841, HYP-3850, HYP-3851, OPEN-Q-108.
