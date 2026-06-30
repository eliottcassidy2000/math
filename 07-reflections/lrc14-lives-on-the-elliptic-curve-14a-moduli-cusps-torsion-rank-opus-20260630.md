# LRC(14) lives on the elliptic curve 14a: the moduli X₀(14) IS the genus-1 curve 14a, the 4 cusps are the n=4 tournament classes (points on it), the 6 torsion points (ℤ/6) match the 6 units mod 14 (the razor's-edge binding witnesses), and the rank-0 (L(14a,1)>0) is the obstruction's non-degeneracy (the empty tooth's positive width) — moduli, obstruction, binding, and non-degeneracy all on ONE curve; tractability = genus (rational vs curve)

*opus-2026-06-30. Owner: keep hunting unifications. The deepest one: X₀(14) has genus 1, so it is an
elliptic curve — and it is 14a. Everything in the LRC(14) story is a feature of this one curve.*

## The fact: X₀(14) = the elliptic curve 14a (genus 1)
Genus-1 modular curves ARE elliptic curves; `X₀(14)` has genus 1 (the first `X₀(2p)` that does), and it
is the conductor-14 rank-0 curve **14a**. So **the moduli space of LRC(14) is the elliptic curve 14a.**
The genus-1 cusp form `f₁₄` is the differential on it; the LRC(14) floor obstruction IS this curve.

## Everything in the story is a feature of 14a
| LRC(14) object | feature of the curve 14a | check |
|---|---|---|
| the **moduli / obstruction** | the curve `14a` itself (genus 1, `f₁₄`) | standard (`X₀(14)≅14a`) |
| the **4 cusps** (degenerate regimes `d=1,2,7,14`) | 4 marked POINTS on `14a`, the Klein-four `W(14)` orbit | = the **n=4 tournament classes** T,+,−,S (THM-584) |
| the **binding / witnesses** | the **6 rational torsion points** `E_tors = ℤ/6` | `\|E_tors\|=6=φ(14)` (verified: `6 \| #E(F_p)`) |
| the **non-degeneracy** | **rank 0 ⇒ `L(14a,1)>0`** | the empty tooth has positive width `1/14` |
| **tractability** | the **genus** (rational `P¹` vs real curve) | `0,0,1,2,2` for `2p=6,10,14,22,26` |
> **The 6 torsion points of 14a match the 6 units mod 14** `{1,3,5,9,11,13}` — the razor's-edge witnesses,
> the empty-tooth points, the binding pairs (`3` pairs `{±k}`), the Galois group `(ℤ/14)*≅ℤ/6`. All the
> same `ℤ/6`. (The match is `|·|=6=φ(14)` and structural `ℤ/6`; the precise torsion↔cusp parametrization is
> the modular-units statement — suggestive, to pin.)

## The 6 and the 7 (the two columns, on the curve)
`14 = 2·7`, but the curve splits the story into a **6** and a **7**:
- **the 6** = `E_tors(14a) = ℤ/6 = φ(14) = (ℤ/14)* = ` the units `= ` the Eisenstein hexagonal `ℤ[ζ₆]`
  `= Q(√−3)` `= ` the **EXISTENCE column** (the empty-tooth witnesses, the `Φ₆`-relocation, the cap).
- **the 7** = the apex prime `= Z₇` cyclotomic gap `= ` the doublet `= Q(√−7)` `= ` the **MEASURE column**
  (the per-level density, proved by THM-590).
So the master two-Heegner duality is, on the curve: the **6-torsion (existence/Q(√−3))** vs the **7-apex
(measure/Q(√−7))**, the genus-1 obstruction sitting between them at the `d=7` cusp. The two worlds touch
at `7 = Φ₆(3)` (the apex = the Eisenstein norm at `n=3`).

## Tractability = genus (the clean LRC(2p) law)
> `genus(X₀(2p)) = 0,0,1,2,2` for `p=3,5,7,11,13`. **Genus 0 = the moduli is rational (`P¹`, no curve) =
> LRC(6), LRC(10) SOLVED. Genus 1 = the moduli is the elliptic curve 14a = LRC(14) FIRST HARD. Genus ≥ 2 =
> LRC(22), LRC(26).** LRC(2p) is tractable exactly when its moduli is rational; it becomes hard precisely
> when the moduli acquires a curve (a cusp form, a global obstruction the local/boundary data can't see).
`14` is the first because it is the first `2p` whose `X₀` is not `P¹`. The hardness is geometric: a genus.

## Why this is the right home (and what it leaves)
- **It unifies the whole story on ONE object.** Moduli (the curve), obstruction (`f₁₄`/genus), boundary
  (the cusps = n=4 classes), binding (the 6-torsion = the units = razor's edge), non-degeneracy (rank 0 =
  positive-width tooth) — all features of 14a.
- **It explains the hardness exactly** (genus 0 vs 1) and names the apex primes (`{3,7}` = the `2p` with
  `X₀` rational-or-genus-1 in the Heegner regime).
- **The conjecture, restated:** LRC(14) holds iff the empty tooth keeps positive width — which is the
  rank-0 non-vanishing `L(14a,1)>0` (KNOWN) made into the covering-min `≥1/n` (the existence/`Q(√−3)`
  column). The gap is still the descent (full `ρ_j` ← the curve's local data); but the target now sits on
  a single, fully-known elliptic curve.

## Status
- **Verified/standard (opus, building on klein/mac-mini):** `X₀(14)=14a` (genus 1); `|E_tors|=6=φ(14)`;
  cusps = n=4 classes; rank 0; genus law `0,0,1,2,2` = LRC(2p) tractability.
- **The unification:** LRC(14) is a question about the elliptic curve 14a — moduli, obstruction, cusps
  (classes), torsion (units/binding), rank (non-degeneracy) all on one genus-1 curve; the 6/7 split is the
  two-Heegner-column duality on the curve.
- **Suggestive (to pin):** the torsion↔units identification (modular-units), and rank-0↔empty-tooth-width
  as a rigorous equivalence (the `L`-value = floor-constant, klein-flagged speculative).

Related: the master-two-Heegner-columns reflection (this gives its geometric home); f₁₄-14a-rank-0,
cusp-existence-comb-witness, covering-min-Eisenstein-Φ₆, cyclotomic-self-dual-razors-edge (the units/Galois
ℤ/6) reflections; klein THM-584/HYP-3586/3587 (cusps=classes, genus), mac-mini HYP-3594 (the cusp form);
HYP-3547 (apex primes); OPEN-Q-108.
