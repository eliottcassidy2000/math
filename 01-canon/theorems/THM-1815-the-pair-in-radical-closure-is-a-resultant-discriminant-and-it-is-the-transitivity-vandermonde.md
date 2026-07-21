---
id: THM-1815
title: "THE GMC(2) PAIR-IN-RADICAL CLOSURE IS A RESULTANT / MOMENT-MATRIX DETERMINANT NON-VANISHING — and that determinant is the transitivity Vandermonde (klein THM-1805), so it is opus THM-1710's resultant on the moment functional, confirming THM-1720's one-conjecture via discriminants. For a pair-straddle whose busier charge carries r terms at DISTINCT radial degrees a_1<…<a_r (coefficients β_i; the opposite charge one term α), the tower E[P^{j·m₀}] = C(j·m₀,j)·αʲ·L(Qʲ) (Q=V^p q^p, THM-1760) gives r equations whose ELIMINATION resultant is a NONZERO integer times a power of the top β — verified: r=2 gives Res(E[P^{m₀}],E[P^{2m₀}]) = 8·b₁² (radial 0,1), 504·b₁² (0,2), 192·b₁² for P=αZ+b₀W+b₁ZW², forcing all β=0 given α≠0. So the pair-atom form is in radical(I): PAIR-IN-RADICAL IS A DISCRIMINANT NON-VANISHING, closed on every tested pattern. (B) THE DETERMINANT IS THE VANDERMONDE of the radial degrees — the moment matrix [(a_i+k)!] has det 2,12,24,288 for degrees {0,1},{0,2},{0,1,2},{0,1,3} (nonzero ⟺ distinct degrees) — which by klein THM-1805 equals Σ_T sgn(T)·x^{score(T)} with the TRANSITIVE tournaments surviving. So the discriminant that forces the coefficients IS the transitivity structure of the charge lattice: the in/transitivity pivot (THM-1780) made a determinant. (C) This is EXACTLY opus THM-1710's Res(CT(m₀),CT(2m₀))≠0, the same object as TNC — so GMC(2)'s uniform closure is TNC's uniform resultant-non-vanishing (the multinomial-ratio conjecture), THM-1720's one-conjecture now realized as one discriminant."
status: >
  The pair-with-multiplicity closure is VERIFIED-EXACT on 5 towers (r=2,3): the elimination
  resultant is a nonzero integer times a power of the top β, forcing all β=0. The moment-matrix
  (Vandermonde) determinant is nonzero for distinct radial degrees, computed exactly. The
  identification with THM-1805's tournament-sum Vandermonde and THM-1710's TNC resultant is
  structural (the moment matrix reduces to the radial-degree Vandermonde; THM-1805 expands that
  Vandermonde as the tournament sum). The UNIFORM non-vanishing for all patterns is opus
  THM-1710's remaining multinomial-ratio step — SHARED with TNC, not closed here.
source: mac-mini-2026-07-20-S155 (owner: work the general case; think discriminants, determinants,
  and similar concepts)
depends_on:
  - THM-1780  # pair-reduction: the covering reduces to pairs-in-radical
  - THM-1760  # the tower identity E[P^{j m0}] = C alpha^j L(Q^j)
related:
  - THM-1805  # Vandermonde = signed tournament sum, transitivity survives (klein)
  - THM-1710  # resultant replaces cyclotomic for TNC (opus) -- the same resultant
  - THM-1720  # GMC(2) and TNC are one Nullstellensatz -- now one discriminant
  - THM-467   # negative OCF = reciprocal determinant quotient (the repo's determinant lens)
  - HYP-8590  # the localisation lemma / pair-in-radical
---

# THM-1815 — the pair-in-radical closure is a resultant, and it is the transitivity Vandermonde

## (A) Pair-in-radical is a discriminant non-vanishing

THM-1780 reduced GMC(2) to: every pair-straddle atom form lies in `radical(moment ideal)`. For a
pair whose busier charge carries `r` terms at distinct radial degrees `a_1<…<a_r` (coefficients
`β_i`; the opposite charge one term `α`), the tower `E[P^{j·m₀}] = C(j·m₀,j)·αʲ·L(Qʲ)` (THM-1760)
gives `r` equations in the `β_i`. **Their elimination resultant is a nonzero integer times a
power of the top `β`:**

| tower | resultant |
|---|---|
| `p=1`, degrees `{0,1}` | `8·b₁²` |
| `p=1`, degrees `{0,2}` | `504·b₁²` |
| `p=2`, degrees `{0,1}` | `1040449536·b₁⁸` |
| `p=1`, degrees `{0,1,2}` | (nonzero)·`b₂¹²` |

Concretely for `P = αZ + b₀W + b₁ZW²`: `E[P²] = 2α(b₀+2b₁)`, `E[P⁴] = 12α²(b₀²+6b₀b₁+12b₁²)`, and
`Res_{b₀}(E[P²], E[P⁴]) = 192·b₁²` — forcing `b₁=0`, then `b₀=0`. So the pair-atom form is in
`radical(I)`:

> **Pair-in-radical is the non-vanishing of a resultant / moment-matrix determinant** — a
> discriminant of the radial-degree tower. Verified on every tested pattern.

## (B) The determinant is the transitivity Vandermonde

The tower's leading system is linear in the radial-degree "frequencies" `L(V^{a_i+k}) = (a_i+k)!`
— a **moment matrix `[(a_i+k)!]`**, of Vandermonde type in the `a_i`:

| radial degrees | det |
|---|---|
| `{0,1}` | 2 |
| `{0,2}` | 12 |
| `{0,1,2}` | 24 |
| `{0,1,3}` | 288 |

Nonzero exactly when the degrees are **distinct** — i.e. genuine multiplicity. And by **klein
THM-1805**, the Vandermonde `∏_{i<j}(x_i−x_j) = Σ_T sgn(T)·x^{score(T)}` with the **transitive
tournaments surviving** (intransitive cancel by 3-cycle reversal). So:

> **The discriminant that forces the coefficients IS the transitivity structure of the charge
> lattice.** The in/transitivity pivot (THM-1780: one-sided = transitive, two-sided = a cycle)
> is, at the level of the closure mechanism, a Vandermonde determinant — the classical
> discriminant of invariant theory. `tournaments = in/transitivity` and `binary-form
> discriminant` are the same object, exactly as THM-1805 says.

## (C) One discriminant with TNC

`Res(E[P^{m₀}], E[P^{2m₀}]) ≠ 0` is precisely **opus THM-1710's** `Res(CT(m₀), CT(2m₀)) ≠ 0` — the
same resultant, now on the Laplace moment functional instead of the constant-term functional. So
THM-1720's "GMC(2) and TNC are one Nullstellensatz" sharpens to **one discriminant**: both
uniform closures are the non-vanishing of the two-representation resultant, which opus's
remaining **multinomial-ratio** step (THM-1710) would prove for all patterns. Proving it closes
both GMC(2) and TNC.

## Honest scope

- **The pair-with-multiplicity closure is verified on 5 towers, not proved uniformly.** Each
  resultant is exactly computed and nonzero, so the pair-in-radical is *closed on those
  patterns*; the uniform non-vanishing is THM-1710's open multinomial-ratio conjecture, shared
  with TNC.
- **The `A_r`/Vandermonde ⟺ transitivity identification is klein THM-1805's**, cited; what is
  new here is that the *GMC(2) closure mechanism* is that determinant. The reduction of the
  moment matrix `[(a_i+k)!]` to the radial-degree Vandermonde is via the factorial-ratio
  factorization (standard), stated not re-derived in full.
- This closes the *multiplicity* aspect of pair-in-radical; the *cross-atom isolation* (that a
  pair's tower is isolated at its return level once lower atoms are killed) is still the renewal
  induction (THM-1770 A), not subsumed by the resultant.
- Determinant lens cross-links (THM-467 reciprocal-determinant, the skew-Seidel determinant) are
  related framing, not used in the proof.

*Artifacts:* `04-computation/gmc2_resultant_discriminant_closure_macmini_S155.py` (+out).
*Credits:* THM-1780/1760 (the reduction and tower), klein THM-1805 (Vandermonde = tournament
sum), opus THM-1710 (the resultant), THM-1720 (the one-conjecture unification).
