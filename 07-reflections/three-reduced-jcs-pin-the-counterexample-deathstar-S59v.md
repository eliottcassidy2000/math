# Three reduced Jacobian conjectures, proven in one day, jointly pin the counterexample

**death-star-2026-07-20-S59v** (HYP-8205; owner: long high-leverage session, pull
often, invent a reduced JC that holds). What actually happened: the owner's prompt
went fleet-wide and **three different reduced JCs were each proven by a different
agent within hours**, and my own attempt (the weight-sign one) converged with
mac-mini's THM-1370 while I was computing it (MISTAKE-199, 4th recurrence — credited,
not re-filed). This reflection is the non-competing synthesis: the three theorems are
three independent invariants, and the owner counterexample is the **unique minimal
object that escapes all three at once.** Everything here credits the primary authors;
my only original content is the intersection reading + the independent verification.

## The three reduced JCs (each a theorem, each by a different agent, all 2026-07-20)

| # | Reduced JC that HOLDS | invariant it controls | author |
|---|---|---|---|
| **SIGN** | **elliptic / definite-weight JC**: a ℂ*-equivariant Keller map whose weights are all one sign is an automorphism, every dimension (properness: F⁻¹(0) is ℂ*-invariant + étale = {0} ⟹ finite étale over simply-connected ℂⁿ ⟹ degree 1). | the weight **signature** (definite vs indefinite) | **mac-mini THM-1370** (I independently verified it; = my intended weight-sign result) |
| **PARITY** | **equivariant fixed-locus JC**: for σ,τ-equivariant F, the restriction F|Fix(σ) is a constant-Jacobian map in lower dimension; when dim Fix ≤ 1 that is JC₁ (true), so F is injective on Fix(σ) unconditionally. | the weight **parity** (Fix(σ) = Fix of λ=−1 = the even-weight coords) | **opus THM-1350** |
| **DEGREE** | **JC≤2**: every polynomial local diffeomorphism of ℂⁿ of geometric degree ≤ 2 is injective, all n (Smith rule; the conjecture fails first and only at geometric degree 3). | the **geometric degree** of the map | **klein-S333** |

Three theorems, three structural features of a polynomial map — its ℂ*-weight
*sign*, its weight *parity*, and its geometric *degree*. Each says: if this feature
is "small," the map is forced injective.

## The intersection: why the counterexample is exactly what it is

A counterexample must violate all three — it must have **indefinite weights**
(else SIGN makes it invertible), a **fiber of odd size ≥ 3** (PARITY: the σ-action
on a τ-fixed fiber leaves an odd number fixed, and F|Fix gives exactly one, so
|fiber| is odd; a genuine collision needs ≥ 3), and **geometric degree ≥ 3** (else
DEGREE makes it injective). The owner's map F is the object that meets each bound with
equality:

- weights **(1, −1, −2)** — indefinite, and minimal: dimension 2 is excluded by the
  equivariant JC₂ theorem (THM-1345), so dim 3 is the floor, and (1,−1,−2) is the
  forced grading (mac-mini's "sharp/unique");
- **fiber = 1 + 2 = 3** — the minimal odd collision (one σ-fixed point on the z-axis +
  one doubled free orbit λ ↦ λ²), exactly opus's forced triple;
- **geometric degree = 3** — the minimal viable degree, exactly klein's / the engine
  trichotomy's cuspidal-cubic floor.

So the three reduced JCs are not three separate facts about the counterexample; they
are three faces of its minimality. **F is the unique smallest object that is
indefinite, odd-fibered, and cubic at once** — and each "smallest" is a distinct
proven reduced JC. The counterexample sits precisely at the triple point where all
three reductions fail together, and it fails each of them by the minimum margin.

## Why the three invariants are genuinely independent (not the same theorem thrice)

They decompose the weight vector by different features and one non-weight feature:
- SIGN partitions the coordinates by sign(w_i) — a properness / attracting-fixed-point
  property (definite ⟹ the origin is the whole zero fiber).
- PARITY partitions by w_i mod 2 — a σ = (λ=−1)-involution property (even-weight coords
  are σ-fixed). (1,−1,−2) has sign-signature (1,2,0) and parity-signature (odd,odd,even):
  the two partitions differ, so SIGN and PARITY see different subspaces of the same map.
- DEGREE is not a weight property at all — it is the covering degree, controlled by the
  Smith rule / Campbell, and it is what the equivariant frame does *not* see.

That the same minimal object saturates all three independent bounds is the content:
the counterexample is over-determined, pinned from three directions, which is exactly
why it is rigid (THM-1305), sporadic (S59n), and unique-in-its-box (kps-S128c98).

## What this hands the fleet

- A **clean statement of the frontier's shape**: to find a *second* essentially new
  counterexample, you must move off at least one of the three bounds — a genuinely
  new indefinite signature (the dim-3 landscape beyond (1,−1,−2); (1,−1,−3) is empty,
  S59n), a fiber of size 5 or 7 (the next odd rungs, opus's forcing), or geometric
  degree 4/5 (the realizability program; A₄/S₄ needs z-degree ≥ 2, the Smith rule).
  The reduced-JC trio *is* the classification's coordinate system.
- The honest meta-lesson (MISTAKE-199, now 4×): on a hot fleet-wide prompt, three
  agents proving three complementary theorems in a day is the system working — the
  right individual move, once you see the collision, is to **synthesize the parallel
  results into the joint picture**, which no single one of them states. That is the
  only thing this session adds, and it is worth adding.

## Cross-links

mac-mini-S123 THM-1370 (elliptic/definite-weight JC) · opus-S399 THM-1350
(equivariant fixed-locus JC) · klein-S333 (JC≤2, the dihedral dictionary) ·
THM-1345 (equivariant JC₂, the dim-2 floor) · THM-1305 (rigidity) · THM-1300
(the counterexample) · S59n (the (1,−1,−k) landscape) · the PROBLEM-LEDGER (§A5/A8).
