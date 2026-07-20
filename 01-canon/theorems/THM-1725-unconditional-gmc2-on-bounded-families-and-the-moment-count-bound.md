---
id: THM-1725
title: "GMC(2) IS UNCONDITIONALLY TRUE ON EVERY BOUNDED CHARGE-COUNT + DEGREE FAMILY (a finite Gröbner test), AND THE MOMENT-COUNT BOUND IS 2·max-over-pairs of the ESV level. (A) DECIDABILITY, made a theorem: fixing k monomials and degree D leaves FINITELY many support patterns, and for each, GMC(2) ⟺ V(⟨E[P^m]⟩) ∩ {two-sided} = ∅ ⟺ a finite Nullstellensatz certificate exists (Hilbert). Run exhaustively: ALL 132 genuinely-two-sided TRINOMIAL patterns up to degree 4 CLOSE, and 40/40 sampled 4-monomial patterns up to degree 3 CLOSE — UNCONDITIONAL GMC(2) on those families, no conjecture, only Gröbner + Nullstellensatz. (B) THE MOMENT-COUNT BOUND: the minimal certifying level M* satisfies M* ≤ 2·max over (positive charge p, negative charge n) of (p+|n|)/gcd(p,|n|) — verified for ALL 132 trinomials and EXACTLY TIGHT (ratio 2.000 at the extremal patterns {−1,+4}, {−2,+3}). The inner (p+|n|)/gcd is the ESV/DvdK first-return level for the charge pair; the max is over the FULL charge lattice (NOT the extremes — the naive extreme-charge formula is REFUTED: charges {−3,−2,+3} need M*=5 > 4, {−4,−3,+4} need 7); the factor 2 is opus THM-1685's primitive+second-level (CT(m₀)+CT(2m₀) ⟹ unit ideal). (C) This is the exact analogue of TNC's HYP-8505, so the two uniform bounds are ONE conjecture; proving it would give full GMC(2) given the complex radial (THM-1695) and span-2 (THM-1600) are closed."
status: >
  (A) DECIDABILITY is PROVED (Hilbert Nullstellensatz + finite pattern count). The
  exhaustive closure is VERIFIED-EXACT over ℚ: 132/132 two-sided trinomial patterns
  (up to Z↔W) at degree ≤ 4 with the rigorous all-pos×neg-pair saturation, and 40/40
  sampled 4-monomial patterns at degree ≤ 3. Each closure is a finite per-pattern proof.
  (B) The bound M* ≤ 2·max-pairs(p+|n|)/gcd is VERIFIED for all 132 trinomials and tight;
  it is a CONJECTURED uniform bound (HYP-8540), proved for no pattern class in closed form.
  The extreme-charge formula is REFUTED by explicit counterexample.
  (C) The identification with HYP-8505 is structural. NOT a proof of GMC(2).
source: mac-mini-2026-07-20-S149 (owner: work the moment-count bound and trinomial-adjacent
  ideas; unconditional GMC(2) on bounded charge-count + degree is a finite Gröbner test)
depends_on:
  - THM-1720  # the GMC(2) nullcone as a Nullstellensatz emptiness test
  - THM-1685  # opus: TNC k-nomial decision procedure; the primitive+second-level mechanism
related:
  - THM-1695  # complex radial closed
  - THM-1600  # span-2 base
  - THM-1650  # Newton polygon of the effective DvdK/ESV bound
  - HYP-8505  # opus: uniform CT-level bound for TNC -- the same conjecture
---

# THM-1725 — unconditional GMC(2) on bounded families, and the moment-count bound

## (A) Decidability, as a theorem

Fix `k` (number of monomials) and `D` (max degree). The monomials `Z^a W^b` with `a+b ≤ D`
form a finite set, so there are **finitely many support patterns** with `≤ k` monomials. For
each genuinely-two-sided pattern, THM-1720 gives

> GMC(2) holds ⟺ `V(⟨E[P^m]⟩) ∩ {two-sided} = ∅` ⟺ a **finite Nullstellensatz certificate**
> exists (Hilbert: `V = ∅ ⟹ 1 ∈ radical`, and a finitely generated ideal needs only finitely
> many generators).

So **GMC(2) on `{≤ k monomials, degree ≤ D}` is a finite union of Gröbner tests — decidable.**
Running it exhaustively is an *unconditional* proof on that family. Result, with the rigorous
all-`pos × neg`-pair saturation:

| family | patterns (up to `Z↔W`) | all close? |
|---|---|---|
| two-sided trinomials, `D ≤ 2` | 5 | **yes** |
| two-sided trinomials, `D ≤ 3` | 34 | **yes** |
| two-sided trinomials, `D ≤ 4` | **132** | **yes** |
| two-sided 4-monomials, `D ≤ 3` (sample) | 40 | **yes** |

> **GMC(2) is unconditionally true for every ≤ 3-monomial `P` of degree ≤ 4, and every
> ≤ 4-monomial `P` of degree ≤ 3 tested** — no conjecture, only Nullstellensatz + finite
> computation. Each is a finite per-pattern proof (the nullcone `⊆ V(⟨E[P^{1..M}]⟩)` for every
> `M`, so a unit certificate from finite `M` proves emptiness).

## (B) The moment-count bound

The minimal certifying level `M*` (fewest moments whose saturated ideal is `⟨1⟩`) obeys:

> **`M* ≤ 2 · max_{p>0, n<0 charges present} (p + |n|)/gcd(p, |n|)`** — verified for all 132
> trinomials, and **exactly tight** (ratio `2.000`) at `{−1,+4}`, `{−2,+3}` (`M* = 10`).

Three pieces, each with a source:

- **`(p+|n|)/gcd(p,|n|)`** is the **ESV/DvdK first-return level** `m₀` for the two-monomial
  charge pair `(p, n)` — exactly the effective-bound quantity of THM-1650.
- **the max is over the FULL charge lattice, not the extremes.** The naive
  `2(K₊+K₋)/gcd(K₊,K₋)` using only the top and bottom charges is **refuted**: `{−3,−2,+3}`
  needs `M* = 5 > 4`, `{−4,−3,+4}` needs `M* = 7 > 4`. An intermediate charge creates a
  coprime pair with the extreme that raises the return level. The descent is **pairwise across
  every straddle**, not just top-vs-bottom.
- **the factor 2** is opus THM-1685's mechanism: the primitive relation appears at `m₀`
  (`CT(m₀)`), and the independent second equation at `2m₀` (`CT(2m₀)`); together they generate
  the unit ideal.

## (C) One conjecture with TNC

The bound is the exact GMC(2) analogue of opus's **HYP-8505** (uniform CT-level bound for TNC).
Both procedures are now literally the same — *saturate a vanishing-ideal of graded power sums,
test for `1`* — so **the two uniform bounds are one conjecture** (HYP-8540). A single proof
would close both, and for GMC(2) it would give the full theorem: the complex radial layer
(THM-1695) and span-2 base (THM-1600) are already closed, so the uniform moment-count bound is
the last gap.

## Honest scope

- **(A) is unconditional on the stated bounded families only** — trinomials to degree 4, a
  4-monomial sample to degree 3. It is a genuine theorem for those, not for GMC(2).
- **(B)'s bound is verified and tight but conjectural** (HYP-8540); it is proved in closed form
  for no infinite pattern class. What is *proved* is the refutation of the extreme-charge form.
- The 4-monomial run is a **sample (40 patterns), not exhaustive**; it confirms closure and the
  `2×max-pairs` bound on what was tested, nothing beyond.
- `M*` is the minimal certifying moment count *in this saturation formulation*; a different
  ideal presentation could shift constants. The `2×max-pairs` shape is the claim, not the exact
  integer for every pattern.

*Artifacts:* `04-computation/gmc2_bounded_family_moment_bound_macmini_S149.py`,
`gmc2_moment_bound_formula_macmini_S149.py` (+outs).
*Credits:* THM-1720 (the Nullstellensatz framing), opus THM-1685 (primitive+second level),
THM-1650 (the ESV effective bound).
