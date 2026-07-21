---
id: THM-1840
title: "THE BOTH-SIGNS SINGLE-CHARACTER NULLCONE NON-VANISHING IS PROVED FUNCTIONAL-AGNOSTICALLY, and it is the CYCLOTOMIC-CLEAN base case of the whole moment/covering program. (A) For a single-character both-signs element P = a·χ₊ + b·χ₋ (χ₊ charge +p, χ₋ charge −n, gcd(p,n)=1, a,b≠0), its ONLY atom is n copies of +p and p copies of −n (size m₀ = p+n), so for ANY functional L with return weight W(m₀) ≠ 0, E_L[P^{m₀}] = C(m₀; n,p)·aⁿ·bᵖ·W(m₀) — a SINGLE non-degenerate term, nonzero because no other balanced tuple of size m₀ exists to cancel it. Verified: (2,3) → 10·W·a³b², (3,5) → 56·W·a⁵b³, generic weight W. This is FUNCTIONAL-AGNOSTIC (factorial/Gaussian, sinc/LRC with W(m₀)≠0, or any weight) and CYCLOTOMIC: a single minimal vanishing sum of the characters w^p, w^{−n} is one term; a cyclotomic cancellation needs ≥ 2 relations. So the single character is EXACTLY the regime where the cyclotomic single-shot is clean — the multi-character case is where it fails (opus THM-1710's non-cyclotomic tuned points) and the resultant/discriminant (THM-1815) takes over. (B) That failure begins with TWO 3-CYCLE ATOMS (minimal 3-charge vanishing sums): in the symmetric charge model they first coexist at n ≥ 6 (counts 0,0,0,2,2,4,4,8,8,12,12,18 for n=3..14); the owner's n≥13 is a model-specific threshold (LRC's 13 speeds), noted not re-derived. (C) THE BLUE PARITY (1 at odd n, 0 at even n) = (1−(−1)ⁿ)/2 is the nontrivial ℤ/2 cyclotomic character, and it is the SAME parity as the bicycle-space triviality / skew-Seidel forced-zero eigenvalue (THM-1440 C) and the complement S↦−S self-conjugate fixed direction — the SC/complement-parity lore in cyclotomic (degree-2) form."
status: >
  (A) PROVED — the single-character return is one non-degenerate term (elementary, general in
  the functional), verified symbolically for factorial and generic weights on 4 pairs. This is
  the functional-agnostic base case; it is not new as a fact (two-sided binomials are never in
  the nullcone) but the functional-agnostic + cyclotomic characterization is the contribution.
  (B) The two-3-cycle-atom threshold is COMPUTED in a symmetric charge model (n≥6, exact
  counts); the n≥13 claim is model-dependent and reported, not verified in the owner's model.
  (C) The parity identifications are exact (bicycle dim / forced-zero eigenvalue are THM-1440,
  cited); the 'blue parity = this parity' is a framing linking the tiling-blue lore to the
  ℤ/2 cyclotomic character.
source: mac-mini-2026-07-20-S158 (owner: prove the abstract both-signs single-character nullcone
  non-vanishing by the functional-agnostic method; n≥13 admits two 3-cycle atoms; the blue parity
  is a clean instance of the SC/complement-parity lore; think cyclotomic)
depends_on:
  - THM-1780  # the pair-reduction; this is its single-character base case
  - THM-1815  # the resultant/discriminant that takes over at multi-character
related:
  - THM-1835  # boxeph: SAME owner directive, GIT/atom-stratification + complement-reversal side
  - THM-1710  # opus: cyclotomic single-shot REFUTED for multi-character (trinomials)
  - THM-1440  # bicycle-space / forced-zero-eigenvalue parity = the blue parity
  - THM-1805  # Vandermonde = tournament sum (the 3-cycle = intransitivity)
  - HYP-8625  # the functional-agnostic barrier (single-character is where it works)
---

# THM-1840 — the single-character nullcone non-vanishing is cyclotomic-clean

> **Concurrent-work note.** boxeph's **THM-1835** answers the *same* owner directive
> (functional-agnostic single-character non-vanishing; n≥13 two 3-cycle atoms; blue parity as
> SC/complement lore) from the **GIT / tournament** side — the 0-unstable atom stratification, a
> correction to THM-1830 (4-atom forms enter at n=9), and the complement-reversal parity lemma.
> This file (mac-mini, ceded 1835 → 1840, boxeph pushed the stub first) is the **complementary
> analytic side**: the functional-agnostic *moment* proof (one term, any weight) and the explicit
> **cyclotomic** reading (single relation clean vs ≥2 relations non-cyclotomic, opus THM-1710).
> Same three-part answer, two lenses; cross-read them together.

## (A) The functional-agnostic single-character non-vanishing

A **single-character both-signs** element is a binomial `P = a·χ₊ + b·χ₋`, `χ₊` of charge `+p`,
`χ₋` of charge `−n`, `gcd(p,n)=1`, `a,b ≠ 0`. Its only atom (minimal balanced tuple) is `n`
copies of `+p` and `p` copies of `−n` (size `m₀ = p+n`, since `n·p = p·n`).

> **For ANY functional `L` with return weight `W(m₀) ≠ 0`,
> `E_L[P^{m₀}] = C(m₀; n,p)·aⁿ·bᵖ·W(m₀) ≠ 0`.**
>
> *Proof.* At size `m₀` there is exactly one balanced tuple type (`n` of `+p`, `p` of `−n`), so
> `E_L[P^{m₀}]` is a **single term** — nothing else of size `m₀` exists to cancel it. The
> multinomial coefficient `C(m₀; n,p)` is a positive integer, `aⁿbᵖ ≠ 0` (both signs), and
> `W(m₀) ≠ 0`. ∎

Verified symbolically: `(2,3) → 10·W·a³b²`, `(1,4) → 5·W·a⁴b`, `(3,5) → 56·W·a⁵b³` (generic
weight `W`); factorial weights give the same monomial times `(m₀-block)!`.

This is **functional-agnostic** — it holds for the Gaussian/factorial weights (GMC/NC2), for the
LRC/sinc weights whenever `W(m₀) ≠ 0`, and for any weight system — because the single term never
cancels. And it is **cyclotomic**: the atom is a minimal vanishing sum of the characters `w^p`,
`w^{−n}` (roots of unity on `|w|=1`); a single minimal relation is one non-degenerate term, and
cyclotomic cancellation requires `≥ 2` relations.

> **The single character is exactly the regime where the cyclotomic single-shot is clean.** It
> is the base case of the pair-reduction (THM-1780): every pair-straddle, restricted to its two
> characters, is a single-character binomial and so non-vanishes here.

## (B) Where it stops: two 3-cycle atoms

The clean cyclotomic non-vanishing fails as soon as there are **≥ 2 three-cycle atoms** (minimal
3-charge vanishing sums, e.g. `{+2,+3,−5}`): their returns can interfere, and opus **THM-1710**
showed the cancellation points are **non-cyclotomic** (the tuned modulus is a ratio of
multinomials, generically off the unit circle). There the cyclotomic single-shot is dead and the
**resultant / discriminant** (THM-1815) takes over.

Computed counts of 3-cycle atoms in the symmetric charge model `{−h,…,h}∖{0}`:

| n | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| #3-atoms | 0 | 0 | 0 | 2 | 2 | 4 | 4 | 8 | 8 | 12 | 12 | 18 |

So two 3-cycle atoms first coexist at **n ≥ 6** in this model. **The owner's `n ≥ 13`** is a
different, model-specific threshold (the LRC 13-speed set), reported and deferred, not re-derived
here — see Honest scope. A 3-cycle atom is a directed 3-cycle = intransitivity (THM-1805), so
"two 3-cycle atoms" is "two independent intransitive triples," the first genuinely multi-character
obstruction.

## (C) The blue parity is the ℤ/2 cyclotomic character

The **blue parity** `b(n) = 1` (odd `n`), `0` (even `n`) equals `(1−(−1)ⁿ)/2` — the nontrivial
character of `ℤ/2`, the degree-2 cyclotomic character. It coincides with:

| `n` | `b(n)` | bicycle dim (THM-1440 C) | skew-Seidel forced-zero (odd n) |
|---|---|---|---|
| odd | 1 | 0 | yes (corank ≥ 1) |
| even | 0 | `n−2` | no |

So the tiling-blue parity is the **same parity** as the bicycle-space triviality, the skew-Seidel
forced-zero eigenvalue, and the complement `S ↦ −S` self-conjugate fixed direction (THM-1440).
The SC/complement-parity lore is, at its degree-2 core, this one ℤ/2 cyclotomic character — the
generator of the whole odd/even axis the project runs on.

## The cyclotomic thread

The three pieces are one story told in the cyclotomic language:

- **single character** → one minimal relation → cyclotomic non-vanishing is clean (this theorem);
- **two 3-cycle atoms** → multiple relations → cyclotomic cancellation possible but non-cyclotomic
  in fact (opus THM-1710) → the discriminant/resultant (THM-1815) replaces it;
- **blue parity** → the ℤ/2 cyclotomic character → the SC/complement-parity odd/even axis.

## Honest scope

- **(A) is elementary as a fact** (a two-sided binomial's `m₀`-th moment is a nonzero monomial
  times the return weight); the contribution is the *functional-agnostic + cyclotomic framing*
  making it the clean base case of the pair-reduction and identifying exactly where it stops.
- **(B)'s `n≥13` is not matched by the computed model** (which gives `n≥6`); I report the honest
  computation and treat the `13` as a claim about a specific model (LRC's 13 speeds), not verified
  here. Do not cite `n≥13` from this file.
- **(C) is a framing** linking the tiling-blue parity to the ℤ/2 cyclotomic character and to the
  THM-1440 parities (cited); "blue parity = bicycle parity" is an identification of two 0/1
  sequences that agree, presented as lore-unification, not a new invariant.
- The functional-agnostic claim requires `W(m₀) ≠ 0`; for the LRC sinc weights this fails at
  special `δ`, which is exactly the S157 measure barrier — the single-character result is clean
  only where the weight does not vanish, consistent with that.

*Artifacts:* `04-computation/single_character_nullcone_cyclotomic_macmini_S158.py` (+out).
