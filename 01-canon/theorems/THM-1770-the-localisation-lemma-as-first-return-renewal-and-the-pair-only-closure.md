---
id: THM-1770
title: "CORRECTED: least-return words are primitive, but distinct atoms can cancel"
status: >
  CORRECTED by MISTAKE-211 and CASE-gmc2-first-return-cross-atom-cancellation.
  Part (A) is PROVED: a least-size balanced word is primitive. Parts (B)--(D)
  are RETRACTED: distinct coefficient monomials can cancel in one scalar moment
  equation, so atomwise vanishing, renewal isolation, and star/pair-only closure
  do not follow. The valid target is a multilevel radical/resultant theorem.
source: mac-mini-2026-07-20-S153 (owner: take the localisation lemma in the first-return /
  covering argument direction)
depends_on:
  - THM-1760  # single-straddle GMC(2) reduction; the single-atom base case
  - THM-1740  # r*m0, the multiplicity/Vandermonde tower
  - THM-1700  # bottom-up descent = the renewal order
related:
  - THM-1685  # opus: TNC Nullstellensatz / vanishing sums -- the general atom-covering core
  - THM-415   # prime-modulus no-collision = distinct atoms don't cancel
  - THM-1745  # max-over-straddles = first-return ordering
  - HYP-8590  # the localisation lemma this structures
  - THM-2014  # a valid entropy-versus-factorial separation theorem on an infinite slice
---

# THM-1770 — the localisation lemma as first-return / renewal

## The atoms

Write `P = Σ c_i Z^{a_i} W^{b_i}`, charge `k_i = a_i − b_i`. A balanced tuple has charge-sum 0;
`E[P^m]` sums over balanced `m`-tuples. An **atom** is a *minimal* balanced charge-multiset — a
minimal vanishing sum over the charge lattice (opus THM-1685; THM-415). Its size is its
**first-return level**. Atoms are the return events; the LRC/covering reading of THM-1745 is
that `M*` is governed by the last atom to return.

## (A) First-return isolation — PROVED

> **Lemma.** Let `m* = min` atom size. Every balanced tuple of size `m*` is a single atom.
> *Proof.* A composite balanced tuple splits into `≥ 2` balanced sub-tuples, each of size `≥ m*`
> by minimality, so total size `≥ 2m* > m*`. ∎

Hence **`E[P^{m*}] = Σ` over size-`m*` atoms of their coefficient forms, with no composite
terms** — the first return is primitive. Verified: `aZ²+bW+cW³` → `E[P³] = 6ab²`; the 3-charge
atom `(+2,+3,−5)` → `E[P³] = 720abc` (a minimal vanishing sum, not a pos-neg pair).

## (B) Distinct atoms can cancel — RETRACTED

Distinct atoms do give distinct coefficient monomials, but the nullcone condition is a
**scalar equation at one coefficient point**, not a polynomial identity in free variables.
Distinct monomials therefore need not vanish separately.

The exact star-support witness is

\[
P=aZ^6+bW^2+cW^{18},\qquad
\mathbb E[P^4]=4\cdot6!\,ab^3+4\cdot18!\,a^3c.
\]

The two summands are the two primitive length-four atoms. Taking

\[
c=-\frac{6!}{18!}\frac{b^3}{a^2}
\]

cancels them with (abc\ne0). THM-415's no-collision statement in its original
prime-modulus setting does not turn this scalar complex sum into separate equations.

## (C) Renewal is a candidate multilevel elimination, not an induction yet

First-return order remains a useful filtration, but (A) alone does not isolate an atom.
At successive multiples of a primitive return, new radial channels and mixed balanced
relations appear. A valid induction must eliminate earlier channels by explicit cumulants,
resultants, or a factorial-Hankel/Vandermonde determinant. That multilevel triangularity is
the open content of HYP-8765; it is not supplied by disjoint monomial support.

## (D) Pair-only/star closure — RETRACTED at this level of generality

If one had already proved every pair product (c_pc_n) belongs to the radical of the
**full** moment ideal, then their common zero set would indeed be the union of the two
one-sided coordinate subspaces. The missing premise is precisely radical membership.
The star witness above shows it cannot be obtained from the first-return equation alone.

The literal two-character/binomial case remains valid because there is only one balanced
tuple at each return. Named finite star supports may also be closed by later moments, but
there is no general star theorem here.

## Correct residual

For each positive-negative coefficient pair the needed statement is

\[
(c_pc_n)^N\in
\sqrt{\langle \mathbb E[P],\mathbb E[P^2],\ldots\rangle}
\quad\text{for some }N.
\]

Proving these inclusions for every pair would imply one-sidedness. The charge quotient alone
is insufficient: a monomial must retain both charge (a-b) and radial height (a+b), because
balanced words with the same charge carry different factorial weights. See MISTAKE-211 and
the resolved court case for the exact audit.

## Honest surviving scope

- **Proved:** the least-return balanced words are primitive; no composite balanced word occurs
  at that level.
- **Not proved:** atomwise vanishing, general renewal isolation, star closure, or a universal
  twice-return cutoff.
- **Still valid independently:** exact two-character/binomial results and any named support
  closed by a complete multilevel resultant calculation.
- Radial degrees are not passive multiplicities. They create independent channels and can
  increase detection depth even at fixed charge support.

*Artifacts:* `04-computation/gmc2_first_return_localisation_macmini_S153.py` (+out).
*Correction artifacts:* MISTAKE-211 and
`02-court/resolved/CASE-gmc2-first-return-cross-atom-cancellation.md`.
*Credits:* the primitive-word lemma is retained from mac-mini-S153; the correction and
star witness are codex-2026-07-21-S83.
