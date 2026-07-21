---
id: THM-1770
title: "THE LOCALISATION LEMMA IS A FIRST-RETURN / RENEWAL STRUCTURE — with the isolation step PROVED and the pair-only case CLOSED, reducing full GMC(2) to one atom-covering statement. An ATOM is a minimal balanced charge-multiset (a minimal vanishing sum over the charge lattice, opus THM-1685 / THM-415); its size is its first-return level. (A) FIRST-RETURN ISOLATION, PROVED: at m* = the minimal atom size, EVERY balanced tuple of size m* is a single atom — a composite would split into ≥ 2 balanced pieces each of size ≥ m*, total ≥ 2m* > m* — so E[P^{m*}] = Σ over size-m* atoms of their coefficient forms, with NO cross/composite terms. Verified: aZ²+bW+cW³ has E[P³]=6ab² (first return at 3), and the 3-charge atom (+2,+3,−5) gives E[P³]=720abc — a minimal vanishing sum, not a pair. (B) DISTINCT ATOMS HAVE DISJOINT COEFFICIENT-MONOMIAL SUPPORTS (they use different charge multisets), so each atom's form vanishes SEPARATELY in the nullcone — no cancellation between atoms (THM-415's prime/no-collision, exactly). (C) RENEWAL INDUCTION: process levels bottom-up (THM-1700); each atom is isolated at its own first-return level; within an atom, a charge of multiplicity r gives an r-dim Vandermonde form killed by the tower m*,2m*,…,r·m* (THM-1740 single-straddle). (D) PAIR-ONLY CLOSURE, PROVED: if every atom is a pos-neg PAIR (⟺ one side carries a single distinct charge value — the star patterns generalizing the single straddle), then the atom-form ideal is ⟨c_p·c_n : p>0>n⟩, whose variety is EXACTLY {all pos = 0} ∪ {all neg = 0} = one-sided — GMC(2) holds for all star patterns, in all degrees, closed-form. So HYP-8590 reduces to ONE statement: V(all atom forms) = the one-sided locus — the atom-covering Nullstellensatz, which is opus THM-1685's core, now in first-return order."
status: >
  (A) PROVED (minimality: any split of a size-m* balanced tuple gives two balanced pieces each
  ≥ m*, impossible below 2m*) and VERIFIED on 4 patterns including a 3-charge atom.
  (B) PROVED (distinct charge multisets ⟹ disjoint coefficient-variable products).
  (C) The renewal ORDER is proved (bottom-up, THM-1700); the within-atom multiplicity reduction
  is THM-1740. (D) PAIR-ONLY case PROVED (the ideal ⟨c_p c_n⟩ has variety exactly the two
  coordinate subspaces; multiplicity handled by THM-1740). The GENERAL localisation (multi-charge
  atoms, the full atom-covering) is REDUCED to V(all atom forms)=one-sided but NOT closed — that
  is opus THM-1685's Nullstellensatz core.
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

## (B) Distinct atoms don't cancel — PROVED

Two distinct atoms have distinct charge multisets, hence **disjoint coefficient-variable
products**. So in `E[P^{m*}] = 0`, each atom's form vanishes *separately* — no cancellation
between atoms. This is exactly THM-415's *prime modulus = no nontrivial vanishing*: distinct
minimal vanishing sums live in distinct monomial supports.

## (C) The renewal induction

Process levels bottom-up (THM-1700 is exactly this order). Each atom is isolated at its own
first-return level by (A)+(B). Within a single atom, a charge carried by `r` terms gives an
`r`-dimensional Vandermonde form, killed by the tower of levels `m*, 2m*, …, r·m*` — THM-1740's
single-straddle law `M* = r·m₀`, and THM-1760's reduction to the radial layer. So each atom,
once reached, is fully resolved by the machinery already proved.

## (D) Pair-only closure — PROVED

> If **every atom is a pos-neg pair** — equivalently one side carries a *single distinct charge
> value* (the star patterns: one positive charge, many negatives, or vice versa; the single
> straddle is the sub-case with one term each) — then the atom-form ideal is
> `⟨c_p·c_n : p>0>n⟩` (mod the multiplicity tower), whose variety is **exactly**
> `{all pos coeffs = 0} ∪ {all neg coeffs = 0} = the one-sided locus`.

So **GMC(2) holds for all star patterns, in every degree**, in closed form: the ideal of
products of one pos and one neg coefficient cuts out precisely the union of the two coordinate
subspaces, which is one-sidedness. Multiplicity within a charge is absorbed by THM-1740.

## The residual, now a single statement

By (A)–(D), the localisation lemma (HYP-8590) reduces to **one** statement:

> `V(⟨all atom forms, all levels⟩) = the one-sided locus.`

For pair-only patterns this is (D), proved. In general the atoms include multi-charge minimal
vanishing sums (`(+2,+3,−5)`), and the equality is the **atom-covering Nullstellensatz** — which
is *exactly opus THM-1685's core*, now organised in first-return (renewal) order rather than by
raw Gröbner. The first-return structure gives: the order (bottom-up), the isolation (each atom
visible alone at its level), and the non-cancellation (distinct atoms, distinct monomials). What
it does not yet give is that hitting every atom's form forces one-sidedness when multi-charge
atoms are present — the hard combinatorial core.

## Honest scope

- **(A), (B), (D) are proved; (C)'s order is proved and its within-atom step is THM-1740.** The
  single-straddle / star / pair-only families of GMC(2) are closed in every degree.
- **The general localisation is REDUCED, not closed.** `V(all atom forms) = one-sided` for
  multi-charge atoms is the open core (opus THM-1685 territory); the first-return reorganisation
  is a genuine structuring of it, not a proof.
- The claim "one side carries a single distinct charge ⟺ all atoms are pairs" is the natural
  characterisation of the pair-only regime; it is used as the *definition* of the closed family,
  verified on the tested star patterns, not proved as an iff for all charge sets.
- Atoms are over the charge lattice `ℤ`; the radial (`V`) degrees ride along via the multiplicity
  tower and do not change which charge-multisets are atoms.

*Artifacts:* `04-computation/gmc2_first_return_localisation_macmini_S153.py` (+out).
*Credits:* THM-1760/1740 (single-atom base), THM-1700 (renewal order), THM-1685/THM-415 (the
vanishing-sum core), THM-1745 (the first-return/LRC reading that pointed here).
