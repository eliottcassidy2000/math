---
id: THM-1830
title: "Historical DvdEZ implication map: lattice distinction survives; claimed LRC transfer does not"
status: >
  SUPERSEDED / PARTLY REFUTED. The conditional NC2->GMC(2) charge argument
  and the distinction between Gaussian charge and LRC resonance lattices
  survive. THM-2022 now proves NC2/GMC(2) independently. The claim that an
  amoeba/finite-place DvdEZ method is functional-agnostic and therefore proves
  LRC(14) was never established: no predicate-preserving transfer or LRC
  noncancellation theorem is supplied. Sections C--F are historical proposals,
  not proved consequences.
author: opus-2026-07-20-S436
depends_on: [THM-1535 (charge lemma + DvdEZ conjecture + halfspace criterion), THM-1540 (sign-coherent over C), THM-1825 (LRC = moment nullcone), THM-1685/1735/1710 (the functional-agnostic proof method), HYP-8365 (DvdEZ), HYP-8620 (LRC resonance nullcone)]
---

# THM-1830 — What DvdEZ implies: the honest conditional map

> **CURRENT CORRECTION.** Keep only the elementary conditional charge
> implication and the lattice-rank obstruction below. The claimed
> “functional-agnostic” transfer to LRC is conjectural and does not prove
> LRC(14). NC2/GMC(2) are now closed by THM-2022; LRC(14) and JC(2) remain
> separate. This file is provenance, not current canon.

Owner asked: assuming the Derksen–van den Essen–Zhao 2D nullcone conjecture, what can we
prove — can we chain `DvdEZ ⟹ GMC(2) ⟹ LRC(14)`? **The first arrow holds; the second does
not; but a sibling structure and a transferable proof survive.** The honest map:

## A. `DvdEZ ⟹ GMC(2)` — SOLID (proven implication)

The **DvdEZ 2D nullcone conjecture** (HYP-8365, THM-1535 §4): the `n=2` nullcone `N₂` is
*exactly* the charge-definite polynomials (Newton polygon misses the diagonal). Then:

- **Charge lemma** (THM-1535 §1, PROVED): a charge-definite nullcone element can never refute
  GMC, in any dimension.
- **Sign-coherent case** (THM-1540 §3, PROVED over `ℂ`): the non-both-signs part is settled.
- **DvdEZ** empties the *both-signs* residual (says it contains only charge-definite = nothing
  new).

Therefore **`DvdEZ ⟹ GMC(2)` in full.** This arrow is real and already in canon.

## B. `GMC(2) ⟹ LRC(14)` — FALSE (proven non-implication)

The tempting chain breaks at the second arrow, because **the two moment functionals are
different lattices**:

| | moment functional | grading | for `{1,…,13}` |
|---|---|---|---|
| **GMC** | `E_Gauss[z^a\bar z^b] = a!\,δ_{ab}` (per variable) | **per-variable** charge `(a_i−b_i) ∈ ℤ^{n/2}` — kernel of `r` **independent** characters | rank-`0` charge-0 lattice |
| **LRC** | `E_circle[e^{2πikt}] = δ_{k,0}` | **single weighted** character `(k_i) ↦ Σk_iv_i ∈ ℤ` | rank-**12** resonance lattice `Λ(v)` |

The resonance lattice `Λ(v) = {(k_i) : Σk_iv_i = 0}` is the kernel of **one** character; GMC's
charge-0 is the kernel of **many**. They are genuinely different objects (rank 12 vs 0 for the
13 speeds). **So GMC(2) does not logically imply LRC(14)** — no functorial reduction turns the
per-variable Gaussian nullcone into the single-character resonance nullcone. The naive chain is
**not valid**, and asserting it would be wrong.

## C. The sibling structure that *does* survive

Both are **moment nullcones for a functional killing the nonzero values of a lattice
character.** The abstract theorem I proved (THM-1535 §1 / HYP-8370), **"charges in an open
halfspace ⟹ nullcone-trivial ⟹ the conjecture holds,"** holds **in every dimension for any
such functional** — it covers the *easy* (definite) part of **both** GMC and LRC. And:

- GMC(2)'s hard case is **both-signs** (charges `±1`).
- LRC(14)'s hard case is **both-signs** (`\hat h(k) = \sin(2πkλ)/(πk)` alternates — THM-1825).

> **The shared hard core is: a both-signs moment-nullcone sum cannot vanish identically unless
> the generator is a single character.** DvdEZ is the *Gaussian instance*; LRC(14) is the
> *resonance instance*. They are **siblings**, not a chain.

## D. The transfer: it is the METHOD, not the theorem

DvdEZ's *statement* is Gaussian-specific. But the **proof method** I have been building —
**amoeba/multinomial-radius separation (THM-1710), the finite-place mod-`p` certificate
(THM-1735), the `k`-nomial Nullstellensatz emptiness (THM-1685), positivity (THM-1705)** — is
**functional-agnostic**: it uses only (i) the charge-representation combinatorics of the lattice
and (ii) non-vanishing of the coefficient products. Both are present for the LRC resonance
functional (with `\hat h(k)` as coefficients).

> **Historical proposal, not proved:** an amoeba/finite-place proof of DvdEZ
> might suggest tests for the resonance lattice. No map preserving the LRC
> covering predicate or coefficient-weighted noncancellation was constructed,
> so it cannot be transferred as a theorem and does not imply LRC(14).

## E. The three routes, ranked

1. **Sibling-method (historical proposal).** Prove the abstract "both-signs single-character nullcone
   non-vanishing" theorem by the amoeba/finite-place method. Its Gaussian instance is GMC(2)
   (via DvdEZ), its resonance instance is LRC(14). One method, two theorems.
2. **Assume-DvdEZ-for-GMC(2)-only.** `DvdEZ ⟹ GMC(2)` is free (§A). This settles the Gaussian
   Moment / Mathieu-subspace / Zhao corollary stack (the repo's Gamma-bridge program) — but
   **not** LRC(14). Honest scope.
3. **Direct high-dim encoding (creative, likely dead).** Encode LRC(14) as nullcone-membership
   for a specific polynomial in ~13 complex Gaussians. But the required functional is the
   *single-character* (diagonal) restriction, not standard GMC — so this recovers §B's
   obstruction, not a new reduction. Recorded as a dead end so it is not re-attempted.

## F. What the conditional map recorded historically

- **GMC(2)** follows logically from the stated NC2 classification (§A), but
  its current unconditional proof is THM-2022, not the refuted Gamma bridge.
- The **charge-definite / sign-coherent part of LRC** (the abstract halfspace theorem, §C) —
  which is *already* unconditional (THM-1535 §1), so DvdEZ adds nothing there.
- **Not LRC(14)** — that needs the resonance-instance of the same method (HYP-8620), which
  DvdEZ's *statement* does not provide but its *proof* would template.

## Verification

`04-computation/dvdez_implication_map_opus_S436.py` — the functional/lattice-rank distinction
(resonance rank 12 vs GMC charge-0 rank 0 for `{1..13}`); the arrow analysis. Output in
`05-knowledge/results/`.
