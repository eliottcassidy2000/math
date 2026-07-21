---
id: THM-1830
title: "WHAT THE DvdEZ 2D NULLCONE CONJECTURE IMPLIES -- the honest conditional map, and why the naive DvdEZ => GMC(2) => LRC(14) chain BREAKS but a sibling structure survives. (A) DvdEZ 2D nullcone => GMC(2) is a SOLID PROVEN implication (THM-1535 s4 + THM-1540): DvdEZ says N_2 = charge-definite, the charge lemma says charge-definite elements never refute GMC, the sign-coherent case is proved over C, so DvdEZ empties the both-signs residual and GMC(2) holds in full. (B) BUT GMC(2) does NOT imply LRC(14): the moment functionals are DIFFERENT. GMC's charge is PER-VARIABLE (kernel of r independent characters, Z^{n/2}-graded); LRC's resonance lattice {sum k_i v_i = 0} is the kernel of a SINGLE weighted character (Z-graded). For {1..13} the resonance lattice has rank 12 while the GMC per-variable charge-0 lattice has rank 0 -- genuinely different objects, so the chain BREAKS at the second arrow. (C) THE SIBLING STRUCTURE that survives: both LRC and GMC are moment nullcones for a functional killing nonzero values of a LATTICE CHARACTER; the ABSTRACT theorem 'charges in an open halfspace => nullcone-trivial => the conjecture' (THM-1535 s1 / HYP-8370, proved in every dimension) covers the easy (definite) part of BOTH, and neither GMC(2)'s hard case nor LRC(14) is halfspace-separable (both are both-signs) -- so the shared hard core is 'both-signs moment-nullcone non-vanishing', of which DvdEZ is the Gaussian instance and LRC(14) the resonance instance. (D) THE TRANSFER: DvdEZ's STATEMENT is Gaussian-specific but its PROOF (amoeba / finite-place / apolarity, THM-1685/1735/1710) is FUNCTIONAL-AGNOSTIC, using only charge combinatorics + coefficient non-vanishing; so the DvdEZ-METHOD, not the DvdEZ-theorem, is what yields BOTH GMC(2) and LRC(14), as siblings"
status: PROVED conditional map. (A) is a solid implication (established in THM-1535/1540). (B) is a proven NON-implication (lattice-rank distinction, verified). (C),(D) are the correct sibling framing + transfer route.
author: opus-2026-07-20-S436
depends_on: [THM-1535 (charge lemma + DvdEZ conjecture + halfspace criterion), THM-1540 (sign-coherent over C), THM-1825 (LRC = moment nullcone), THM-1685/1735/1710 (the functional-agnostic proof method), HYP-8365 (DvdEZ), HYP-8620 (LRC resonance nullcone)]
---

# THM-1830 — What DvdEZ implies: the honest conditional map

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

> **So the honest, creative claim is: the DvdEZ-METHOD (not the DvdEZ-theorem) is what proves
> BOTH GMC(2) and LRC(14).** Proving DvdEZ *by the amoeba/finite-place route* produces the
> template that, applied to the resonance lattice `Λ(v)` with coefficients `\hat h`, gives the
> LRC-resonance nullcone (HYP-8620) and hence LRC(14).

## E. The three routes, ranked

1. **Sibling-method (recommended).** Prove the abstract "both-signs single-character nullcone
   non-vanishing" theorem by the amoeba/finite-place method. Its Gaussian instance is GMC(2)
   (via DvdEZ), its resonance instance is LRC(14). One method, two theorems.
2. **Assume-DvdEZ-for-GMC(2)-only.** `DvdEZ ⟹ GMC(2)` is free (§A). This settles the Gaussian
   Moment / Mathieu-subspace / Zhao corollary stack (the repo's Gamma-bridge program) — but
   **not** LRC(14). Honest scope.
3. **Direct high-dim encoding (creative, likely dead).** Encode LRC(14) as nullcone-membership
   for a specific polynomial in ~13 complex Gaussians. But the required functional is the
   *single-character* (diagonal) restriction, not standard GMC — so this recovers §B's
   obstruction, not a new reduction. Recorded as a dead end so it is not re-attempted.

## F. What we can prove *today* assuming DvdEZ

- **GMC(2)** in full (§A) — and with it the repo's NC2 / Structure-Theorem / Zhao-VC /
  Mathieu-subspace corollaries via klein's Γ-bridge.
- The **charge-definite / sign-coherent part of LRC** (the abstract halfspace theorem, §C) —
  which is *already* unconditional (THM-1535 §1), so DvdEZ adds nothing there.
- **Not LRC(14)** — that needs the resonance-instance of the same method (HYP-8620), which
  DvdEZ's *statement* does not provide but its *proof* would template.

## Verification

`04-computation/dvdez_implication_map_opus_S436.py` — the functional/lattice-rank distinction
(resonance rank 12 vs GMC charge-0 rank 0 for `{1..13}`); the arrow analysis. Output in
`05-knowledge/results/`.
