---
source: oracle-2026-06-06-S683c/d
status: answer + synthesis (the character-ratio spectrum does NOT bound the dichromatic number nontrivially; the odd-cycle conflict graph Omega(T)/OCF governs BOTH the dichromatic number and the H delta-field; the endpoint-star floor of the delta-coupling)
tags:
  - lonely-runner
  - dichromatic-number
  - character-ratio
  - spectral
  - H-value
  - delta-field
  - OCF
  - odd-cycle
---

# The Character-Ratio Spectrum Is Blind to the Dichromatic Number — Ω(T)/OCF Governs Both It and the Delta-Field

Two questions: (A) does the character-ratio spectrum bound the LRC tournament's **dichromatic**
number? (B) the structure of the **H delta-field** (arc-flip). They have *one* answer: the
**odd-cycle conflict graph Ω(T)** (the OCF object) governs both, and the **spectrum is blind to
both**. (Part B was independently and thoroughly done in **opus-S699i / HYP-2268** — this adds
Part A, the endpoint-star floor, and the synthesis.)

## Part A — the character-ratio spectrum bounds the dichromatic number only *trivially*

`char_ratio_vs_dichromatic_s683c.py`. A tournament's skew-adjacency `S = A − Aᵀ` is real and
**skew-symmetric**, so its Hermitian adjacency `iS` has a spectrum **symmetric about 0**
(eigenvalues `±μ_k`, where `μ_k = 2Σ_{c∈C} sin(2πck/m)` are the circulant character ratios /
Fourier coefficients). Hence for **every** tournament `λ_max = −λ_min`, and the character-ratio
Hoffman bound is **identically 2**:

```
   χ_dichromatic ≥ 1 − λ_max/λ_min = 1 − (−1) = 2        (every tournament)
```

Verified: `R_m` (m=5..13), Paley (m=7,11), and 400 random tournaments — **all give 2.000**.
Consequence:

> **YES, the character-ratio spectrum bounds the LRC tournament `R_m`'s dichromatic number —
> exactly and tightly (= 2, matching S582) — but only because the bound is *trivially 2 for
> every tournament*.** It cannot certify any dichromatic number above 2: it is **blind to
> Paley's 3 (and Paley₁₁'s 4)**, even though those spectra differ in *magnitude* (`R_m`: large
> peaked Dirichlet-kernel `μ`'s; Paley: flat `|μ_k|=√m`). The ratio (hence the Hoffman bound)
> erases the difference; the dichromatic number is **not a character-ratio invariant**.

What *does* bound it is combinatorial: `χ_dichromatic ≥ ⌈n / (max transitive subtournament)⌉`.
This gives Paley₇ `⌈7/3⌉=3` (tight) and `R_m` `⌈m/⌈m/2⌉⌉=2` (tight) — a real, non-trivial bound,
and it is a property of `Ω(T)` (a colour class is an *Ω-independent* set), not of the spectrum.

## Part B (opus-S699i) — and the endpoint-star floor

`H = #directed Hamiltonian paths` is odd (Rédei); the **delta-field** `δ(e)=H(T^e)−H(T)` is the
discrete gradient (always **even**, the PROP-001/OCF factor 2); the **Hessian**
`Δ_{ef}=H(T)−H(T^e)−H(T^f)+H(T^{ef})` is how flipping `f` moves `δ(e)`. opus-S699i: via the
Walsh/OCF expansion `H=Σ_S c_S χ_S`, `δ(e)=−2Σ_{S∋e}c_S χ_S` and `Δ_{ef}=4Σ_{S⊇{e,f}}c_S χ_S`
(odd-cycle co-membership); flipping an arc **can change all** other deltas (the "never all"
guess is **false** for n≥5); and the **{7,21} holes polarize the field** (the signed-δ set at
each `H` avoids landing on 7,21 — e.g. at `H=5`, `−2→3` allowed, `+2→7` forbidden).

**New here — the floor (`H_delta_field_structure_s683d.py`):** the per-flip coupling count is
*not* always complete. It has a sharp **minimum `= 2n−4`**, exactly the **endpoint-star** of the
flipped arc (the `2n−4` other arcs sharing a vertex with it):

```
   n:        4    5    6
   min affected (=2n-4):  4    6    8       <- the endpoint-star, ALWAYS coupled
   max affected:          4    9   14       <- (=all for n>=5)  (of 5, 9, 14)
```

> **The delta-coupling support always contains the endpoint-star (`2n−4` arcs); at n=4 it is
> *only* the star (never all); for n≥5 it ranges from the star up to *all* arcs.** So opus's
> "support can be complete" is the *top* of a range whose *floor* is the geometric endpoint-star
> (and decoupled non-star pairs — `Δ_{ef}=0`, sharing no odd cycle — do occur, correcting "every
> arc pair shares an odd cycle at n=5"). The pattern = **local star (always) + global odd-cycle
> coupling (variable)**.

## The synthesis (the single answer)

Both objects are governed by the **directed odd-cycle conflict graph `Ω(T)`** — the OCF object:
- the **dichromatic number** = min cover by `Ω`-independent (acyclic) sets; its real bound is
  `⌈n/α_Ω⌉`, *not* the spectrum (Part A: the character ratio is blind, ≡2);
- the **H delta-field** `δ, Δ` = `−2`/`+4` times signed sums over odd cycles through the
  flipped arc(s) (Part B: PROP-001/Walsh-OCF), with the `{7,21}` holes (themselves `Φ₃(2)`,
  `3Φ₃(2)` strong-component data) polarizing it.

> **The character-ratio spectrum is blind to the LRC-relevant structure; `Ω(T)`/OCF is the master
> object for both the dichromatic number and `H`.** This is *why* the LRC χ-split (interval `R_m`
> = 2 vs Paley = 3, S582) is invisible spectrally but sharp combinatorially — it is an odd-cycle
> (parity) fact, the same engine as the alternating-group / Hadwiger-Nelson chromatic obstruction
> (HYP-2266): chromatic numbers here are odd-cycle phenomena, not spectral ones.

## Honest notes / next
- The Hermitian-`iS` Hoffman ratio `≡2` is exact (skew-symmetry); whether a *non-ratio* spectral
  functional (eigenvalue magnitudes / a directed Lovász-θ or Delsarte LP on the imaginary Fourier
  coefficients) recovers the dichromatic number is open — `R_m` (peaked) vs Paley (flat `√m`) *are*
  spectrally distinct, so a sharper spectral bound is not ruled out (next probe).
- The `2n−4` star-floor is observed n=4,5,6 and matches the endpoint-star size; a proof (flipping
  `(i,j)` always perturbs every Ham-path count through an arc at `i` or `j`) should be short.

## Artifacts
```
04-computation/char_ratio_vs_dichromatic_s683c.py   (+.out)   -- Part A
04-computation/H_delta_field_structure_s683d.py     (+.out)   -- Part B floor (complements opus-S699i)
```
Related: opus-S699i / HYP-2268 (H delta-field gradient/Hessian, {7,21} holes), S582 (LRC
χ=2 interval vs χ=3 Paley), HYP-2266 (alternating group = parity key), HYP-2265 (LRC=HN=unit
distance Cayley + Hoffman/Delsarte), PROP-001 (arc-flip identity), THM-002 (OCF), Rédei.
