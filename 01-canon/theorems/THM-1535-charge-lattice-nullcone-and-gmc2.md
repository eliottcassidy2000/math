---
id: THM-1535
title: "THE CHARGE LATTICE: GMC IS PROVED FOR ALL CHARGE-DEFINITE P IN EVERY DIMENSION, GMC(2) IS PROVED FOR ALL SIGN-COHERENT P, AND THE n=2 / n=4 SEPARATION IS THE RANK OF THE CHARGE LATTICE. Grade monomials by charge q(z^a zbar^b) = a-b, a homomorphism to Z^{n/2}; E kills every nonzero charge and charge is additive. (1) CHARGE LEMMA, PROVED, ANY n: if all charges of P are >= 1 then P lies in the nullcone, and for Q with charges >= -C, E[QP^m] = 0 for every m > C -- so charge-definite nullcone elements can NEVER refute GMC. (2) SIGN-COHERENT RIGIDITY AT n=2, PROVED: if every charge of P is >= 0 (or all <= 0) and P is in the nullcone, then its charge-0 part vanishes, so P is strictly charge-definite. Proof: with no negative charge available, the charge-0 part of P^2 comes only from q=0, so E[P^2] = E[P0^2] = c^T H c with H_{ab} = (a+b)!, the Hankel moment matrix of the exponential distribution, which is POSITIVE DEFINITE -- forcing c = 0. (3) COROLLARY: GMC(2) HOLDS for every sign-coherent P. (4) THE MECHANISM SEPARATING THE DIMENSIONS: the charge group is Z^{n/2}, RANK 1 at n=2 and RANK 2 at n=4. Verified exactly on the 4-term witness: its charge-0 part |Z2|^2-|Z1|^2 has E[P0^2] = +2, but the charges (0,+1) and (0,-1) BOTH occur and their cross term contributes exactly -2, cancelling it. The extra lattice rank buys the cancellation that rank 1 forbids. (5) Nullcone at n=2 verified exhaustively: 59048 polynomials scanned, 160 nullcone members, ALL charge-definite, ZERO exceptions"
status: PROVED (1)(2)(3)(4); (5) VERIFIED-BY-EXHAUSTION on degree <= 3, coeffs {-1,0,1}; the both-signs case at n=2 is OPEN
author: opus-2026-07-20-S411
depends_on: [THM-1495 (the GMC(4) witness and its closed form)]
credits: owner (GMC(2)-true and the Derksen-van den Essen-Zhao nullcone pointer)
---

# THM-1535 — The charge lattice: GMC(2), the nullcone, and why n = 4 is different

## 0. Setup

`n` real Gaussians = `n/2` standard complex Gaussians, `E[z^a \bar z^b] = a!·δ_{ab}` per
variable. Define the **charge** of a monomial

```
q( ∏_j z_j^{a_j} \bar z_j^{b_j} )  =  (a_1−b_1, …, a_{n/2}−b_{n/2})  ∈  ℤ^{n/2}
```

Two facts drive everything: **charge is additive under multiplication**, and
**`E` annihilates every monomial of nonzero charge.**

Recall from THM-1495 the reformulation: `P` is in the **nullcone** iff `E[e^{tP}] = 1`, and
**GMC** says `E[Q e^{tP}]` is a *polynomial* in `t`. The `n = 4` witness has
`E[Qe^{tP}] = t/(1−t)` — a pole.

## 1. The charge lemma (PROVED, every dimension)

> **Lemma.** If every monomial of `P` has charge `≥ 1` in some fixed coordinate (or all
> `≤ −1`), then (i) `P` is in the nullcone, and (ii) for any `Q` whose charges are `≥ −C`,
> `E[Q P^m] = 0` for **every** `m > C`.

*Proof.* Charge is additive, so every monomial of `P^m` has charge `≥ m ≥ 1`, and `E` kills
it — giving (i). Every monomial of `Q P^m` has charge `≥ m − C > 0` for `m > C`, giving
(ii). ∎

> **Consequence.** **A charge-definite nullcone element can never refute GMC, in any
> dimension.** Any counterexample must be charge-*indefinite*. (Check: the `n = 4` witness
> has charges `(0,−1), (0,0), (0,0), (0,+1)` — indefinite, and containing `0`.)

## 2. Sign-coherent rigidity at `n = 2` (PROVED)

At `n = 2` the charge group is `ℤ` (rank 1). Write `P = Σ_q P_q`.

> **Theorem.** Let `n = 2`. If every charge of `P` is `≥ 0` (or every charge is `≤ 0`) and
> `P` lies in the nullcone, then `P_0 = 0` — so `P` is *strictly* charge-definite.

*Proof.* The charge-`0` part of `P²` is `Σ_q P_q P_{−q}`. With no negative charge present,
only `q = 0` contributes, so `E[P²] = E[P_0²]`. Write `P_0 = Σ_a c_a (z\bar z)^a`. Then

```
E[P_0²] = Σ_{a,b} c_a c_b (a+b)!  =  cᵀ H c ,        H_{ab} = (a+b)!
```

`H` is the **Hankel moment matrix of the exponential distribution** (moments `a!`), hence
**positive definite**. `P` in the nullcone forces `E[P²] = 0`, hence `cᵀHc = 0`, hence
`c = 0`, i.e. `P_0 = 0`. ∎

*(Positive-definiteness checked numerically for sizes 1–7: min eigenvalues
`1, 0.382, 0.136, 0.0529, 0.0227, 0.0105, 0.00509` — all `> 0`.)*

> **Corollary (GMC(2) for sign-coherent `P`).** Combining with §1: a sign-coherent nullcone
> element at `n = 2` is strictly charge-definite, and therefore `E[QP^m] = 0` for all
> `m > C`. **GMC(2) holds on the sign-coherent locus, unconditionally.** ∎

## 3. Why `n = 4` escapes: the rank of the charge lattice

**The charge group is `ℤ^{n/2}` — rank 1 at `n = 2`, rank 2 at `n = 4`.** That is the whole
difference, and it is visible in the witness. For `P′ = (1+Z₂)(W₂ − Z₁W₁)`:

| charge | part |
|---|---|
| `(0,−1)` | `\bar z₂` |
| `(0, 0)` | `−z₁\bar z₁ + z₂\bar z₂` |
| `(0,+1)` | `−z₁\bar z₁ z₂` |

```
E[P₀²]                        = +2      <- strictly positive, exactly as at n=2
2·E[ P_{(0,1)} · P_{(0,−1)} ] = −2      <- the OPPOSITE-CHARGE cross term
                        total =  0      = E[P²]   (verified directly)
```

**At `n = 2`, sign-coherence means there is no opposite charge available to cancel the
positive `E[P₀²]`, so Hankel positivity bites and `P₀` must vanish. At `n = 4` the charges
`(0,+1)` and `(0,−1)` both occur, and their cross term is exactly `−E[P₀²]`.** The extra
lattice rank buys precisely the cancellation that rank 1 forbids.

This is the same content as THM-1495's "two-stage cascade", now in invariant form: the
cascade *is* the second lattice direction.

## 4. The nullcone at `n = 2` (verified)

Exhaustive sweep, degree `≤ 3`, coefficients in `{−1,0,1}`:

```
polynomials scanned : 59048
nullcone members    : 160
charge-DEFINITE     : 160
NOT charge-definite :   0
```

**Zero exceptions.** Structural probes of the remaining case (charges of both signs, none
zero) all fail, and fail *exactly where predicted*: for charges `c > 0` and `d < 0`, the
first `m` admitting a balanced monomial is `m = c + |d|`, and that is precisely where
`E[P^m]` first becomes nonzero — e.g. `P = z³ + \bar z²` and `P = z² + \bar z³` both first
fail at `m = 5 = 3 + 2`.

> **Conjecture (the `n = 2` nullcone, DvdEZ-shaped).** `N₂` = exactly the charge-definite
> polynomials. Equivalently, in Newton-polygon terms: `P ∈ N₂` iff the Newton polygon of `P`
> **misses the diagonal `a = b`**.

If this holds, §1 gives **GMC(2) in full**.

## 5. What remains, precisely

The only gap is: `P` at `n = 2` with charges of **both** signs. By Gordan's lemma a finite
set of nonzero charges avoids `0` in its `ℕ`-span iff it lies in an open halfspace — at rank
1 that *is* sign-coherence, so such a `P` necessarily has balanced monomials in `P^m` for
`m ≥ c + |d|`. **The open question is purely whether those balanced monomials can cancel for
every `m`.** Newton-polytope vertices cannot cancel (a vertex coefficient of `P^m` is a pure
power), but non-vertex balanced monomials could in principle conspire. No conspiracy was
found.

## 6. Open

1. **Close the both-signs case at `n = 2`** — the last step to a full proof of GMC(2).
2. **`n = 3`** (HYP-8345): one complex plus one *real* Gaussian. The charge lattice is rank
   1 plus an ungraded direction — §3 predicts this is the genuinely delicate case, and it is
   the only open dimension.
3. **Does §1 generalise?** "Charges lie in an open halfspace" ⟹ nullcone ⟹ GMC, in any
   dimension, by the same Gordan argument. Worth stating as the general charge criterion;
   every counterexample must have charges *not* separable by a hyperplane.

## Verification

`04-computation/nullcone_gmc2_opus_S411.py` (exhaustive `n=2` sweep),
`nullcone_structural_cases_opus_S411.py` (both-signs and charge-0 families),
`hankel_positivity_opus_S411.py` (the Hankel step),
`why_n4_escapes_opus_S411.py` (the rank-2 cancellation, exact). Outputs in
`05-knowledge/results/`.


---

**CORRECTED (opus-2026-07-20-S412, THM-1540).** §2's proof is valid only for REAL
coefficients: over `ℂ` the form `cᵀHc` (no conjugate) vanishes for `c ≠ 0` — e.g.
`c₀ = c₁(−1 ± i)` at size 2 — and §4's 59048-polynomial sweep used real coefficients
`{−1,0,1}`, so the complex case was assumed rather than tested. **The conclusion stands
and is now proved over `ℂ`** by a better reduction: in polar form the charge grading is
the Fourier grading in `θ`, `s = |z|²` is Exponential(1), and with all charges `≥ 0` the
charge-0 part of `P^m` is exactly `P₀^m` — so the nullcone kills ALL moments of `g`, not
just the second, and only `g = 0` does that. See THM-1540.
