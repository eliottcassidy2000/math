---
id: THM-650
title: TWO CLOSED FORMS FOR THE CROSSING TRILOGY — (A) the pair-space rank formula PROVED: the mod-2 quadratic part of the cylinder crossing form equals the share-an-endpoint adjacency of E(K_{m,n′}) at every generic twist, and its GF(2) rank is m + n′ − 1 − [m ≡ n′ (mod 2)]; (B) the n=8 σ-fixed crossing minimum equals 20 = 6 + 2·7, PROVED by a human-checkable certificate: the 6 mirror-self crossings are forced, and 7 explicit edge-disjoint odd cycles pin the quotient MaxCut at 25
status: PROVED (both; complete proofs below). (A): three-case coefficient computation + kernel identification; verified in EXACT Fraction arithmetic at 3 twists × 9 shapes (an earlier float artifact at boundary-integral endpoints is part of the record). (B): forced-crossing lemma + exhaustive maxcut = 25 + the explicit 7-cycle packing certificate (5 triangles + 2 pentagons, listed).
source: mac-mini-2026-07-07-S53 (HYP-5127; owner: prove the rank formula and the n=8 quotient-min = 20 closed form)
depends_on:
  - THM-649   # the crossing trilogy (both statements were its shaped follow-ups)
related:
  - HYP-5052  # opus-S142 (book-side parity law)
  - HYP-4961  # klein-S172 (constructive self-loop existence — the parallel even-n program)
external: MaxCut/odd-cycle-packing duality (the certificate gadget); Zarankiewicz/Kleitman parity context.
---

# THM-650 — Two closed forms

## (A) The pair-space rank formula

**Theorem.** For the winding cube of the two-circle model on K_{m,n′} at generic twist,
the mod-2 crossing form is exactly degree ≤ 2, its quadratic part B satisfies
B[e,f] = 1 ⟺ e and f share an endpoint, and
> **rank_{GF(2)}(B) = m + n′ − 1 − [m ≡ n′ (mod 2)].**

**Proof.** *Coefficients.* The pair (e,f)-coefficient is the mod-2 second difference of
the pair's crossing count c(w_e, w_f) = #(ℤ ∩ open(x, x+δ)), which depends on the
windings only through δ. Since c(1,1) = c(0,0) (same δ), the coefficient is
c(1,0) + c(0,1) mod 2, i.e. c(δ₀−1) + c(δ₀+1).
- *Disjoint pairs:* all endpoints non-integral (x = (i₁−i₂)/m, x+δ₀ = (j₁−j₂)/n′ + ℤ,
  both with distinct indices), so c(δ) = ⌊x+δ⌋ − ⌊x⌋ up to sign and
  c(δ₀+1) + c(δ₀−1) = 2(⌊x+δ₀⌋ − ⌊x⌋) ≡ 0. Coefficient 0.
- *Shared inner (x = 0):* c(δ₀) counts ℤ ∩ (0, δ₀)-oriented = 0 for |δ₀| < 1;
  c(δ₀−1) = [δ₀ < 0], c(δ₀+1) = [δ₀ > 0]; sum = 1. Coefficient 1.
- *Shared outer (x + δ ∈ ℤ structurally):* base interval (x, 0)-oriented with
  |x| < 1: c = 0; the two single shifts give intervals (x, ±1) containing 0 exactly
  when its side matches the sign of x: [x<0] + [x>0] = 1. Coefficient 1. ∎(support)
*Rank.* Hence B = (J_m−I)⊗I + I⊗(J_{n′}−I) ≡ J_m⊗I_{n′} + I_m⊗J_{n′} (mod 2), i.e.
(BX)_{ij} = rowsum_i(X) + colsum_j(X) on X ∈ GF(2)^{m×n′}. Kernel: all row sums and all
column sums equal a common constant c; the c = 0 space (even row/col sums) has dimension
(m−1)(n′−1), and the c = 1 coset exists iff m ≡ n′ (mod 2) (total-sum consistency),
adding 1. So rank = mn′ − (m−1)(n′−1) − [m ≡ n′] = m + n′ − 1 − [m ≡ n′ (mod 2)]. ∎
*Verification.* Exact Fraction arithmetic, twists {137/1000, 61/100, 923/1000},
(m,n′) up to (4,5): the support statement holds at every twist, ranks
2,4,4,4,6,6,6,6,8 all match. (Record: a first float run produced spurious non-star
coefficients at some twists — the shared-outer endpoints are structurally INTEGRAL and
float noise corrupts the boundary test; exact arithmetic is mandatory here.)

**Corollary.** The pair space never vanishes for m + n′ ≥ 4 (rank ≥ 2), so cylindrical
crossing parity is never invariant over the winding cube — unlike the book model's
odd-n affine law. The rank is also the number of independent "parity channels" any
cylindrical Kleitman-type argument must control.

## (B) The n=8 σ-fixed minimum = 20, closed form

**Theorem.** min{Q(t) : t ∈ Fix(σ)} = 20 = 6 + 2·7 (> Z(8) = 18): the mirror break of
THM-649(II) with an explicit certificate.

**Proof.** *Forced part.* The 70 crossing pairs of the n=8 tile chords split under σ
into 32 two-element orbits and 6 fixed pairs; every fixed pair is a MIRROR-SELF pair
{c, σc} with c crossing σc (verified: exactly 6 such; and the three σ-fixed chords
(1,8),(2,7),(3,6) are nested, contributing 0). On Fix(σ), c and σc carry the same bit,
hence the same page — each mirror-self pair contributes a crossing ALWAYS: f = 6 forced.
*Quotient part.* Q|Fix = 6 + 2q where q = same-page count over the 32 orbit
representatives, a function of the 12 quotient bits: q = 32 − cut(w) on the quotient
crossing multigraph G_q (12 vertices, 32 edges). Exhaustively, maxcut(G_q) = 25.
*Human-checkable lower bound.* Seven EDGE-DISJOINT ODD CYCLES of G_q (each forces ≥ 1
uncut edge in any cut):
  Δ₁ (4,8)(4,10)(8,10); Δ₂ (5,8)(5,9)(8,9); Δ₃ (2,6)(6,9)(2,9); Δ₄ (1,4)(4,5)(1,5);
  Δ₅ (3,4)(3,6)(4,6); C₆ (6,8)(6,11)(7,8)(7,10)(10,11); C₇ (3,5)(3,7)(5,11)(7,9)(9,11).
Hence maxcut ≤ 32 − 7 = 25, so q ≥ 7; the σ-fixed optima attain q = 7.
Therefore min Q|Fix = 6 + 14 = 20. ∎

**Remark (the pairing principle completed).** THM-649(II)'s mechanism — σ-paired
crossings contribute doubly, quantizing the symmetric landscape — is now a closed-form
proof: the forced-6 comes from symmetry alone (pairing-with-sign-flip: same bit ⟹ same
page ⟹ crossing), and the 2·7 from odd-cycle combinatorics of the quotient. No census
remains in the statement; the exhaustive runs are confirmations.

## Files
`rank_formula_and_quotient20_macmini_S53.py` + `.out` (exact verification, the packing
search, and the certificate list).
