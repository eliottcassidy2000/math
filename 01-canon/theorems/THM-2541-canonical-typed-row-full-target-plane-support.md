---
id: THM-2541
title: "Canonical typed-row 169-twist target variance positivity and full target-plane support"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY REFEREED (three ways: klein's
  exact inverse DFT [MSG-2172-to-opus], mac-mini's independent rerun
  from origin/main [MSG-2174-to-opus], and the artifact's own section
  [9b] reproducing klein's 169/169 census with three exact controls).
  SCOPE GUARDRAILS: the row is the canonical TYPED row THM-2309 (25) --
  typed, NOT asserted to be a scalar cover (no covering row can be
  exhibited without refuting LRC(14)); the aggregates are the
  unrestricted A(q), NOT the all-91-unit projector B(q) of THM-2334
  (49); the variance is not a rational number (it is
  [delta_hat(m)/(4 pi^2 XY)]^2 times a nonnegative real cyclotomic);
  no scalar row is removed and LRC(14) remains OPEN.
source: >
  opus-2026-07-26/27 (computation and certificate); strengthened to full
  target-plane support at klein's request (MSG-2171/2172); independently
  rerun by mac-mini (MSG-2174).
depends_on:
  - THM-2334 (twist bank, target aggregates, eqs (28)-(44), (42)-(43))
  - THM-2309 (the canonical typed row, eq (25))
  - THM-2349 (delayed-word universe)
  - THM-2305 (canonical word order)
related:
  - THM-2479-two-colour-trichotomy-marked-91-unit-edge-incidence
  - HYP-9031-d5-h1-dictionary-lrc-word-current-vs-jc-flux (prediction 1)
  - THM-2519 (klein's root-side conditional variance; comparison OPEN)
script: 04-computation/lrc14_169_twist_variance_opus_20260726.py
output: 05-knowledge/results/lrc14_169_twist_variance_opus_20260726.out
script_sha256: c0aa9c55c3e7d38dc28b4705e58c776a979f17d2874d1b762f96b97d2364e5e9
output_sha256: 8c2cc1df0ff9e2ccc8b985daae1ac7041256911ebf49179b66c5210446a08bff
hash_basis: working-tree bytes (LF)
---

# THM-2541 -- the first nonzero-target survival object, at full support

## Statement

On the canonical typed row `THM-2309 (25)`,
`w = (1, 14, 27, 40, 53, 66, 13, 2197, 742586)` (guard `H=1` odd,
valuation profile `(1,3,5)`), with owner `j=1`, the first
positive-measure delayed word `sigma = {a}` in THM-2305's fixed order at
clock `k=2` (`R = 169`), and the first fully certified marked triangle
`X = 13`, `m = 1`, `Y = 742599`:

1. **(Variance positivity)** The THM-2334 (42) nonzero-target mass is
   strictly positive:
   `sum_{q != 0} |A(q)|^2 = (1/169) sum_ell |H(ell) - Hbar|^2 > 0`.
   All `168` nonzero twists satisfy `H(ell) != H(0)`.
2. **(Full target-plane support)** Strictly stronger: EVERY one of the
   `169` target aggregates `A(q)`, `q in F_13^2`, is a nonzero complex
   number; in particular all `168` nonzero targets survive.

## Certificate (exact; floats descriptive only)

`gamma(ell) = e_13(m ell_8) (2 pi i X)hat(E_Q^ell)(X) conj((2 pi i Y)hat(E^ell)(Y))`
is a cyclotomic integer in `Z[zeta_NN]`,
`NN = 2^4 3^3 5 7^2 11 13^8 53 = 50334435734703120`. Ring homomorphisms
`zeta_NN -> h` into `F_p` at `p_1 = 7 NN + 1` and `p_2 = 19 NN + 1`
(orders verified) decide nonvanishing: a nonzero image proves the
cyclotomic number nonzero. Positivity: `gamma(ell) - gamma(0)` nonzero
mod `p_1` for all `168` nonzero twists (first witness distinguished mod
`p_2` as well). Full support: the exact inverse DFT
`Anum(q) = sum_ell gamma(ell) z13^{-(ell_TA q_a + ell_TB q_b)}`,
`z13 = h^{NN/13}`, has all `169` values nonzero mod `p_1`, with controls
(A) `sum_q Anum = 169 gamma(0)`, (B) sign-reversal census invariance,
(C) forward reconstruction at two twists -- all exact PASS. Gauge
(representative independence) `H(ell + w) = H(ell)` exact PASS at four
cosets. The `ell -> (ell_TA, ell_TB)` pairing is verified to biject the
coset space onto `F_13^2`.

## Why this matters (and does not overclaim)

This is the first evaluated instance of THM-2334's designated cheapest
decisive object, and it comes out ALIVE at full support: the marked
current of a canonical typed row sees every nonzero target. By the
`Phi_91`-Galois lever (Gal(Q(zeta_91)/Q) transitive on primitive
characters), rational NONCONSTANCY statements of this shape are exactly
what the `7 (x) 13` mixed-character crux consumes; the remaining
distance to the ledger is the all-91-unit projector `B(q)` (THM-2334
(49)), word-support masking, bounded visible height, terminal-component
phase transport, and scalar-row exclusion -- none of which this theorem
touches. klein's THM-2519 root-side conditional variance measures a
DIFFERENT conditioning (comparison typed OPEN; see MSG-2173).
