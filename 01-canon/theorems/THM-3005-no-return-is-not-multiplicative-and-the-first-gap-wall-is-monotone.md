---
id: THM-3005
title: "Global no-return is not multiplicative, and the first-gap wall is itself monotone"
status: PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT
source: klein-S428
depends_on:
  - THM-3001-newton-circuit-reversal-involution-and-two-end-curvature-law
  - THM-3003-antipodal-circuit-rigidity-and-the-multipole-spread-criterion
related:
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
  - THM-3004-circuit-sign-change-cluster-law-and-classifier-refutation
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
  - THM-2991-pf-infinity-arbitrarily-delayed-newton-ratio-return
script: 04-computation/gmc_no_return_is_not_multiplicative_and_the_wall_is_monotone_thm3005.py
output: 05-knowledge/results/gmc_no_return_is_not_multiplicative_and_the_wall_is_monotone_thm3005.out
script_sha256: 2906fe03fffa458289bf84c02cf9abec8b697070ff43b330a41944cf4a65ca81
output_sha256: 9605c10fba54f32c9a1bb575f76f7866f866aa48e4a4e5c692958f68796ae6c1
hash_basis: LF-normalized bytes
---

# THM-3005 -- global no-return is not multiplicative, and the first-gap wall is monotone

**PROVED + VERIFIED-EXACT.**

Call `N` **ratio-monotone** (equivalently: `N` has global no-return) when its
circuit `c_k=log(R_k/R_(k-1))`, `k=2..d-1`, is single-signed.  Notation is
THM-3000 section 1.

## 1. The `N times N*` machine

**Theorem.**  For every positive-coefficient `N` of degree `d`, the product
`N N*` has an **antipalindromic** circuit:

    c_k(N N*) = -c_(2d+1-k)(N N*)  for every k.                 (1)

Consequently `N N*` is never ratio-monotone unless its circuit vanishes
identically.  Since `N` ratio-monotone implies `N*` ratio-monotone (THM-3001:
reversal negates the circuit), the **ratio-monotone class is not closed under
multiplication**.

*Proof.*  The root multiset of `N N*` is `{r_i} union {1/r_i}`, so its log-root
multiset is `{ell_i} union {-ell_i}`, which is symmetric about `0` -- hence
symmetric about its own mean.  THM-3003 section 1 (the exact iff) then gives
`R_k(N N*)=R_(2d-k)(N N*)` for all `k`, and the second-difference palindrome is
equivalent to the third-difference antipalindrome (THM-3003 section 1), which is
(1).  An antipalindromic sequence is single-signed only if it is zero.  QED.

So this is not a sporadic counterexample: `N -> (N, N*)` is a **machine** that
converts any ratio-monotone polynomial into a pair of ratio-monotone factors
whose product fails.

**Minimal integer witness** (exhaustive over integer roots `<=6`, each factor of
degree `<=5`):

    A=(n+1)^2(n+2),   B=A*=(n+1)(n+2)^2,   A B=(n+1)^3(n+2)^3,

with circuits `-`, `+`, and `++--`.  `A`'s reciprocal roots `{1,1,1/2}` are
`(1/2){2,2,1}=B`'s roots, so the minimal witness *is* the machine at its
smallest.  Total degree `6`.

The identity (1) is verified independently of monotonicity on four `N`,
including one whose own circuit is `++-` (not monotone): the antipalindrome
holds regardless, only the "both factors monotone" reading needs a monotone `N`.

## 2. The first-gap wall is itself ratio-monotone

Let `W_M` be the explicit wall of THM-2997 eq (9).  Its roots are
`1/2,1,2,...,M` with the five-band multiplicity profile.  Exact computation:

| `M` | `deg W_M` | distinct scales | circuit sign changes | `C(mu_W)` | `C(mu_W*)` |
|---|---|---|---|---|---|
| 6 | 129 | 7 | **0** | `+0.137011` | `-2.974661` |
| 8 | 178 | 9 | **0** | `+0.151536` | `-4.730729` |
| 10 | 229 | 11 | **0** | `+0.165714` | `-6.626403` |
| 12 | 281 | 13 | **0** | `+0.178905` | `-8.720673` |
| 14 | 330 | 15 | **0** | `+0.193836` | `-10.926369` |

The wall's circuit is **strictly positive throughout** -- zero reversals at
degree `330` -- so `W_M` satisfies global no-return outright.  THM-3001's proved
limiting screen holds at every listed width with a margin that **grows** on
both sides; in particular it is far outside MISTAKE-341's `O(1/d)` ambiguity.

**The wall is therefore not the obstruction.**

## 3. Why wall-stripping is necessary, not cosmetic

Combining sections 1 and 2: `R_M=W_M N_M` with `W_M` ratio-monotone, but
monotonicity is not multiplicative.  Hence

- monotonicity of `W_M` and of `N_M` would **not** imply monotonicity of `R_M`;
- monotonicity of `R_M` would **not** imply monotonicity of `N_M`.

So THM-2989/2997's decision to strip the wall and argue about the core is forced
by the algebra, not chosen for convenience.  It also closes off the tempting
shortcut of proving no-return for the resultant directly and inheriting it.

## 4. Dichotomy: separation, not scale count, creates reversals

THM-3004's cluster law needs **well-separated** clusters.  The wall sits at the
opposite extreme -- a near-continuum, roots `1/2,1,2,...,M` with ratios tending
to `1`.  At equal scale count the behaviour is opposite:

    11 distinct scales, near-continuum (wall, M=10):  0 sign changes
    11 distinct scales, ratio 10^6 apart:            19 sign changes (=2m-3)

So the sign-change count interpolates between `0` and `2m-3`, and the controlling
parameter is the **separation**, not the number of scales.  This is the reason
THM-3004's law is not immediately alarming for the first-gap core: the core's
root measure inherits a continuum-like spread from the wall's `-1..-M` support.
That inheritance is a plausible expectation, **not** a proof, and remains an open
obligation.

## 5. The three structural non-properties, and which architectures are dead

Together with THM-3001 and THM-3004, the ratio-monotone class is now known to be:

1. **not characterizable by any reversal-closed hypothesis class** (THM-3001
   section 2) -- positivity, positive real-rootedness, PF-infinity, Hurwitz
   stability and strict ULC are all inert;
2. **not determined by any bounded set of moments** (THM-3004) -- the
   sign-change count is a support-structure property;
3. **not closed under multiplication** (this file).

Dead architectures, therefore: root-class hypotheses alone; low-moment or
low-jet certificates alone; and factor-by-factor induction.  What survives is
what THM-2989/2997 actually do -- a **family-specific, reversal-breaking,
whole-object** invoice -- plus THM-3001's two-scalar screen as a necessary
condition and THM-3003's spread criterion as the cheapest sufficient route for a
fixed edge.

## 6. Boundaries

- Section 1 is exact for every `d`; it needs positivity of the coefficients and
  nothing else.  The "not closed" conclusion needs one ratio-monotone `N` to
  exist, which the witness supplies.
- Section 2 is FINITE-EXACT at `M=6,8,10,12,14` only.  It is **not** proved for
  all `M`, and the wall formula is THM-2997 eq (9) as encoded, which is itself an
  assumption for `M>=35` (the continuation wall invoice is untouched here).
- Section 4's inheritance remark is an expectation, not a theorem.
- Nothing here proves no-return for the core, GMC(2), ULC, or any all-width
  family conclusion.

## 7. Reproduction

    python3 04-computation/gmc_no_return_is_not_multiplicative_and_the_wall_is_monotone_thm3005.py

Three parts, all reporting `True`: the `N N*` machine with the exhaustive minimal
witness, the wall's exact circuit and curvatures through `M=14`, and the
separation dichotomy at equal scale count.
