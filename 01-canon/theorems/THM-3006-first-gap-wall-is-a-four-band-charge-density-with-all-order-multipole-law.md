---
id: THM-3006
title: "The first-gap wall is a four-band charge density, with an all-order multipole law"
status: PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT
source: klein-S428 (with a lane finding independently re-derived and confirmed here)
depends_on:
  - THM-3003-antipodal-circuit-rigidity-and-the-multipole-spread-criterion
related:
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
  - THM-3005-no-return-is-not-multiplicative-and-the-first-gap-wall-is-monotone
script: 04-computation/gmc_wall_multipole_density_all_orders_thm3006.py
output: 05-knowledge/results/gmc_wall_multipole_density_all_orders_thm3006.out
script_sha256: b166efe06ab9d991854df65b93db44ca149fe2b3790004257b3f66eb53b79f25
output_sha256: 43fc6af8220cd4be2102ea1c003453e168e6086ce7b46891828f32c37c8c02e9
hash_basis: LF-normalized bytes
---

# THM-3006 -- the first-gap wall is a four-band charge density

**PROVED + VERIFIED-EXACT.**

THM-2997 (24) records four wall-moment limits.  In the multipole language of
THM-3003 section 3 the wall is a **charge distribution**, and because its
multiplicity profile is piecewise constant its moments are available at **every**
order in closed form, not four at a time.

## 1. The band profile

From THM-2997 eq (9), writing `x=s/M`, the wall's multiplicity density is the
step function

    rho(x) = 28 on (0,1/3],  26 on (1/3,1/2],  24 on (1/2,2/3],  23 on (2/3,1),

together with `O(1)` endpoint atoms -- the `(n+1/2)^6 (n+1)^6 (n+2)^24` head, the
`(n+M)^20` tail, and the mod-`10` root -- which do not affect the leading term.

## 2. All-order multipole law

**Theorem.**  For every `k>=0`,

    lim_(M->infinity) w_k/M^(k+1)
      = [23 + (2/3)^(k+1) + 2(1/2)^(k+1) + 2(1/3)^(k+1)]/(k+1).   (1)

*Proof.*  `w_k/M^(k+1) -> int_0^1 rho(x)x^k dx`, and

    28/(k+1)(1/3)^(k+1)
  + 26/(k+1)[(1/2)^(k+1)-(1/3)^(k+1)]
  + 24/(k+1)[(2/3)^(k+1)-(1/2)^(k+1)]
  + 23/(k+1)[1-(2/3)^(k+1)]

collapses to (1).  The coefficients `1,2,2` are exactly the **downward jumps** of
the multiplicity profile at `x=2/3,1/2,1/3` (`24-23`, `26-24`, `28-26`).  QED.

## 3. Controls and new values

`k=0,1,2,3` reproduce THM-2997 (24) exactly: `76/3`, `145/12`, `2551/324`,
`1681/288`.  Direct summation over eq (9) at `M=600,1200,2400` converges
monotonically to (1) for every `k<=6`.

**New, extending (24) by three indices:**

    w_4/M^5 -> 90211/19440    = 4.640483539...
    w_5/M^6 -> 179795/46656   = 3.853630830...
    w_6/M^7 -> 3229771/979776 = 3.296438165...

## 4. Scope, and the obligation this does NOT discharge

- (1) is about the **wall** only, which is explicit; it is conditional on
  THM-2997 eq (9) as encoded (the continuation wall invoice is untouched).
- It supplies **one half** of the third-edge invoice of THM-3000 section 7 /
  THM-3003 section 3: with jets additive over `R_M=W_M N_M`, the core's fourth
  power sum is `p_4(N)=P_4-w_4`, and `w_4` is now known exactly.  `P_4` -- the
  resultant's own fourth log-jet -- is **not** supplied and remains the open
  obligation.
- **The root-modulus route is now known to be blocked from the published data.**
  A Hurwitz positive-coefficient polynomial exists with the same degree and the
  same first three power sums as the core, whose largest root modulus is
  `Theta(M^2)`; every classical coefficient bound (Cauchy, Kojima, Fujiwara) is
  at least `a_(d-1)/a_d ~ (131/12)M^2`.  So THM-3003's spread criterion
  `kappa=o(d^(1/4))` cannot be discharged from `d, sigma_1, sigma_2, sigma_3`
  alone -- a genuinely new jet (or a root-location argument outside the
  coefficient package) is required.  Empirically the true `max|r|` for the core
  is `Theta(M)` (about `34.5+2.53M` over `M=6..31`), comfortably inside
  `o(M^(5/4))`; that is FINITE-EXACT evidence, not a proof.
- Nothing here proves no-return, ULC, or GMC(2).
