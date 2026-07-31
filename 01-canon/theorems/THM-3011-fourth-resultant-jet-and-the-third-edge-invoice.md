---
id: THM-3011
title: "The fourth resultant log-jet, and discharge of the third-edge invoice"
status: FINITE-EXACT (M=6..62) + VERIFIED-NUMERIC (asymptotics) / NOT PROVED / AWAITING INDEPENDENT HOSTILE AUDIT
source: klein-S428 (lane finding, adversarially verified and independently re-derived in-lane; arithmetic re-checked here)
depends_on:
  - THM-3006-first-gap-wall-is-a-four-band-charge-density-with-all-order-multipole-law
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
  - THM-3003-antipodal-circuit-rigidity-and-the-multipole-spread-criterion
related:
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
  - THM-2982-first-gap-wall-stripped-norm-core-strict-ulc-through-thirty-four
  - THM-2994-first-gap-hurwitz-hermite-biehler-prefix
---

# THM-3011 -- the fourth resultant log-jet and the third-edge invoice

**FINITE-EXACT (`M=6..62`) + VERIFIED-NUMERIC (asymptotics).  NOT PROVED.**

This discharges, at the stated statuses, the obligation left open by THM-3000
section 7 and THM-3003 section 3.  It does **not** produce the symbolic
`P_4(M,U,V)` those files invoice; see section 5.

## 1. The order-four jet engine

Each Macaulay row carries one form, so the `36` selected rows scale by
`n^(deg form)` with total `20(M-1)+10(2M-1)+6(3M-1)=58M-36`, exactly the repo's
chart degree.  Hence `chart(n)=n^(58M-36) det(m_0+m_1 t+m_2 t^2+...)`, `t=1/n`,
and with `X_k=m_0^(-1)m_k` the log-determinant jets are THM-2997 (15) plus the
new fourth line

    L_4 = tr X_4 - tr(X_1X_3) - tr(X_2^2)/2 + tr(X_1^2X_2) - tr(X_1^4)/4.   (1)

The extraneous factor is the explicit `q200^6 * c300 * curvature` of degree
`12M-10=(58M-36)-(46M-26)`, so `L_j(R)=L_j(chart)-6L_j(q200)-L_j(c300)-L_j(curvature)`.

**Control (decisive).**  This engine **re-derives the frozen THM-2997 jet table
from scratch** -- `P_1/D^2`, `P_2/(16D^4)`, `-P_3/(128D^6)` -- exactly at every
width `M=6..24`, fixing `c_1,c_2,c_3=1,16,-128`.  A fully independent rebuild
(full interpolation of the raw `36x36` minor at `58M-36+1` depths, exact flag
division, Newton on the top five coefficients) agrees exactly at `M=6,7,8,12,16,20`.

## 2. The dominant-U slice (VERIFIED-NUMERIC, reconstructed)

    A_4(M) = -(23/10)M^5 + 3M^4 - (47/24)M^3 + (707/3840)M - 15797937/128,
    p_4(R_M) = -4A_4 = (46/5)M^5 - 12M^4 + (47/6)M^3 - (707/960)M + 15797937/32,  (2)

up to a tail of size about `-(3/4)M^8 (3/4)^M`.

**Status discipline.**  (2) was obtained by exact-rational **fit** of an ansatz
through exact large-width values, not by symbolic derivation.  The identical
procedure returns the known THM-2997 (20) polynomials at `j=1,2,3` with error
exactly zero, and the residual `L_4-A_4` falls from `-1.93e3` at `M=100` to
`-4.74e-14` at `M=260` with ratio converging to `-3/4` and no polynomial floor.
That is strong evidence, **not a proof**: a finite reconstruction cannot exclude a
term exponentially small at every width.

**(2) must not be evaluated at small or moderate `M`.**  Below about `M=60` the
exponential tail is comparable to or larger than the polynomial part (at `M=20`
the tail is `~1.0e8` against `p_4(R)=1.3e8`).  Every sign statement below is
computed from exact values, never from (2).

## 3. Equidistribution, and the invoice discharged

(2) confirms the conjecture `sigma_j(R_M)=46M^(j+1)/(j+1)+O(M^j)` at `j=4`, and
sharpens it: the `M^j` coefficient is `-12` for every `j=1,2,3,4`, and the
`M^(j-1)` coefficient is `+47j/24` for `j=2,3,4` -- a dipole of weight `-47/24`
at `s=M` in the multipole reading of THM-3003/THM-3006.

Combining (2) with THM-3006's `w_4/M^5 -> 90211/19440` and jet additivity
`p_4(N)=P_4-w_4`:

    p_4(N_M)/M^5 -> 46/5 - 90211/19440 = 88637/19440 = 4.55951646... > 0,   (3)
    w = m_4/m_1^4 = d^3 p_4/u^4 -> (62/3)^3 (88637/19440)/(131/12)^4
                                = 337994862976/119272468005 = 2.83380455...  (4)

Both re-checked here in exact rationals.  **The invoice is cleared by an
unbounded margin**: THM-3000 (18) needs `w >= -(923/60000)d`, a threshold tending
to `-infinity`, while the true `w` is a positive `O(1)` constant.  Equivalently,
THM-3003's graded hypothesis at `j=4` needs `m_4/m_1^4=o(d)` and gets `O(1)`.

## 4. FINITE-EXACT sign census, and a surprise

Exact rational arithmetic on the THM-2969/2973/2978/2982 atlas cores, cross-checked
two ways (Newton's identities on the exact core coefficients; chart jets minus the
encoded wall):

- `p_4(N)>0` at **every** width `M=6..62`, no exceptions;
- the **third edge** `R_3>R_2` holds exactly for `M=33..62` and fails for `M<=32`
  (`R_3-R_2 = -1.5459e-07` at `32`, `+2.1803e-08` at `33`, `+1.4975e-07` at `34`);
- the second edge `R_2>R_1` first holds at `M=34`
  (`-2.0451e-08` at `33`, `+1.1691055537e-07` at `34`, matching THM-2997 (35));
- `w(M=34)=5.609478` at `d=701`, clearing THM-3000 section 7's box-minimised
  threshold `-7.4613` with margin `13.07`.

> **The third edge turns positive one width EARLIER than the second.**  The
> ordering `1<R_1<R_2<R_3` therefore begins at `M=34` because of the *second*
> edge, not the third; the third edge was never the binding constraint.

## 5. What is NOT done, and one route refuted

- The **symbolic** `P_4(M,U,V)` (expected `~729` monomials, `M`-degree `8`,
  `U,V` weight `<=16`, extending the frozen `27/122/333` counts) was not produced.
  Until it is, (2)--(4) stay VERIFIED-NUMERIC and, for `M>=35`, additionally
  conditional on THM-2997's continuation-wall hypothesis.
- **The sector route is REFUTED.**  `p_4=sum r^4` is termwise nonnegative only
  inside `|arg r|<=pi/8`.  Certified root isolation on the exact atlas cores gives
  `max|arg r| = 26.02 deg` at `M=6`, `56.39` at `M=10`, `66.80` at `M=14`; even
  `|arg r|<=pi/4` fails from `M=9`.  THM-2994's Hurwitz prefix gives only
  `|arg r|<pi/2` and only for `6<=M<=16`, and THM-2999 is a RESERVED stub.
- **No crude route exists.**  `p_4=e_1p_3-e_2p_2+e_3p_1-4e_4` has its two largest
  terms of order `M^8` while `p_4=Theta(M^5)`, a three-order cancellation; strict
  ULC gives only `p_4>=-Theta(M^8)`, missing the required `-0.511M^6` by a factor
  `M^2`.  So THM-3000 section 7 is right that a crude `O(M^5)` bound would suffice,
  but there is no crude way to get one: `p_4` pins `e_4` to relative accuracy
  `M^(-3)`, exactly one more genuine order of the Macaulay expansion.
- Nothing here proves no-return, ULC, GMC(2), or removes the wall invoice.

**Reuse note.**  At `M=12` the THM-2997 encoded wall is one root short of the
proved atlas wall (the `(12,5)` Smith sporadic; exact valuation `27` vs `26`), so
"chart jets minus encoded wall" is off by `5^j` there unless the atlas rule is used.
