---
id: THM-2550
title: "Double nondegeneracy of the canonical typed row: positive owner-loop drift and a non-replica lawful table"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT (two independent exact companions,
  each with internal cross-checks anchored byte-exactly to the audited
  THM-2541 artifact values; python/-O identical); INDEPENDENT HOSTILE AUDIT
  ACCEPTED BOTH COMPUTATIONS WITH THE MISTAKE-283 SCOPE REPAIR.  Resolves
  both computational arms of HYP-9032's decision tree in the nondegenerate
  direction on the canonical typed row.  The two arms share row and word
  labels, not a clock, packet, marked triangle, or supplied common-root
  service.  The row is TYPED (THM-2309 (25)), not an asserted scalar cover;
  no scalar row is removed; the ledger stays 165; LRC(14) remains OPEN.
source: opus-2026-07-27 (the two decisive single-row computations
  named by THM-2368 Sec 8 / THM-2365 (32) and by THM-2512 (26)-(28))
depends_on:
  - THM-2365 (drift dichotomy; eqs (19b), (31), (32))
  - THM-2449 (lawful owner response table; clock law (25)-(28))
  - THM-2512 (replica dichotomy; ignition (17)-(20), (27)-(28))
  - THM-2461 (source-refined packet, eq (9))
related:
  - THM-2541 (full target-plane support, same row/word/triangle)
  - THM-2466 / THM-2471 (nearby chain; their missing common-root hypotheses are not supplied here)
  - THM-2367 (the zero-drift hostile shape that did NOT occur)
  - HYP-9032 (the transplant trichotomy; both arms resolved here)
scripts:
  - 04-computation/lrc14_owner_loop_drift_typed_row_opus_20260727.py
  - 04-computation/lrc14_replica_dichotomy_typed_row_opus_20260727.py
outputs:
  - 05-knowledge/results/lrc14_owner_loop_drift_typed_row_opus_20260727.out
  - 05-knowledge/results/lrc14_replica_dichotomy_typed_row_opus_20260727.out
---

# THM-2550 -- the canonical typed row is doubly nondegenerate

> **SCOPE CORRECTION (MISTAKE-283).**  The two nondegeneracy computations are
> exact, but Part (A)'s `k=2` drift and Part (B)'s large-clock interaction are
> not one packet or clock.  The former does not directly feed THM-2466 or
> THM-2471: their common oriented root/service and collision-root
> identification hypotheses remain open.

Row: `THM-2309 (25)`, `w = (1,14,27,40,53,66,13,2197,742586)`, owner
`j = 1`, word stratum `sigma = {a}` (the THM-2541 setting).

## (A) Owner-loop drift is POSITIVE (both tensors)

The THM-2365 (19b) successor masses of the source-refined owner packet
`f = 1_Q P^2 1_{E_1}` are decisively nonconstant: `168/168` nonzero
lawful twists differ from `S(0,0)`; witness
`S(0,0) = 254882231/116398499035` vs
`S(0,1) = 547500161177/241992479493765` (exact). Drift energies:

```text
D_{H_F} = 99817232068971546871095390313 / 1391555314836478909802943689850619200  (packet)
D_{H_E} = 9497453727128823229273 / 7328279741345331835978099188                  (bare owner)
```

both `> 0` exactly; `1192/2197` resp. `2016/2197` circulant-law
violations; the diagonal law `H(t,s,t) = 0` and the (19b) row identity
hold at all `169` twists (structural sign check). The THM-2368 Sec 8
escalation ladder terminates at its first rung: **no zero-drift
circulant hostile arises on this row.** Since the BARE tensor is
drift-positive, THM-2365 (31) branch 1 applies: at every sufficiently
large delayed clock some `gcd(m,91)=1` and ordinary `X` have a nonzero
target fibre on this row.  This conclusion does not supply THM-2466's
common oriented root base and positive service density, nor identify
THM-2471's first-collision root with the endpoint/deep root used here.

## (B) The lawful table is NON-REPLICA (host array ignited)

The THM-2449 response table `A^R` was computed exactly at four clocks
(`k = 2` control, and `6, 16386, 32766` in the guaranteed class
`k = 6 mod 16380`, via an exact prefix-count identity for
`R = 13^16386`); `(M_rho, C_rho)` extracted from two clocks and
verified at a third plus the exact product-measure formula. BOTH
additive defects are nonzero -- `d_{M_rho} != 0` and `d_{C_rho} != 0`
-- with ALL `91/91` interaction entries nonzero in each matrix
(witnesses at `(0,0)`: `~5.245e-04` exact rational, resp. `~62.08`
exact rational). By THM-2512 (27)-(28), all `5,184` primitive cut
coefficients fire at every sufficiently large lawful clock in the
class (at most one exception): **the three transplant inputs (T1)
row-zero, (T2) dichotomy, (T3) magnitude floor now hold on a
live-branch array**, with the previously missing denominator floor
supplied by `lcm(den) = 6823601362619167115089900800`.

## Sharp side-findings

- The THM-2449 (18) delta anchor FAILS on this row (`A_{ell,0} > 0`
  for all `ell`) -- consistent with the no-scalar-cover guardrail;
  the dichotomy is unaffected (only the replica-form normalization
  (10) would need the anchor).
- The clock law `A = M + C/R` fails exactly at the sub-threshold
  clock `k = 2`: the `k >= K` scope of THM-2449 (25)-(28) is sharp.

## What follows and what does not

The canonical typed row now carries all three nondegeneracy
certificates: full target-plane support (THM-2541), positive
owner-loop drift (A), and an ignited live host array (B).  HYP-9032's
prediction P1 is CONFIRMED in the componentwise sense that neither arm
is degenerate.  They share the typed row and word label, but Part (A)
uses the `k=2` packet while Part (B) uses `k=6 mod 16380`; no common
clock, packet, marked triangle, or ancestry fibre has been proved.
NOT implied: any universal bypass, direct entry into THM-2466/2471,
landing on the prior marked triangle, an all-91-unit mask, exclusion of
any scalar row, or LRC(14).  The next obligations are separately to
transport the drift-positive instance through THM-2368 (37)'s
phase/event word with the required common-root service, and to realize
the ignited interaction on one nonnegative Boolean ancestry carrier.
