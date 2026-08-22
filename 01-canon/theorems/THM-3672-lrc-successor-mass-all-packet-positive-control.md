---
id: THM-3672
title: "LRC successor-mass all-packet positive control"
status: >
  FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the
  typed non-cover row of THM-3669, all thirty legal ordered graft charts, and
  hence all 120 owner-pivot packets, have a strictly negative THM-2365
  three-site successor-mass defect.  The calculation uses the actual
  successor weight 2-d(13cx), not plain marked mass.  This validates the
  THM-3670 rational gate on one hostile control.  THM-3674 converts the
  strongest exact defect into drift at least 1.078089085e-12 and
  nonzero-deep, nonzero-target energy at least 8.292992962e-14.  This proves
  no covering-row statement and does not prove LRC(14).
source: kps-s193 / Parfit successor-ledger extraction, 2026-08-21
audit: >
  PASS -- Socrates independently reproduced normal and optimized transcripts,
  all hashes, the parent pin, half-open interval intersections, two-window
  danger-prefix integration, THM-2365 sign convention, thirteen successor
  identities, thirty strict negative defects, the canonical reduction, and
  the 6*5*4 packet count.  No correction was required and the non-cover scope
  was confirmed.  The later energy invoice is exact rational substitution
  into independently hostile-audited THM-3674.
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-3666-lrc-owner-pivot-dual-pair-swap-twist-basis
  - THM-3669-lrc-typed-control-all-packet-three-twist-defects
  - THM-3670-lrc-successor-mass-transfer-and-thirteen-number-gate
  - THM-3674-sharp-successor-variance-drift-and-target-energy-tariff
script: 04-computation/lrc_successor_mass_all_pair_swap_control_thm3672.py
output: 05-knowledge/results/lrc_successor_mass_all_pair_swap_control_thm3672.out
script_sha256: a191f934b20494b98e878fc5504a328c47d4cb6a92ea332246f782fac80c01c8
output_sha256: 0d72ff662c908d9d9934976edce67895135c17a1eb0ba409e0abc382563b52c6
successor_ledger_sha256: a41ec0c3502221eb1e0e17bd1d941817fdba134f93c796050a4b454772fee668
defect_ledger_sha256: 9129d6dc7c7c03db0b99cd55f8ee504f65448b87ebdac008e0cc66356787014b
hash_basis: raw LF bytes
---

# THM-3672 -- the exact successor detector fires on every control packet

**FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3669 proves that every owner-pivot packet on one typed non-cover row has a
nonzero complex target-current defect.  THM-3670 instead offers a cheaper
rational sufficient test.  This companion evaluates that test exactly.

## 1. Control row and orientation

Use the fixed data

```text
w=(1,14,27,40,53,66,13,2197,742586),
owner=c1, target blockers a=c2 and c=c3,
word={a}, R=169,
N=50334435734703120.                                (1)
```

For graft labels `k,l` in `U={0,...,5}`, put

```text
alpha_k=e_a-e_k,
beta_l =e_c-e_l.                                   (2)
```

THM-2365 translates present factors by `-s alpha_k/13-t beta_l/13`.
Consequently the companion's negative target dipoles are exactly the positive
THM-2365 coordinates

```text
A_k=S_(k,l)(1,0),
B_l=S_(k,l)(0,1),                                  (3)
```

not the deep-shift colour `r` in THM-2365's tensor.  Here

```text
S_(k,l)(s,t)
 =integral_T F_(s,t)(x)(2-d(13cx)) dx.             (4)
```

The successor speed, one danger-window length, and danger period on the
common numerator grid are

```text
13c=9653618,
window=372432060,
period=5214048840.                                  (5)
```

## 2. Thirteen exact masses

For each row below the tuple is

```text
(successor numerator, marked-mass numerator,
 danger-intersection numerator, component count).  (6)
```

The exact ledger is

```text
S0  (110219232915792,60084076348296,9948919780800,188056)
A0  ( 99299969997228,54135630512964,8971291028700,169431)
A1  ( 99299969997228,54135630512964,8971291028700,169431)
A2  ( 99299969997228,54135630512964,8971291028700,169431)
A3  ( 99050735597814,54000012243177,8949288888540,169006)
A4  ( 96962693739912,52861450868226,8760207996540,165443)
A5  ( 96754696164126,52747826144013,8740956123900,165088)
B0  (113880033524816,61887542465528,9895051406240,186674)
B1  (113880033524816,61887542465528,9895051406240,186674)
B2  (113880033524816,61887542465528,9895051406240,186674)
B3  (113525486798282,61695142630901,9864798463520,186093)
B4  (111210186545380,60436521934750,9662857324120,182297)
B5  (111203438664916,60432849886708,9662261108500,182285).              (7)
```

Every rational value has denominator `N`.  The identity checked interval by
interval is

```text
successor numerator
 =2*(marked-mass numerator)-(danger numerator).     (8)
```

## 3. All thirty defects have one sign

For each legal ordered graft pair `k!=l`, THM-3670's defect numerator is

```text
delta_(k,l)=S0+A_k-2B_l.                            (9)
```

The complete thirty-entry matrix is stored in the transcript.  Its largest
entry is still negative:

```text
max_(k!=l) delta_(k,l)=-12887674416812<0.           (10)
```

Thus all thirty ordered graft charts have nonzero successor defect.  The
unused omitted label gives four owner packets per chart, so all 120 packets
pass the rational test.  The canonical chart `(k,l)=(1,2)` is

```text
-18240864136612/N
 =-350785848781/967969917975060.                    (11)
```

In particular this control lies far from THM-3670's simultaneous
recirculation locus: neither the six `A_k` nor the six `B_l` are all equal.
The repeated first three values are a real symmetry of this row, not a
rounding artifact.

## 4. Exact energy invoice

THM-3674 applies because every legal pair-swap chart uses three distinct
sites.  For the canonical defect (11), its sharp tariff gives

```text
D >=123050711705006599185961
    /148212992112761067316289860457462400
   =8.3022891550147709424... * 10^-13,

E_dt>=123050711705006599185961
    /1926768897465893875111768185947011200
   =6.3863762730882853403... * 10^-14.              (12)
```

The largest defect magnitude in the thirty-entry ledger is

```text
20786137969714/N
 =799466844989/1935939835950120.                    (13)
```

It occurs at `(k,l)=(5,0),(5,1),(5,2)` and yields the stronger exact bounds

```text
D >=639147236236665754410121
    /592851968451044269265159441829849600
   =1.0780890850486303030... * 10^-12,

E_dt>=639147236236665754410121
    /7707075589863575500447072743788044800
   =8.2929929619125407924... * 10^-14.              (14)
```

Here `D` and `E_dt` are the normalized THM-3674 quantities for the selected
packet tensor.  Equation (14) is a certified energy floor on this non-cover
control, not a lower bound uniform over hypothetical covering rows.

## 5. Verification and scope

The companion pins the parent THM-2334 interval engine by SHA-256, reconstructs
all thirteen marked interval sets from half-open rational endpoints, evaluates
the danger comb by an independent prefix formula, and checks the complete
ledger and defect digests.  Normal and optimized runs produced the stored
transcript and ended in `PASS`.

This is a **typed non-cover positive control**.  It does not show that a
hypothetical covering row has any nonzero successor defect, does not exclude
the rigid thirteen-number pattern there, and does not preserve the separately
needed all-coordinate `91`-unit, visible-height, or terminal-phase sidecars.
LRC(14) remains open.

**QED (finite exact computation).**
