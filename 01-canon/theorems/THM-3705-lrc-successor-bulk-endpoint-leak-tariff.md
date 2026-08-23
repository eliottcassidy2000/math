---
id: THM-3705
title: "LRC successor bulk-endpoint leak tariff"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT.  Every THM-3670 successor mass is
  exactly 13/7 times its plain marked mass minus a signed rational endpoint
  residual.  The centred danger primitive has sharp range [-3/49,3/49], so
  an r-interval-piece representation has endpoint leak at most 6r/(49C) at
  successor frequency C.  Orthogonal projection away from the thirteen-mass scalar
  line and THM-3701 give a one-packet drift tariff from the excess of plain
  mass separation over endpoint leakage.  The THM-3672 control retains
  97.489 percent of its actual scalar-line distance under this reverse
  triangle gate.  The unsigned component bound is too weak even there;
  cover-specific signed endpoint cancellation remains open.  No LRC(14)
  conclusion is claimed.
source: codex-lrc14-20260822 / endpoint primitive composed with THM-3701
depends_on:
  - THM-3701-lrc-radial-successor-mass-gate-and-star-frame
  - THM-3672-lrc-successor-mass-all-packet-positive-control
related:
  - THM-732-disc-v-bernoulli-edge-pair-dedekind-form-exact-certificates-far-element-tail
  - THM-744-F-telescoping-refinement-of-Xi
script: 04-computation/lrc_successor_bulk_endpoint_leak_thm3705.py
output: 05-knowledge/results/lrc_successor_bulk_endpoint_leak_thm3705.out
script_sha256: 82faa434621b965a4d714b932c078698ec0bff00171f43d9172e25615c1534cf
output_sha256: 035538645e5184ad0bce3625a1d286611968ecd0c901bcd5b85b43e99cc3033b
hash_basis: raw LF bytes
---

# THM-3705 -- successor mass is bulk minus a signed endpoint leak

**PROVED + FINITE-EXACT + VERIFIED-EXACT.**  THM-3701 reduced the live
successor-mass failure locus to one scalar line.  This theorem separates the
input to that gate into a geometric bulk vector and an exact boundary ledger.
The resulting criterion is sufficient, not necessary, and does not yet
exclude the scalar line on a genuine cover.

## 1. The exact one-dimensional endpoint identity

Write

```text
d(y)=1_(||y||<1/14),
P(y)=integral_0^y (d(t)-1/7) dt,       0<=y<=1,      (1)
```

and extend `P` periodically.  On one period,

```text
P(y)= 6y/7,                 0<=y<=1/14,
      1/14-y/7,             1/14<=y<=13/14,
      6(y-1)/7,             13/14<=y<=1.            (2)
```

Thus `P(0)=P(1)=0` and sharply

```text
|P(y)|<=3/49.                                      (3)
```

Let `C` be a positive integer and let `E` be a disjoint union of `r`
half-open intervals `[a_j,b_j)` in the circle, represented in `[0,1]`.  The
chain rule away from finitely many endpoints gives

```text
integral_E d(Cx) dx
 =|E|/7+e_C(E),

e_C(E)
 =(1/C) sum_j [P({C b_j})-P({C a_j})],             (4)

|e_C(E)|<=6r/(49C).                                (5)
```

All endpoint conventions change only null sets.  When the interval endpoints
are rational, (4) is an exact rational finite sum.  The successor observable
therefore has the identity

```text
S_C(E):=integral_E (2-d(Cx)) dx
       =13|E|/7-e_C(E).                            (6)
```

Equation (5) is attained already by the one-interval witness

```text
E=[1/(14C),13/(14C)).                               (6a)
```

Indeed, `d(Cx)=0` almost everywhere on this interval, while
`|E|=6/(7C)`, so `e_C(E)=-6/(49C)` exactly.

## 2. The thirteen-coordinate scalar-line gate

Fix exactly the THM-3701/THM-3670 scope: one owner/word stratum whose
successor blocker is one of the two target blockers.  Its successor frequency
is the common integer `C=13c`.  Let

```text
E_0,E_A0,...,E_A5,E_B0,...,E_B5                    (7)
```

be the thirteen marked sets producing `S_0,A_0,...,A_5,B_0,...,B_5`.  Define
three vectors in `R^13`:

```text
M=(|E_0|,|E_A0|,...,|E_B5|),
e=(e_C(E_0),e_C(E_A0),...,e_C(E_B5)),
S=(S_0,A_0,...,A_5,B_0,...,B_5).                   (8)
```

Componentwise (6) gives the exact vector identity

```text
S=(13/7)M-e.                                       (9)
```

Let `L_sc=R(1,...,1)` and write `dist` for Euclidean distance.  Orthogonal
projection to `L_sc^perp` followed by the reverse triangle inequality gives

```text
dist(S,L_sc)
 >=Gamma_sc
 :=max(0,(13/7)dist(M,L_sc)-dist(e,L_sc)).         (10)
```

THM-3701, equation (22b), now implies that in every six-chart derangement
bank some lawful packet satisfies

```text
D_H >= Gamma_sc^2/474552,
E_dt>= Gamma_sc^2/6169176.                         (11)
```

In particular `Gamma_sc>0` forces a nonconstant successor marginal, positive
drift, and some exact `gcd(m,91)=1` nonzero deep-colour/nonzero-target fibre.
No phase-cone premise is used.

If the thirteen component counts are `r_i`, (5) supplies the cheaper but
usually much weaker sufficient condition

```text
dist(e,L_sc)<=||e||_2
 <=(6/(49C)) sqrt(sum_i r_i^2),                    (12)

(13/7)dist(M,L_sc)
 >(6/(49C)) sqrt(sum_i r_i^2).                     (13)
```

The useful object is therefore the signed endpoint vector `e`, not merely
the number of components.

## 3. Exact calibration on the THM-3672 control

For THM-3672's thirteen-entry ledger, with its common denominator already
included, the exact squared distances are

```text
dist(M,L_sc)^2
 =7267086384084568099667227
  /97444439258883015987041328374400,

dist(e,L_sc)^2
 =31351115173959958275463
  /530530835965029753707225010038400,

dist(S,L_sc)^2
 =6394584339188742740200237
  /24361109814720753996760332093600.               (14)
```

Numerically,

```text
dist(M,L_sc)             =0.00027308738005673034...,
(13/7)dist(M,L_sc)       =0.00050716227724821349...,
dist(e,L_sc)             =0.00000768725311338875...,
Gamma_sc                 =0.00049947502413482475...,
dist(S,L_sc)             =0.00051233924637549495.... (15)
```

Thus `Gamma_sc>0`, the endpoint leak consumes only about `1.516%` of the
scaled bulk distance, and (10) retains

```text
Gamma_sc/dist(S,L_sc)=0.9748912027886265....        (16)
```

Equation (11) alone gives the exact lower bounds whose right sides are
approximately

```text
D_H:  approximately 5.257069820261716*10^-13,
E_dt: approximately 4.043899861739781*10^-14.       (17)
```

This is a calibration of the new route, not the strongest available control
invoice: THM-3701's direct best chart gives a larger value.

For the canonical forward chart `(k,l)=(1,2)`, before division by
`N=50334435734703120`, (9) reads exactly

```text
bulk defect       =-124219914907348/7,
endpoint residual =   3466134048936/7,
successor defect  =      -18240864136612.           (18)
```

The sign is `successor=bulk-endpoint`.  In contrast, the unsigned
component-count estimate in (12) is `0.0081249209...` on this control, over
sixteen times the scaled bulk distance.  It does not even prove
`Gamma_sc>0`.  The strong result in (15) comes from signed endpoint
cancellation.

## 4. Dedekind-adjacent interpretation and exact frontier

The primitive `P` is a periodic piecewise-linear, first-Bernoulli-type
endpoint observable.  Its values at the rational LRC arrangement endpoints
form a finite signed rational ledger.  Squaring or Fourier-expanding such
endpoint data enters the generalized `B_2` edge-pair world of canonical
THM-732.  Corrected THM-744 supplies a related structural warning: endpoint
primitive terms can telescope while short signed residuals remain.

Neither comparison identifies (4) with one classical two-variable Dedekind
sum.  The exact covering-row target exposed by (10) is instead

```text
prove on every genuine scalar cover with an eligible low target that

(13/7)dist(M,L_sc)>dist(e,L_sc),                   (19)
```

or force an unsampled successor-table variation by another route.  A proof
of (19) needs cover-specific pairing among the signed endpoints.  Per-set
absolute bounds lose precisely that information.

The theorem preserves rationality, endpoint provenance, the successor scale,
and THM-3701's lawful packet transfer.  It does not preserve a preselected
frequency triangle, visible height, terminal phase, or all-coordinate
`91`-unit support.  It does not exclude the scalar line and does not prove
LRC(14).

The standard-library-only companion pins the THM-3672 and THM-3701 artifacts,
checks 593 exact primitive values and 15 independent interval identities,
reproduces (14) and (18), and matches under normal and optimized execution.

```bash
python -B 04-computation/lrc_successor_bulk_endpoint_leak_thm3705.py
python -B -O 04-computation/lrc_successor_bulk_endpoint_leak_thm3705.py
```

**QED.**
