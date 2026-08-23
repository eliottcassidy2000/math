---
id: THM-3710
title: "LRC successor endpoint 182-skeleton and bad-prime collapse"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  In the typed THM-3672
  control, every far-blocker word boundary maps to the two-coset successor
  phase grid (14j+/-1)/182.  All thirty lawful charts contain its full 26
  points and the same 24-point exact-denominator-182 stratum, while their
  million-scale raw endpoint ledgers aggregate to only 42--45 phases.  The
  exact-182 and complementary endpoint vectors are nearly antiparallel in
  THM-3701's live scalar-line quotient (cosine -0.999183753...), explaining
  the small total leak.  Characteristic thirteen is bad for period 182:
  although the uniform full grid is a top-radical packet, deleting the two
  denominator-14 corners drops the support indicator to radical order zero,
  while all sixty actual signed **full-grid** chart fibers also have order
  zero.  After the denominator-14 corners are deleted, three exceptional
  exact-182 fibers instead have order two.  This refutes a support-only
  Ramanujan/Bockstein shortcut on the pinned non-cover control; it proves no
  genuine-cover statement and no LRC(14) conclusion.
source: codex-lrc14-20260822 / far-word endpoint pullback and THM-3672 ledger
audit: >
  PASS -- two independent audits reconstructed the rational successor-phase
  aggregation, all thirty endpoint ledgers, scalar-line split, and Dedekind
  value.  A separate characteristic-thirteen referee verified the radical
  algebra, all sixty full-grid order-zero fibers, and the corner-deleted
  histogram {0:57,2:3} with exactly the three stated order-two exceptions.
  Normal and optimized transcripts match the stored output and both raw-LF
  hashes match this metadata.  No cover-row or LRC(14) overclaim was found.
depends_on:
  - THM-3672-lrc-successor-mass-all-packet-positive-control
  - THM-3705-lrc-successor-bulk-endpoint-leak-tariff
  - THM-2041-frobenius-stability-of-exact-period-projectors
related:
  - THM-2058-primitive-phase-packets-and-deck-fan-intervals
  - THM-3706-lrc-dedekind-checksum-unbounded-bockstein-depth
  - THM-732-disc-v-bernoulli-edge-pair-dedekind-form-exact-certificates-far-element-tail
script: 04-computation/lrc_successor_endpoint_182_bad_prime_thm3710.py
output: 05-knowledge/results/lrc_successor_endpoint_182_bad_prime_thm3710.out
script_sha256: a555f2dbfd5ee342db8e5f7ed08b1b864fd24de98fa4baa89f0f40f094548688
output_sha256: a45658125b28bc56f80e4e767c4fa4c31e87efab54618e27e425e3b524542683
hash_basis: raw LF bytes
---

# THM-3710 -- the 182 skeleton does not force uniform radical depth

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.**  THM-3705 exposed the
signed endpoint vector as the small correction to the plain marked-mass bulk.
This theorem resolves its exact phase anatomy on the canonical THM-3672
positive control.  The phase map and bad-prime algebra are symbolic; all sign,
support, multiplicity, and distance claims are finite-exact statements about
that one typed non-cover row.

## 1. Endpoint ledger and fixed universe

Use THM-3672's exact constants

```text
p=13,       n=14,
R=p^2,      w_B=2p^5,      C=p w_B=2p^6.            (1)
```

Here `R` is the delayed-word clock, `w_B` the far/deep blocker in the word,
and `C` the successor frequency.  For each ordered chart `(k,l)` with
`0<=k,l<6` and `k!=l`, aggregate all endpoints of the signed marked-set
combination

```text
E_0+E_Ak-2E_Bl                                      (2)
```

by successor phase `theta={Cx}`.  Let `c_(k,l)(theta)` be the resulting
integer multiplicity.  With `P` the THM-3705 centred danger primitive, the
chart endpoint residual is exactly

```text
e_(k,l)=(1/C) sum_theta c_(k,l)(theta)P(theta).      (3)
```

The finite universe below is exactly these thirty charts reconstructed from
the pinned THM-3672 interval engine.  It is not a sample over arbitrary rows.

## 2. The far-word boundary map produces 182

Put `u={Rx}`.  Since `C/R=w_B/p=2p^4` is an integer, the successor phase of a
word endpoint does not depend on the choice of word branch.  A zero-shift
far-blocker boundary obeys

```text
w_B u=q+/-1/n.                                      (4)
```

Therefore

```text
{Cx}={(C/R)u}
     =(q+/-1/n)/p
     =(nq+/-1)/(np)       (mod 1).                  (5)
```

As `q` varies modulo thirteen, these phases form the full two-coset grid

```text
G_182={14j+/-1 mod 182:j in Z/13Z}.                 (6)
```

It has 26 points.  The residues `13` and `169` reduce to the threshold phases
`1/14` and `13/14`; the remaining 24 phases have exact denominator 182, with
numerators

```text
1,15,27,29,41,43,55,57,69,71,83,85,
97,99,111,113,125,127,139,141,153,155,167,181.      (7)
```

The same integer appears in the older Dedekind lane for a structural reason.
For `D=np+1=183`, one reciprocity step gives

```text
s(14,183)=-91/1098=-np/(12D),
-12D s(14,183)=np=182.                              (8)
```

Thus the previously observed scaled Dedekind integer `-12D s` is literally
the successor phase modulus obtained by pulling the far-word boundary through
the clock.  No modular-form interpretation is used.

## 3. Exact compression of all thirty control charts

Before phase aggregation, the weighted endpoint `l1` counts range over

```text
1435476..1461670.                                   (9)
```

After exact aggregation, every chart has only

```text
42..45 nonzero phases,
effective coefficient l1 41856..67762.              (10)
```

Every chart contains every point of the full grid (6), and every chart has
exactly the common primitive support (7).  The actual coefficients are not
uniform: there are thirteen coefficient profiles modulo thirteen across the
thirty charts.

For the canonical chart `(k,l)=(1,2)`, the three ledger sizes are

```text
(support,effective l1,naive l1)=(42,59408,1461670). (11)
```

Its exact residual split is

```text
e_total
 =11109404003/1129298237637570
 =+9.837440308275...*10^-6,

e_exact182
 =586695/6149354666
 =+9.540757231712...*10^-5,

e_nonexact182
 =-3716699972/43434547601445
 =-8.557013200884...*10^-5.                          (12)
```

The unaggregated absolute primitive bound is `942.3284...` times the true
residual.  Aggregating equal phases first improves this to `38.2999...`
times the true residual, still far from the signed answer.

The sign pattern in (12) holds on every chart.  The exact ranges are

```text
e_total:
  83117797033/9034385901100560
  .. 12230170063/1129298237637570,

e_exact182:
  78691/878479238 .. 598097/6149354666,

e_nonexact182:
  -61198645499/694952761623120
  .. -3420318877/43434547601445.                    (13)
```

Thus the small positive leak is uniformly a cancellation between a positive
exact-denominator-182 phase piece and a negative complementary piece on this
control.

## 4. Near-antiparallel cancellation in the live scalar quotient

Let `e^182,e^rest` be real 13-vectors of the individual endpoint residuals
before chart differencing, split by exact denominator 182, and let

```text
L_sc=span over the reals {(1,...,1)}.               (14)
```

This is THM-3701/3705's live scalar-mass line, not THM-3670's older
two-dimensional forward-mask kernel.  Exact projection gives

```text
A=dist(e^182,L_sc)^2
 =2004390109503/245794658253663815114,

B=dist(e^rest,L_sc)^2
 =21868283125702350430567
  /3139235715769406826669970473600,

I=<proj(e^182),proj(e^rest)>
 =-418382806563352927
  /55555643095377343171212960.                      (15)
```

Consequently

```text
dist(e^182,L_sc)       =9.030356701099...*10^-5,
dist(e^rest,L_sc)      =8.346326706245...*10^-5,
dist(e^182+e^rest,L_sc)=7.687253113389...*10^-6,

I/sqrt(AB)=-0.999183753026638....                   (16)
```

Both pieces are more than ten times farther from the scalar line than their
sum.  Separating exact phase denominators and applying absolute bounds before
recombining therefore destroys the cancellation that makes THM-3705 well
calibrated.

## 5. The bad-prime radical shortcut collapses

THM-2041's exact-period projector theorem assumes that the characteristic is
coprime to the period and proves that this hypothesis is sharp.  Here the
period is `182=14*13`, so characteristic thirteen is bad.  Put

```text
Alg=F_13[C_182],
J=rad(Alg)=(z-1)Alg,             z=g^14,             (16a)
```

where `g` is the canonical generator of `C_182`.  Then `J^13=0` and
`Alg/J` is the semisimple algebra `F_13[C_14]`.  In `Alg`,

```text
(z-1)^13=0,
sum_(j=0)^12 z^j=(z-1)^12,

sum_j (g^(14j+1)+g^(14j-1))
 =(g+g^(-1))(z-1)^12.                               (17)
```

The factor `g+g^(-1)` is a unit modulo `J`: equivalently,
`gcd(X^2+1,X^14-1)=1` in characteristic thirteen.  Hence the uniform full
grid has radical order exactly twelve, not merely membership in `J^12`.
This attractive observation does not survive either exact-denominator
attribution or actual multiplicity.  Removing the two denominator-14 corner
phases leaves twelve points in each sign fiber, whose augmentation is

```text
12=-1 mod 13.                                       (18)
```

Hence each primitive-support indicator has radical order zero.

More decisively, for an actual chart and sign define the **full-grid** fiber

```text
F_+/- (z)
 =sum_(j=0)^12 c_(k,l)({(14j+/-1)/182})z^j
   in F_13[z]/((z-1)^13).                           (19)
```

Its radical order is the least `r` for which the Hasse jet

```text
sum_j c_j binom(j,r)                                (20)
```

is nonzero modulo thirteen.  All sixty full-grid chart/sign fibers (with the
denominator-14 corners retained) have order zero.  For the canonical chart
their coefficient vectors modulo thirteen are

```text
plus  =(6,0,0,0,4,10,10,12,5,12,12,12,12),
minus =(7,6,0,0,9, 3, 3, 1,8, 1, 1, 1, 7).         (21)
```

Deleting the plus index `12` and minus index `1` produces the actual
exact-denominator-182 coefficient fibers.  Their order distribution is

```text
57 fibers of radical order 0,
 3 fibers of radical order 2:
   minus fibers (k,l)=(4,0),(4,1),(4,2).             (21a)
```

Thus corner deletion need not preserve depth zero: it leaves 57 fibers at
zero and raises three exceptional fibers to order two.

The bounded claim

```text
"the presence of the 182 grid forces every actual full-grid or
 corner-deleted exact-182 sign fiber to inherit a positive uniform
 13-primary radical/Bockstein depth"                         (22)
```

is therefore **REFUTED** on the exact control.  Even the uniform full grid
has zero endpoint-primitive contribution by antipodal oddness of `P`, so
support alone never determined the signed leak.

## 6. Preserved data, destroyed data, and the next test

The phase map preserves successor phase, exact denominator, antipodal sign,
chart multiplicity, and the thirteen-primary radical filtration.  Replacing
the ledger by its support or by a scalar Ramanujan/Dedekind statistic destroys
signed multiplicity, the two lower-period threshold atoms, endpoint owner,
the non-182 complement, bulk mass, and scalar-cover semantics.

The supports here are exact **phase-denominator basis supports**, as in
THM-2058; they are not THM-2041 exact-frequency character idempotents.
THM-2041 is used only for its sharp bad-prime/nonreduced boundary.  THM-2058
anticipated the accompanying loss: reduced phase-order packets require an
orientation/residue sidecar.  Here the 24 primitive phases are twelve
antipodal pairs, and their contribution depends on the oriented coefficient
differences `c(theta)-c(-theta)`.  THM-3706 supplies the p-adic analogue: one
Dedekind digit does not control the higher checksum tower.  These old threads
meet at the same missing coordinate -- signed multiplicity through a
bad-prime exact-phase-denominator split.

The actionable genuine-cover target is now narrower:

```text
retain, for both sign fibers,
  the thirteen signed Hasse jets,
  the denominator-14 corner attribution,
  the complementary endpoint ledger,
  and the endpoint owner/current;

then prove that cover semantics prevents the near-antiparallel
payment in (15), or that the bulk still exceeds it.                 (23)
```

Part (5)--(8) is symbolic for the typed clock and far boundary.  Everything
in (9)--(21a) is finite-exact on THM-3672's typed positive **non-cover** row.
No theorem here says a genuine cover has the same support, sign split,
scalar-line cosine, or Hasse profile.  It supplies no preselected triangle,
visible height, all-coordinate `91`-unit current, or owner intertwiner, and
does not prove LRC(14).

The standard-library-only companion pins THM-3672, its THM-2334 interval
engine, and THM-3705.  It rebuilds all thirteen marked sets, checks the
endpoint primitive independently against danger-prefix integration, audits
all thirty charts, and matches under normal and optimized execution.

```bash
python -B 04-computation/lrc_successor_endpoint_182_bad_prime_thm3710.py
python -B -O 04-computation/lrc_successor_endpoint_182_bad_prime_thm3710.py
```

**QED.**
