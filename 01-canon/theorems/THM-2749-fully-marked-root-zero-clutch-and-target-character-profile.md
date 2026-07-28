---
id: THM-2749
title: "Fully marked root-zero clutch and target-character profile"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On every one of
  the fourteen canonical source-one rails, the equal-weight root-zero overlap
  remains a primitive unit after both endpoint present complements, both
  translated E3/full-two-target sections, and the delayed
  Q_(3,{1,2}) terminal fork are placed inside the coefficient integral.  The
  uniform lawful label (s,t)=(0,4) works on all fourteen rails.  Translation
  by 7/13^6 makes the two raw clock vectors equal; root normalization at 12
  and 1 changes that equality by exactly -1.  Freezing the source section at
  clock one gives the exact rail-eight target window t=3,...,11; its window
  polynomial is an integral cyclotomic unit.  In the distinct clock-coindexed
  family the raw target scan is supported on t=2,...,11, while its explicitly
  transported-cylinder-attached linear and bilinear profiles have support
  t=3,...,11 and all twelve primitive C13 characters.  This is a partial
  semantic clutch, not a THM-2334 endpoint current, row exclusion, or LRC(14).
source: root/fully-marked-root-zero-target-profile-2026-07-28 + coordinate-first-audit/fully-marked-root-zero-clutch-2026-07-28
audit: >
  root-zero-two-sided-clutch-hostile-audit-2026-07-28 (independent fixed-window
  interval, 239-piece coefficient, unit, Bezout, and hash audit: ACCEPT);
  root-2026-07-28 (independent all-rail, target-bank, and hostile replay);
  semantic-clutch-referee-2026-07-28 (independent carrier rebuild,
  lifted-interval integration, determinant, and witness audit: ACCEPT)
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2742-full-two-target-present-sheet-deepest-source-semantic-current
  - THM-2744-relative-present-unit-repair-and-root-zero-overlap-clutch
related:
  - THM-830-b3-deletion-deck-mirror-current-calculus
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
  - THM-2750-arm-blind-clutch-no-go-and-minimal-marked-leakage
script: 04-computation/lrc14_fully_marked_root_zero_target_profile_thm2749.py
output: 05-knowledge/results/lrc14_fully_marked_root_zero_target_profile_thm2749.out
script_sha256: 12f150dc8e0fc543cc36fafaa2b84dd57a2dde6e40ce3cbadd8d057817bce3dc
output_sha256: 72fed42be733aca63fc0ccd0a907eadcb02d224c8832d2cd5c42208b34a18048
secondary_script: 04-computation/lrc14_fully_marked_root_zero_clutch_thm2749.py
secondary_output: 05-knowledge/results/lrc14_fully_marked_root_zero_clutch_thm2749.out
secondary_script_sha256: 93b46b2701db8f72d00fa2ae131f9d9afd3200f32998959af3bb2e1fa2f56841
secondary_output_sha256: 04ae77696e91e71be0c4aa4981fdbec720d5b4040d4c29610024f2f9d393271b
hash_basis: LF-normalized bytes
---

# THM-2749 -- the root-zero clutch survives all common markings

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2744 left one precise test: put the semantic fork and a lawful target
section inside the overlap integral, rather than checking them only at one
point.  A first one-sided implementation gave unequal coefficients.  The
missing coordinate is which endpoint owns the mask.  Pulling both endpoint
masks to one common physical carrier restores translation covariance and
gives a uniform all-rail unit.  The result is still coefficient-side and
clock-indexed; no endpoint current is being renamed.

## 1. The common fully marked carrier

Use the canonical period and scales

```text
p=13,     R=13^6=4826809,     S=R/13=371293,
T=297836897838480,            tau=7/R.                 (1)
```

Let `E3` be the exclusive deepest-source cell, let
`Q=Q_(3,{1,2})` be the ordinary deepest terminal fork, and let

```text
S_(ell,s,t)=E3 intersect F_(ell,s,t)                   (2)
```

be THM-2742's full two-target source section.  For rail `j`, pull the target
rail profile back by `tau`, retain only its equal-weight overlap with the
source profile, and write `P_ell^c` for the fixed-label-`7` relative-present
complement.  In source coordinates the common carrier is

```text
C_(j,ell)^(s,t)
 =rho_j intersect T_tau^(-1)rho_j
  intersect P_ell^c intersect T_tau^(-1)P_ell^c
  intersect H^R_12 intersect T_tau^(-1)H^L_1
  intersect S_(ell,s,t) intersect T_tau^(-1)S_(ell,s,t).     (3)
```

Weights on the two copies of `rho_j` are required to agree on `(3)`.  The
target carrier is literally `T_tau C_(j,ell)^(s,t)`.  Finally put the delayed
sector-zero, `h=6`, `kappa=1`, and `D^(-6)Q` digit inside the delayed-carry
functional.  Let

```text
A^-_(j,ell)(s,t) = source-carry-12 coefficient on (3),
A^+_(j,ell)(s,t) = target-carry-6 coefficient on T_tau(3).  (4)
```

The index `ell` occurs in both the carrier and delayed digit.  Thus the seven
entries in `(4)` are **counterfactual clock fibres**
`E3 intersect F_(ell,s,t)`, not seven iterates of one point or one fixed set.

There is a second, inequivalent specialization used below: freeze the source
section at `E3 intersect F_(1,0,t)` while the delayed clock and present
complement vary.  The frozen-section and clock-coindexed banks must not be
identified; their different raw `t`-supports are an exact diagnostic of the
coordinate they retain.

## 2. General marked-transport lemma

The equality in `(4)` is structural once the carrier is common.

> **Marked-transport lemma.**  Let `C` be any weighted source carrier for
> which translation by `tau` gives the target carrier with identical weights.
> Restrict source and target to carry cells `12` and `6`.  Then their
> integrals agree against every retained predicate depending only on
> `{Rx}`.

Indeed,

```text
R tau=7,       S tau=7/13,       c3 tau*182=196=14 mod182. (5)
```

The first identity fixes `{Rx}` and hence every delayed target predicate.
The second maps the carry-`12` open cell bijectively to carry `6`, including
the wrap at one.  The third maps the source deep overlap `(169,181)/182` to
the target overlap `(1,13)/182`.  Translation preserves length and the
assumed weight profile, so change of variables proves the lemma.  Exact
enumeration is needed below only for positivity, content, supports, and unit
determinants.

This lemma also explains the failed one-sided probe: a source mask `S` and a
separately imposed target mask `S` need not be translates.  The common mask
in `(3)` is `S intersect T_tau^(-1)S`; its translate is
`S intersect T_tau S`.

## 3. One uniform section gives fourteen units

Take

```text
(s,t)=(0,4).                                             (6)
```

For every source-one rail `j=0,...,13`, exact integration gives

```text
(A^-_(j,ell)(0,4))_(ell=0)^6
 =(A^+_(j,ell)(0,4))_(ell=0)^6 !=0.                    (7)
```

The exact gcd of all coefficients in this selected fourteen-rail bank is

```text
G_04=413915261760=7224*57297240,          v_13(G_04)=1. (8)
```

Here `g=57297240` is only an inherited broadcast divisor, not the intrinsic
gcd of this subbank.  Likewise the older content `26` is an inherited
lattice normalization.  The ratios

```text
G_04/g=7224=9 mod13,       G_04/26=15919817760=2 mod13 (9)
```

are units after their common `13`-valuation is removed, so all three content
choices give the same unit/nonunit truth.

Divide `(7)` by `G_04`, normalize the source by root `12` and the target by
root `1`, reduce modulo `13`, and subtract the seventh clock coordinate.
The exact six-coordinate profiles and multiplication determinants in
`F_13[C_7]/(Phi_7)` are:

| rail | raw clock support | source/root 12 | target/root 1 | determinant |
|---:|:---:|:---|:---|---:|
| 0 | `1` | `(0,10,0,0,0,0)` | `(0,3,0,0,0,0)` | 1 |
| 1 | `6` | `(1,1,1,1,1,1)` | `(12,12,12,12,12,12)` | 1 |
| 2 | `2,3` | `(0,0,9,2,0,0)` | `(0,0,4,11,0,0)` | 3 |
| 3 | `5` | `(0,0,0,0,0,4)` | `(0,0,0,0,0,9)` | 1 |
| 4 | `2,3` | `(0,0,11,11,0,0)` | `(0,0,2,2,0,0)` | 12 |
| 5 | `4,6` | `(9,9,9,9,3,9)` | `(4,4,4,4,10,4)` | 8 |
| 6 | `2,3` | `(0,0,5,5,0,0)` | `(0,0,8,8,0,0)` | 12 |
| 7 | `5,6` | `(10,10,10,10,10,2)` | `(3,3,3,3,3,11)` | 3 |
| 8 | `1,2,3` | `(0,2,5,4,0,0)` | `(0,11,8,9,0,0)` | 12 |
| 9 | `5` | `(0,0,0,0,0,5)` | `(0,0,0,0,0,8)` | 12 |
| 10 | `1,3` | `(0,3,0,1,0,0)` | `(0,10,0,12,0,0)` | 1 |
| 11 | `5` | `(0,0,0,0,0,5)` | `(0,0,0,0,0,8)` | 12 |
| 12 | `2` | `(0,0,4,0,0,0)` | `(0,0,9,0,0,0)` | 1 |
| 13 | `5,6` | `(10,10,10,10,10,0)` | `(3,3,3,3,3,0)` | 1 |

Every determinant is nonzero.  Moreover every source profile is the negative
of the target profile.  This sign is forced simply because

```text
12^(-1)=-1,                         1^(-1)=1 mod13.    (10)
```

It is a transition in the local coefficient fibre.  It is not a global
physical reflection, a `C2` action, or equality of root-normalized rows.

For orientation, rail `8` has the raw vector

```text
(0,
 339633525654239542165440,
 750593782703678965571520,
 719200126392878704654080,
 0,0,0),                                                (11)
```

with the source/target reduced profiles and determinant shown in the table.

## 4. The open cylinder and the typed target profiles

The explicit adjacent pair is

```text
q =47850889647341/100360982066072,
q'=47851035194197/100360982066072=q+7/R,               (12)
r_cyl=1/100360982066072.                               (13)
```

It lies in clock fibre `ell=1`, rail `8`, and the complete `(0,4)` carrier.
The source and target section slacks are respectively `56447 r_cyl` and
`8854585 r_cyl`; both weighted-carrier slacks are `56447 r_cyl`.  The delayed
fork and digit slacks are exactly `r_cyl`.  Hence every point with
`|delta|<r_cyl` stays in the complete marked carrier and keeps

```text
E3 -> D^6 Q_(3,{1,2}),       carry 12 -> carry 6.      (14)
```

The cylinder is open: its two endpoints lie on the delayed-fork boundary and
are not claimed to satisfy the strict predicate.

### The frozen-section nine-window

First freeze `E3 intersect F_(1,0,t)` for all seven delayed clocks, keep the
same two-sided pullback, and stay on rail `8`.  Exact integration gives

```text
V_t^-=V_t^+=0                              for t=0,1,2,12,
V_t^-=V_t^+=(0,C_0,C_0,C_0,C_0,C_0,C_0)   for t=3,...,11,

C_0=339633525654239542165440.                         (F1)
```

Every supported label has common grid mass `6320326320`.  Moreover
`v_13(C_0)=1` and `(C_0/26) mod13=9`.  After content-`26` division and root
normalization, the source/root-`12` and target/root-`1` classes are the
constants `9` and `4=-9`, both with determinant one.
This is the local mirror sign of the seam calculus in THM-830 and the central
`-1` shadow suggested by THM-2716; it is not the latter theorem's full `C4`
transporter groupoid.

The target window is therefore

```text
W(u)=u^3+u^4+...+u^11.                                (F2)
```

It is not merely nonzero on all primitive characters.  With

```text
V(u)=u^2+u^6+u^10,
W(u)V(u)-1=(u^9+u^5+u-1)Phi_13(u),                    (F3)
```

so `W` is an integral unit in `Z[u]/(Phi_13)` with positive three-term
inverse `V` and resultant norm one.  In `Z[C_13]`, its circular product is
`(3,2,2,...,2)`, hence a delta after quotienting the uniform target-null
line.  This is a coefficient recombination of lawful `t`-sections, not three
physical translates of one Boolean packet.  In characteristic thirteen,
`Phi_13=(u-1)^12`, `W(1)=9`, and `V(1)=3`; both multiplication determinants
are one.

### The clock-coindexed target scan

Now scan `s=0` and all `t in F_13` in the distinct clock-coindexed family of
`(2)--(4)`.  Its **raw** rail-eight integral is
positive for

```text
t=2,3,...,11.                                          (15)
```

Thus a bare “nine-label support” claim is ill-typed: the frozen-section raw
bank `(F1)` has nine labels, while the clock-coindexed raw bank has ten.  The
former unproved stub is superseded by these two typed statements.  A second
nine-label object inside the clock-coindexed bank is the explicit cylinder
attachment.  Let `epsilon_t=1` exactly when both open cylinders in
`(12)--(13)` lie in the `ell=1` common section `F_(1,0,t)`, and zero otherwise.
Exact interval slack gives

```text
support(epsilon)={3,4,...,11}.                          (16)
```

Define the linear pushforward and the clock-paired cross coefficient

```text
B_t=epsilon_t G_04^(-1) sum_ell A^-_(8,ell)(0,t),

C_t=epsilon_t G_04^(-2)
        sum_ell A^-_(8,ell)(0,t) A^+_(8,ell)(0,t).      (17)
```

The clock index makes `C_t` a typed bilinear comparison; raw equality makes
it a sum of squares, but that equality is used only after the two sides have
been constructed.  The exact profiles, in order `t=0,...,12`, are

```text
B=(0,0,0,
   4371492433154,4371492433154,4371492433154,
   4371492433154,4371492433154,
   2633938414646,2633938414646,
   2558092802727,2558092802727,
   0),                                                   (18)

C=(0,0,0,
   6980796083273674034188354,
   6980796083273674034188354,
   6980796083273674034188354,
   6980796083273674034188354,
   6980796083273674034188354,
   3961702116040374827642290,
   3961702116040374827642290,
   3692377863640893849986025,
   3692377863640893849986025,
   0).                                                   (19)
```

Both have exact support `(16)`.  If `zeta` is a primitive thirteenth root,
then for every `k=1,...,12`,

```text
sum_t B_t zeta^(kt) !=0,              sum_t C_t zeta^(kt) !=0. (20)
```

This second proof needs no numerical Fourier estimate.  A rational polynomial of degree
at most twelve vanishing at `zeta^k` is a rational multiple of
`Phi_13=1+X+...+X^12`; all thirteen coefficients would therefore be equal.
Multiplication by `k` merely permutes the coefficients.  The zero/positive
profiles `(18)--(19)` are nonconstant, proving all twelve assertions.

These are `C13` characters of the lawful second target label.  They are not
the seven clock characters, the deepest-target characters of THM-2742, or a
relation-address Fourier current.

## 5. Sharp controls and what is load-bearing

Three controls separate support from covariance.

1. At `t=0`, `E3` and the unshifted `c3`-safe target factor are disjoint, so
   both seven-clock vectors are exactly zero.
2. The clock-coindexed raw `t=2` vector is positive,

   ```text
   (0,0,750593782703678965571520,
       719200126392878704654080,0,0,0),                 (21)
   ```

   but the transported cylinder `(12)` has no `t=2` label, and the frozen
   bank `(F1)` is zero there.  This is the minimal witness separating the two
   target-bank conventions and raw support `(15)` from attached support
   `(16)`.
3. If each endpoint retains only its own `(0,4)` section, rather than the
   source/translated-target intersection in `(3)`, rail `8` gives

   ```text
   source=(0,339633525654239542165440,
            750593782703678965571520,
            722054095148406001101120,0,0,0),

   target=(0,345341652135823400016960,
            756301720214733558465600,
            724908063903933297548160,0,0,0).            (22)
   ```

   They are unequal.  The earlier exact
   `lrc14_semantic_root_zero_clutch_refinement_probe_20260728.py` found the
   same mechanism for its one-sided `U_(0,3)` mask: both endpoint vectors
   remained units but their coefficients sheared.  Its normalized constants
   are `5` and `9`, so the target/source gain is `9/5=7=6/12 mod13`, the carry
   ratio rather than the mirror sign.  That scalar cannot encode the additive
   odometer class: `Hom(C_13,F_13^*)=0`.  Equation `(22)` is the self-contained
   hostile in the present companion.

Therefore the common pulled-back section is load-bearing.  The audit does
not claim every individual present/chart cut remains separately essential
after all the other cuts; some become redundant on particular rails.

## 6. The remaining map is to a fixed THM-2334 triangle

The source object now available is a nonnegative family in `Q[C_7]` with
sidecars

```text
(rail j, clock ell, carries 12/6, roots 12/1,
 target label (s,t), chart tau, normalized sign -1).    (23)
```

THM-2334's target object is different: the complex exact-address orbit
coefficient

```text
C(r;X,m),                 Y=X+m w_c,                    (24)
```

and THM-2625 further refines an endpoint allocation to
`J(l,r)` and determinant sectors `S*(q,Delta)`.  The missing map must insert
the idempotent of `(3)` into one fixed `R=13^6` marked triangle while
retaining the exact relation address before the Radon quotient.

The loss ledger is asymmetric:

| direction | data preserved | data currently lost |
|---|---|---|
| `(23)` to `(24)` | semantic word, target label, translation, chart sign, one nonzero clock coefficient | Fourier allocation `(u,beta,v)`, exact address `r`, triangle `(X,m,Y)`, endpoint determinant sector |
| `(24)` to `(23)` | exact address orbit and complex current | rail, positivity, carry, private root, clock fibre, physical adjacent-chart cylinder |

The cheapest decisive test is therefore:

1. choose one `R=13^6` THM-2334 marked triangle;
2. insert the `(0,4)` carrier idempotent before expanding its endpoint
   factors;
3. retain the exact relation residue and one left/right allocation;
4. apply the THM-2625 dual endpoint transform without first summing over the
   address;
5. test whether one common nonzero `(q,Delta)` sector is transported with the
   normalized sign `-1`.

This is a fixed-triangle transplant, not another marginal-support census.

## 7. Reproduction and independent audit

Run

```bash
python3 04-computation/lrc14_fully_marked_root_zero_target_profile_thm2749.py
python3 -O 04-computation/lrc14_fully_marked_root_zero_target_profile_thm2749.py
python3 04-computation/lrc14_fully_marked_root_zero_clutch_thm2749.py
python3 -O 04-computation/lrc14_fully_marked_root_zero_clutch_thm2749.py
```

Both modes of both companions byte-match their declared outputs.  The
companions use exact integer and rational interval arithmetic, explicit
exceptions, and no truth-bearing `assert`.  The frozen-window companion pins
seven direct dependencies and independently compares its factored scan with
the full Boolean construction at `t=3`.  The all-rail companion pins its two
direct dependencies and checks the common-carrier translation and carry
descent on every coefficient.  A
second route directly enumerates the lifted intersections

```text
R[a,b] intersect (nT+[u,v])                            (25)
```

on rail `8`, without the delayed-prefix antiderivative; its hit counts are
`(0,37128,82056,78624,0,0,0)` and reproduce `(11)`.

An independent replay checked the frozen-window script in both modes against
its stored bytes, rederived `(F3)`, and audited the fixed/coindexed typing.
The independent referee rebuilt `E3`, the target sections, rail/chart
overlap, carry cells, and delayed fork without importing the clutch helper.
It directly enumerated `(25)`, recovered the carrier/carry piece counts
`239/119`, `526/263`, and `504/252` on the three live rail-eight clocks,
obtained determinant `12`, and replayed the open witness and both hostiles.
The `(0,3)` and `(0,4)` rail-eight common carriers are identical, and the
main companion separately performs the direct `(0,4)` integration.  A second
independent replay checked all fourteen unit rows, both target profiles, all
twelve character certificates, and normal/optimized equality.

## 8. Boundary ledger

```text
PROVED HERE:       general common-carrier marked-transport lemma;
                   one uniform (0,4) all-fourteen equal-weight unit bank;
                   exact selected content and normalized sign -1;
                   one strict fully marked open cylinder;
                   frozen-section raw support 3..11 and cyclotomic unit W;
                   clock-coindexed raw t-support 2..11;
                   attached linear/bilinear support 3..11 and all twelve
                   primitive target characters.

DESTROYED:         equality under separate one-sided semantic masks;
                   any unqualified identification of the two target banks;
                   edge/root label under physical rechart.

NOT CONSTRUCTED:   a global clutch action or arm-dependent transition;
                   a fixed THM-2334 triangle/address transplant;
                   endpoint current or pointwise transition amplitude;
                   row exclusion or LRC(14).                           (26)
```

The residual LRC ledger remains `165`.  QED.
