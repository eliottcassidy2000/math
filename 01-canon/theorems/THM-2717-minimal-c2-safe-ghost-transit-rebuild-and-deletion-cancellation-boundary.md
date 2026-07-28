---
id: THM-2717
title: "Minimal c2-safe ghost-transit rebuild and deletion-cancellation boundary"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED; AWAITING FINAL IMMUTABLE REPLAY/PROMOTION.  Deleting the
  failed c2 target gate and reinserting its safe complement rebuilds all 162
  private-root rails with global content 26.  The two forced THM-2706 ghost
  rows are global primitive units.  On the two forced slices, the exact c2
  truth-idempotent split gives an orientation-free finite-etale two-sheet
  coefficient carrier and exposes
  twenty reverse anti-diagonal cancellations; the forced-row quotient has an
  exact central-C4 shadow but supplies no physical or semantic arm selector.
source: root/minimal-c2-safe-ghost-transit-2026-07-28
depends_on:
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2698-central-half-odometer-full-local-cycle-and-semantic-sidecar-boundary
  - THM-2706-relative-segal-macro-cycle-and-minimal-ghost-midpoint-completion
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
  - THM-2720-unshifted-deepest-source-present-packet-global-disjointness
script: 04-computation/lrc14_minimal_c2_safe_ghost_transit_rebuild_thm2717.py
output: 05-knowledge/results/lrc14_minimal_c2_safe_ghost_transit_rebuild_thm2717.out
script_sha256: 81de8ea40fc00b05892c1430937d43dfb57112ed5b8a767282a37177ac4f94f2
output_sha256: bf7caaa73eaf2208cee4e27b0d09671bb34229458c99c7887004657001304437
hash_basis: LF-normalized bytes
---

# THM-2717 -- the minimal `c2`-safe ghost-transit rebuild and the deletion-cancellation boundary

## Status and scope

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED; AWAITING FINAL IMMUTABLE REPLAY/PROMOTION.**  This statement changes the
delayed middle-object grammar.  It rebuilds all 162 canonical source-one
private-root rails after deleting exactly the failed present `c2` gate and
then inserting its safe complement.  The resulting integer lattice has global
content `26`, and the two THM-2706 ghost phases carry global primitive-unit
profiles.  This is not an old THM-2640 sector, a THM-2305 terminal word, or a
canonical endpoint current.  No physical endpoint-current-to-repair cospan,
live-row transfer, row decrement, or LRC(14) conclusion is asserted.

## 1. The three exact middle grammars

Write `d_a` for the strict danger indicator of speed `a`, and
`g_a=1-d_a` for its safe complement.  Let `d_(c1,ell)` denote the translated
`c1` danger gate at source clock `ell in C7`.  For the delayed half digit

```text
d=2h+kappa,                 h in F13, kappa in {0,1},
I_d=[d/26,(d+1)/26),                                      (1)
```

put

```text
U=d_H product_(i=1)^5 g_(q_i) g_(c3),                    (2)

Qdel_(ell,h,kappa)=I_d U g_(c1,ell),                     (3)
Qsafe_(ell,h,kappa)=Qdel_(ell,h,kappa) g_(c2),           (4)
Qold_(ell,h,kappa)=Qdel_(ell,h,kappa) d_(c2).            (5)
```

Thus `Qdel` removes only the failed `c2` target gate from the old local
grammar, `Qsafe` is the `c2`-safe repair of the guard-danger word, and `Qold`
is its old guard-danger/`c2`-danger part.  Since
`g_(c2)+d_(c2)=1` away from the finite strict-open
seams, the words form the disjoint Boolean partition

```text
Qdel = Qsafe disjoint_union Qold                         (6)
```

almost everywhere.  Every row construction below is linear in its middle
word, so `(6)` gives the coefficientwise integer identity

```text
Adel = Asafe + Aold.                                     (7)
```

This identity is the reason the deleted-gate control must not be conflated
with the safe repair.  Projection to the safe arm is available only while the
`c2` truth label is retained.  After forgetting that label the coefficient
map is addition, and primitivity is not preserved by that addition.

## 2. Complete private-root rebuild

Use the THM-2640 carrier, but replace its delayed sector by one of `(3)--(5)`.
For every one of the 162 canonical source-one rails, retain

```text
(edge e, predecessor carry c, future half digit kappa,
 delayed digit h, future clock ell5).                    (8)
```

The exact prefix masses in every grammar are divisible by thirteen, so the
THM-2640 carry descent is integral.  With `d=2h+kappa`, every nonzero row is
still supported at the unique physical root

```text
r=2c+floor(d/13)+1_(e=left)             mod 13.          (9)
```

All other roots vanish before cyclotomic reduction.  For the predicted
nonzero root, write the seven-clock integer profile as

```text
A=(A_0,...,A_6).                                           (10)
```

For a grammar with global content `G`, normalize

```text
Y_ell=(A_ell/G) r^(-1)                    mod 13,
Y(z)=sum_(ell=0)^6 Y_ell z^ell in F13[z]/(Phi_7).         (11)
```

The profile is called primitive-unit when multiplication by `Y` on the
six-dimensional cyclotomic quotient has nonzero determinant.  Subtracting
the common coefficient `Y_6` before taking the determinant is the standard
representative of `(11)` and does not change the class modulo `Phi_7`.

## 3. The global safe lattice and the deletion hostile

The exhaustive four-shard rebuild over all 162 rails gives

```text
                         deleted c2 gate       c2-safe repair
global content G                 26                    26
v_13(G)                           1                     1
shard contents              26,26,26,26           26,26,26,26
positive raw entries             59,424                59,424
nonzero seven-clock profiles     23,872                23,872
primitive-unit profiles          19,689                21,134
nonunit profiles                  4,183                 2,738. (12)
```

For each grammar the referee checks `707,616` carry partitions and
`9,199,008` singleton-root slots.  In particular the provisional factor
`13^5` was an artifact of the much broader whole-circle control, not the
minimal repair lattice.

Neither column is a uniform-unit atlas.  The safe insertion recovers `1,445`
additional primitive profiles, but `2,738` safe profiles remain nonunits.
The deleted-gate column is nevertheless a sharper hostile: having the correct
global content does not recover the safe summand.  The reverse slice below
proves that the inference

```text
(safe-repair unit and old-danger unit)
  ==> deleted-gate unit                                   (13)
```

is invalid even before physical-current typing is considered.  The global
paired cross-tab gives the complete stronger statement:

```text
(deletion unit, safe unit, old-danger unit)       profiles
(0,0,0)                                             2,618
(0,1,1)                                             1,565
(1,0,1)                                               120
(1,1,0)                                               155
(1,1,1)                                            19,414. (13a)
```

Thus `1,565` profiles witness failure of unit preservation under the
codiagonal, while `120` prove that a deletion unit alone does not certify a
safe unit.  These are paired claims, not inferences from the marginal counts
in `(12)`.

## 4. Exact cancellation on the two forced ghost slices

The forward ghost slice is

```text
(j,h,kappa)=(9,0,1),                                     (14)
```

and the reverse ghost slice is

```text
(j,h,kappa)=(2,12,0).                                    (15)
```

Each slice has 22 nonzero `(edge,carry)` profiles.  Their exact local
contents and unit counts are

```text
                       content   v_13   profiles   units
forward deletion       2,225,652   1       22       22
forward safe           2,225,652   1       22       22
forward old danger   472,491,719,940,240  1  22       22

reverse deletion           1,092   1       22        2
reverse safe               1,092   1       22       22
reverse old danger    16,558,902,360  1    22       22.   (16)
```

The two reverse deletion units are exactly the left and right carry-zero
profiles.  On each of the other twenty reverse profiles the safe and old
danger summands are separately primitive units, while their sum is zero in
the cyclotomic quotient.  This assertion uses the **common global
normalization by 26** in all three arms.  The local contents in `(16)` are
used only to census unitness; their differing residual unit scalars are not
added.  Divisibility of the old arm by `26` follows from the raw identity
`(7)` and the two `26`-divisible global banks.  Thus the forgotten aggregate
can destroy the actual safe unit by exact mod-thirteen cancellation.

All six local contents in `(16)` have the same `13`-adic valuation as the
global content `26`; their residual scale factors are units modulo thirteen.
Consequently the zero/nonzero unit classification in `(16)` is unchanged
when the profiles are normalized in the global lattice.

### 4a. The orientation-free truth-idempotent carrier

Let

```text
K=F13[z]/(Phi_7),
K[e]/(e^2-e) = K x K,                                  (16a)
```

where `e` records the old-`c2`-danger arm and `1-e` the safe arm.  The
labelled profile is `(Asafe,Aold)`.  It is a unit in `(16a)` if and only if
both arm profiles are units.  The sheet swap exchanges the two components,
and its invariant multiplicative norm is

```text
N_e(Asafe,Aold)=Asafe*Aold in K.                        (16b)
```

More precisely, the two elementary symmetric invariants

```text
s=Asafe+Aold=Adel,              p=Asafe*Aold,
Delta=s^2-4p=(Asafe-Aold)^2                            (16c)
```

recover the unordered pair as the roots of `T^2-sT+p`.  Whenever `Delta` is
a unit, its two orderings form a finite etale `C2` torsor; a choice of square
root gives

```text
Asafe=(s+sqrt(Delta))/2,       Aold=(s-sqrt(Delta))/2.   (16d)
```

On **all 22** forward and all 22 reverse profiles, the labelled pair and its
norm are units.  On the twenty non-carry-zero reverse profiles one has the
stronger exact identity

```text
Adel=0,              Aold=-Asafe,
Asafe,Aold in F13^* z^2,                                (16e)
```

while every norm and every discriminant has determinant `1`.  Thus the
orientation-free unordered two-sheet object remains finite etale even where
the codiagonal vanishes.  What is missing is precisely the sign/orientation
identifying the `g_(c2)` root.  This is a genuine coefficient carrier which
does not choose the safe arm; neither the product nor its square-root choice
is yet a physical quadratic cross-current.  No symmetry of physical endpoint
currents realizing the abstract sheet swap is asserted.

## 5. The two exact safe ghost points

Let `R=13^6=4,826,809`.  The THM-2706 physical endpoint points are

```text
x0=39123022/82055753,       {R x0}=4/17,
x1=41305372/82055753,       {R x1}=13/17.                (17)
```

The forced middle points can be chosen as

```text
xf=41513423/82055753,       {R xf}=1/17,
   (j,c,h,kappa,e,r)=(9,7,0,1,left,2),

xr=38400313/82055753,       {R xr}=16/17,
   (j,c,h,kappa,e,r)=(2,0,12,0,left,2).                  (18)
```

Both points lie strictly in `Qsafe`; neither lies in `Qold`.  Every rail,
dynamically selected present factor, predecessor carry, half digit, private
half-tooth, and root test is strict.  With macro gains

```text
Kf=4,472,391,                 Kr=1,956,127,               (19)
```

the exact affine subdivisions are

```text
13*1,485,215+4,471,832 = Kf            mod R,
13*4,460,044+1,897,263 = Kr            mod R.             (20)
```

The first split has source/middle bare roots `(7,2)` and the second has
`(7,12)`.  The corresponding safe profiles are primitive units after
normalization by the **global** content `26`.  Thus the two missing
THM-2706 midpoints do admit exact changed-grammar, `c2`-safe, private-root
transit rows.  This is the positive theorem; no deletion-only inference is
used.

### 5a. Exact profile symmetry and the quarter-turn local system

In the common content-`26` lattice, after multiplying by the inverse private
root as in `(11)`, the two safe profiles are

```text
Yf=(7,0,0,0,2,10,0),              det(Yf)=11,
Yr=(0,0,3,0,0,0,0),               det(Yr)=1.             (20a)
```

The action generated by a nonzero `F13` scalar, cyclic clock shift, and
inversion `z->z^(-1)` has orbit sizes `168` and `84`, respectively, and the
two orbits are disjoint.  Thus the rows are not the same unbased label up to
the obvious clock, scale, and Frobenius symmetries.

Choose the first quadratic factor as the `(1,6)` character pair.  Then

```text
Phi_7=(z^2+3z+1)(z^2+6z+1)(z^2+5z+1),                  (20b)
```

where the factors correspond to inverse pairs `(1,6),(2,5),(3,4)`.  Writing
a component as `a+bz`, the exact components of `(20a)` are

```text
pair             Yf       order(Yf)       Yr       order(Yr)
(1,6), m=3      (6,1)        168         (10,4)        21
(2,5), m=6      (1,3)         56         (10,8)        21
(3,4), m=5      (4,2)        168         (10,11)       21. (20c)
```

Let `q=Yr/Yf` in `F13[C7]/(Phi_7)`.  Its three components and orders are

```text
(10,11), order 8;       (2,6), order 168;
(11,4), order 56.                                      (20d)
```

Hence `q` has order `168`, but one common exponent produces the exact
quarter-turn scalar

```text
q^42=5,                q^84=-1,                q^126=8=-5,
5^2=8^2=-1                    in F13.           (20e)
```

In fact the exact scalar intersection and quotient are

```text
<q> intersection F13^* = <q^42>={1,5,-1,8} = C4,
<q>/<q^42> = C42.                                         (20f)
```

This is a precise THM-2716 bridge and boundary.  The ordered pair realizes a
common square root of minus one across all three Frobenius inverse pairs,
matching the mod-17 quarter turn.  Reversing the order replaces `q` by
`q^(-1)` and swaps `5` with `8`.  There are exactly two orientation choices
identifying THM-2716's abstract macro `C4` (equivalently the subgroup
generated by `13` modulo `17`) with the central scalar `C4` in `(20f)`: send
a generator to `5` or to `8`.  Inversion exchanges them.  This pair is again
an unbased `C2` torsor.  The order-168 lift is therefore an exact central-C4
shadow, not itself a C4 arm, and the two rows have not been embedded as the
two elements of one physical cross-Hom.  A physical endpoint-to-repair
cospan or semantic orientation is still needed to base it.

## 6. Broad whole-circle control

If one replaces `(2)` by `1` and retains only the digit and translated
`c1`-safe gate, the broad whole-circle-minus-`c1` rebuild has

```text
global content = 77,971,530 = 210*13^5,
profiles = 89,660,
units = 78,344,
nonunits = 11,316.                                       (21)
```

Its two local ghost slices happen to be unit, but `(21)` shows that neither
the inflated valuation nor the local success is the desired theorem.  This
broad grammar is retained only as a hostile/control; it is not a minimal
lawful repair and is not a dependency of the exact package below.

## 7. Holotopy interpretation and the remaining bridge

THM-2706 identifies the unique deterministic midpoint as the obstruction to
degree-two subdivision.  Equations `(3)--(7)` replace that missing object by
the coarsest local refinement which resolves the `c2` truth value.  This is
not an absolute smallest repair among every possible changed grammar.

Modulo null sets, multiplication by the two complementary Boolean
idempotents gives a Karoubi split of the middle integrand.  After coefficient
integration the labelled object is a pair `(Asafe,Aold)`; forgetting its arm
label is the augmentation

```text
Sigma: K direct_sum K -> K,       (a,b) |-> a+b,
K=F13[z]/(Phi_7).                                      (22)
```

The reverse twenty-profile anti-diagonal lies in `ker(Sigma)`.  Thus the
failure is not point-set descent--the Boolean split is an a.e. bijection--but
loss of a component selector under coefficient pushforward.  Unit extraction
does not commute with `(22)`.

The present calculation supplies neither that selector, a semantic owner
word, a common component phase, nor a live terminal current.  THM-2370's
deletion/no-reference boundary therefore remains relevant: changed-grammar
existence is not canonical current transport.  Absolute Frobenius preserves
each already formed forced safe unit, but preservation creates no physical
safe-arm cospan.

The nonoverlapping audited-request candidate
`THM-2720-unshifted-deepest-source-present-packet-global-disjointness.md`
verifies in its current proof-complete-candidate scope that every unchanged present packet is
disjoint from `E3`; deleting only its final unshifted `c3`-safe factor makes
exactly `78/91` label pairs positive.  The complementary finite-exact `c221`
reflection also has strict scalar `E3 -> Q_(3,{1,2})` cospans at the same
delayed phases.  In the relative `c3` truth split its safe arm is identically
zero and its danger arm equals the deleted packet, with all
`ell=1,...,6`, all thirteen shifts positive and minimum mass
`13828625/7792697484>1/564`; the symmetry is
`M(ell,s)=M(7-ell,-s)`.  These are present-factor, not `c2`, statements.  The
rebuild here retains the dynamic present factor and therefore is not an `E3`
semantic repair.  The clean nonoverlapping next universe is a present-free
`E3` private-row rebuild together with a relative-present truth-idempotent
sheet.

The cheapest next test is to hold one of the exact safe rows in `(18)` fixed
and expand an actual endpoint current through the two truth arms in `(6)`.
One must verify that a nonzero endpoint coefficient lands on `Asafe`, not
merely on their sum `Adel`.  A positive result would physicalize the repair
cospan; a cancellation reproducing the reverse twenty-profile hostile would
prove that an additional charged reference or owner-phase selector is
necessary.

## 8. Exact referee

Run

```text
python 04-computation/lrc14_minimal_c2_safe_ghost_transit_rebuild_thm2717.py
python -O 04-computation/lrc14_minimal_c2_safe_ghost_transit_rebuild_thm2717.py
```

Both runs must byte-match
`05-knowledge/results/lrc14_minimal_c2_safe_ghost_transit_rebuild_thm2717.out`.
The self-contained companion imports only the tracked THM-2698 exact
companion and verifies its LF hash
`45cc393a856c00342fdf84875a0bc5a6d4c3df196ab35bb9ac2aad3cfc966c25`
before use.  Its own LF hashes are
`81de8ea40fc00b05892c1430937d43dfb57112ed5b8a767282a37177ac4f94f2`
and
`bf7caaa73eaf2208cee4e27b0d09671bb34229458c99c7887004657001304437`.

The normal, optimized, and stored transcripts byte-match.  The companion
rebuilds both the deleted-gate and safe-repair banks on every rail, checks the
global contents and primitive determinant census, verifies `(7)` on the four
rails containing the two forced slices, checks all local cancellation data in
`(16)`, and replays the exact point and affine-split geometry in `(17)--(20)`.
